#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance &_instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), epsInt( 0.0 ), epsOpt( 0.0 )
{
	if( k == 0 ) k = instance.n_nodes;
}

// TODO
vector<Instance::Edge> edges;

void kMST_ILP::solve()
{
	try {

		// initialize CPLEX solver
		initCPLEX();

		// add common constraints
		modelCommon();
		// add model-specific constraints
		if( model_type == "scf" ) modelSCF();
		else if( model_type == "mcf" ) modelMCF();
		else if( model_type == "mtz" ) modelMTZ();
		else if( model_type == "dcc" ) cout << "DC-CUT model: no additional constraints\n";
		else if( model_type == "cec" ) cout << "CE-CUT model: no additional constraints\n";
		else {
			cerr << "No existing model chosen\n";
			exit( -1 );
		}

		// build model
		cplex = IloCplex( model );
		cplex.exportModel( "model.lp" );

		// set parameters
		cplex.setParam( IloCplex::Param::Threads, 1 ); // only use a single thread
		cplex.setParam( IloCplex::Param::TimeLimit, 3600 ); // set time limit to 1 hour
		cplex.setParam( IloCplex::Param::WorkMem, 8192 ); // set memory limit to 8 GB

		epsInt = cplex.getParam( IloCplex::Param::MIP::Tolerances::Integrality );
		epsOpt = cplex.getParam( IloCplex::Param::Simplex::Tolerances::Optimality );

		// set cut-callback for cycle-elimination cuts ("cec") or directed connection cuts ("dcc")
		// both for integer (mandatory!) and fractional (optional) solutions
		CutCallback cb( model_type, epsOpt, instance, x, z );
		if( model_type == "dcc" || model_type == "cec" ) {
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation;
			cplex.use( &cb, contextmask );
		}

		// solve model
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished.\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";

	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
		exit( -1 );
	}
}

// ----- private methods -----------------------------------------------

void kMST_ILP::initCPLEX()
{
	cout << "initialize CPLEX ... ";
	try {
		env = IloEnv();
		model = IloModel( env );
	}
	catch( IloException &e ) {
		cerr << "kMST_ILP: exception " << e.getMessage();
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
	}
	cout << "done.\n";
}



// Transform undirected edge vector to directed edge vector ({i,j} -> (i,j), (j,i)
static vector<Instance::Edge> createDirectedEdges(const vector<Instance::Edge> &uEdges)
{
	vector<Instance::Edge> dEdges;
	dEdges.resize(uEdges.size() * 2);

	auto it = copy(uEdges.cbegin(), uEdges.cend(), dEdges.begin());
	transform(uEdges.cbegin(),
			  uEdges.cend(),
			  it,
			  [](const Instance::Edge &e) {
					Instance::Edge ne = {e.v2, e.v1, e.weight}; return ne; 
			  });

	return dEdges;
}

// Create objective function to minimize sum of selected edge weights
static void createObjectiveFunction(IloEnv env, IloModel model, IloBoolVarArray xl, 
									vector<Instance::Edge> edges, u_int numEdges)
{
	IloExpr objFunction(env);
	for (u_int m = 0; m < numEdges; m++) {
		objFunction += xl[m] * edges[m].weight;
	}
	model.add(IloMinimize(env, objFunction));
	objFunction.end();
}

// Create variables x_ij to denote whether edge (i,j) is part of MST
static IloBoolVarArray createVarsX(IloEnv env, vector<Instance::Edge> edges, u_int numEdges)
{
	IloBoolVarArray xVarArray = IloBoolVarArray(env, numEdges);
	for (u_int k = 0; k < numEdges; k++) {
		const u_int i = edges[k].v1;
		const u_int j = edges[k].v2;
		xVarArray[k] = IloBoolVar(env, Tools::indicesToString("x", i, j).c_str());
	}
	return xVarArray;
}

// Create variables z_i to denote whether node i is part of the MST
static IloBoolVarArray createVarsV(IloEnv env, u_int num_nodes)
{
	IloBoolVarArray vVarArray = IloBoolVarArray(env, num_nodes);
	for (u_int i = 0; i < num_nodes; i++) {
		vVarArray[i] = IloBoolVar(env, Tools::indicesToString("v", i).c_str());
	}
	return vVarArray;
}

// Select exactly k nodes
static void createConstraint_selectKNodes(IloEnv env, IloModel model,
									  IloBoolVarArray z, Instance& instance, u_int k)
{
	IloExpr exprNumNodes(env);
	for (u_int i = 1; i < instance.n_nodes; i++) {
		exprNumNodes += z[i];
	}
	model.add(k == exprNumNodes);
	exprNumNodes.end();
}

// Select exactly (k-1) edges
static void createConstraint_selectedEdgeCount_equals_kMinus1(IloEnv env, IloModel model, IloBoolVarArray x,
											vector<Instance::Edge> edges, u_int numEdges, u_int k)
{
	IloExpr exprNumEdges(env);
	for (u_int m = 0; m < numEdges; m++) {
		const u_int i = edges[m].v1;
		const u_int j = edges[m].v2;
		if (i > 0 && j > 0) {
			exprNumEdges += x[m];
		}
	}
	model.add(exprNumEdges == k - 1);
	exprNumEdges.end();
}

// Artificial root node has exactly one selected outgoing edge
static void createConstraint_rootNodeHasOneOutgoingEdge(IloEnv env, IloModel model, IloBoolVarArray x,
														vector<Instance::Edge> edges, u_int numEdges)
{
	IloExpr exprRoot(env);
	for (u_int m = 0; m < numEdges; m++) {
		const u_int i = edges[m].v1;
		if (i == 0) {
			exprRoot += x[m];
		}
	}
	model.add(exprRoot == 1);
	exprRoot.end();
}

// Artificial root node has zero selected ingoing edges
static void createConstraint_rootNodeHasNoIncomingEdge(IloEnv env, IloModel model, IloBoolVarArray x,
													   vector<Instance::Edge> edges, u_int numEdges)
{
	IloExpr exprRoot(env);
	for (u_int m = 0; m < numEdges; m++) {
		const u_int j = edges[m].v2;
		if (j == 0) {
			exprRoot += x[m];
		}
	}
	model.add(exprRoot == 0);
	exprRoot.end();
}

// Selected nodes have at most (k-1), unselected nodes no selected outgoing edges
static void createConstraint_boundOutgoingEdges(IloModel model, IloBoolVarArray z,
												IloExprArray& exprOutDegree, Instance& instance, int k)
{
	for (u_int i = 0; i < instance.n_nodes; i++) {
		model.add(z[i] * (k - 1) >= exprOutDegree[i]);
	}
}

static void createConstraint_selectedNodeHasOutDegree_atLeastOne(IloModel model, IloBoolVarArray z,
																 IloExprArray& exprInDegree,
																 IloExprArray& exprOutDegree,
																 Instance& instance)
{
	for (u_int i = 0; i < instance.n_nodes; i++) {
		model.add(z[i] <= exprOutDegree[i] + exprInDegree[i]); 
	}
}


// Selected nodes have one, unselected nodes no incoming selected edges
static void createConstraint_selectedNodehasInDegreeOne_unselectedNone(IloModel model,
																	   IloBoolVarArray z, 
																	   IloExprArray& exprInDegree,
																	   Instance& instance)
{
	for (u_int i = 1; i < instance.n_nodes; i++) {
		model.add(exprInDegree[i] == z[i]);
	}
}


// Create InDegree expression array
static IloExprArray createInDegreeExprArray(IloEnv env, vector<Instance::Edge> edges,
											u_int numEdges, IloBoolVarArray x,  Instance& instance)
{
	IloExprArray exprInDegree(env, instance.n_nodes);
	for (u_int i = 0; i < instance.n_nodes; i++) {
		exprInDegree[i] = IloExpr(env);
	}
	for (u_int m = 0; m < numEdges; m++) {
		const u_int j = edges[m].v2;
		exprInDegree[j] += x[m];
	}
	return exprInDegree;
}

// Create OutDegree expression array
static IloExprArray createOutDegreeExprArray(IloEnv env, vector<Instance::Edge> edges,
											 u_int numEdges, IloBoolVarArray x, Instance& instance)
{
	IloExprArray exprOutDegree(env, instance.n_nodes);
	for (u_int i = 0; i < instance.n_nodes; i++) {
		exprOutDegree[i] = IloExpr(env);
	}
	for (u_int m = 0; m < numEdges; m++) {
		const u_int i = edges[m].v1;
		exprOutDegree[i] += x[m];
	}
	return exprOutDegree;
}


void kMST_ILP::modelCommon()
{
	try {

		edges = createDirectedEdges(instance.edges);
		const u_int numEdges = edges.size();


		// create x_ij variables
		this->x = createVarsX(env, edges, numEdges);

		// create v_i variables
		this->z = createVarsV(env, instance.n_nodes);

		// objective function
		createObjectiveFunction(env, model, this->x, edges, numEdges);

		// restrict selected edge count to k - 1
		createConstraint_selectedEdgeCount_equals_kMinus1(env, model, this->x, edges, numEdges, this->k);

		// Exactly one node is chosen as the tree root.
		createConstraint_rootNodeHasOneOutgoingEdge(env, model, this->x, edges, numEdges);
	
		// No edge is selected that goes back to artificial root node
		createConstraint_rootNodeHasNoIncomingEdge(env, model, this->x, edges, numEdges);

		IloExprArray edgesInDegrees = createInDegreeExprArray(env, edges, numEdges, this->x, instance);
		IloExprArray edgesOutDegrees = createOutDegreeExprArray(env, edges, numEdges, this->x, instance);

		// Unselected nodes have no outgoing selected edges, selected ones at most (k-1)
		createConstraint_boundOutgoingEdges(model, this->z, edgesOutDegrees, instance, this->k);

		// Selected nodes have at least one selected edge
		createConstraint_selectedNodeHasOutDegree_atLeastOne(model, this->z, edgesInDegrees, edgesOutDegrees, instance);
		
		// Selected nodes except for artificial root node have one, 
		// unselected ones have no incoming edge(s) selected
		createConstraint_selectedNodehasInDegreeOne_unselectedNone(model, this->z, edgesInDegrees, instance);
		
		// Select exactly k nodes
		createConstraint_selectKNodes(env, model, this->z, instance, this->k);
		edgesInDegrees.endElements();
		edgesOutDegrees.endElements();

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelCommon: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelCommon: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelSCF()
{
	try {

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build single commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelSCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelSCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMCF()
{
	try {

		// ++++++++++++++++++++++++++++++++++++++++++
		// TODO build multi commodity flow model
		// ++++++++++++++++++++++++++++++++++++++++++

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMCF: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMCF: unknown exception.\n";
		exit( -1 );
	}
}

void kMST_ILP::modelMTZ()
{
	try {

		
		// create d_i variables to store order
		this->d = IloIntVarArray(env, instance.n_nodes);
		for (u_int i = 0; i < instance.n_nodes; i++) {
			this->d[i] = IloIntVar(env, 0, k, Tools::indicesToString("u", i).c_str());
		}


		IloExpr exprLevelRootNode(env);
		exprLevelRootNode += this->d[0];
		// Artificial root node has level 0
		model.add(exprLevelRootNode == 0);
		exprLevelRootNode.end();
		
		u_int numEdges = instance.edges.size() * 2;

		for (u_int k = 0; k < numEdges; k++) {
			const u_int i = edges[k].v1;
			const u_int j = edges[k].v2;

			IloExpr exprHierarchy(env);
			exprHierarchy = this->d[i] + this->x[k] - this->d[j] - (-this->x[k] + 1) * this->k;
			// Set hierarchy 
			model.add(exprHierarchy <= 0);
			exprHierarchy.end();
		}

		for (u_int i = 0; i < instance.n_nodes; i++) {
			// Set order of unselected nodes to 0
			model.add(this->d[i] <= this->z[i] * (int) instance.n_nodes);
		}	

	}
	catch( IloException &e ) {
		cout << "kMST_ILP::modelMTZ: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cout << "kMST_ILP::modelMTZ: unknown exception.\n";
		exit( -1 );
	}
}

kMST_ILP::~kMST_ILP()
{
	// free CPLEX resources
//	x.end();
//	z.end();
//	if( model_type == "scf" || model_type == "mcf" ) f.end();
//	else if( model_type == "mtz" ) d.end();
	cplex.end();
	model.end();
	env.end();
}
