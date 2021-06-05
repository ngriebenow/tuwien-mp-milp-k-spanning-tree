#include "CutCallback.h"

CutCallback::CutCallback( string _cut_type, double _eps, Instance &_instance, IloBoolVarArray &_x, IloBoolVarArray &_z, vector<Instance::Edge> &_dEdges ) :
	cut_type( _cut_type ), eps( _eps ), instance( _instance ), context( NULL ), x( _x ), z( _z ), dEdges( _dEdges)
{
	arc_weights.resize( 2 * instance.n_edges );
}

CutCallback::~CutCallback()
{
}

void CutCallback::invoke( const IloCplex::Callback::Context &_context )
{
	context = &_context;
	env = context->getEnv();
	xsol = IloNumArray( env, 2 * instance.n_edges );
	zsol = IloNumArray( env, instance.n_nodes );

	// integer solution -> mandatory check for feasibility
	if( context->inCandidate() ) {
		if( !context->isCandidatePoint() ) throw IloCplex::Exception( -1, "Unbounded solution" );
		context->getCandidatePoint( x, xsol );
		context->getCandidatePoint( z, zsol );
		if( cut_type == "dcc" ) connectionCuts();
		else if( cut_type == "cec" ) cycleEliminationCuts();
	}
	// fractional solution -> optional to improve dual bounds
	else if( context->inRelaxation() ) {
		context->getRelaxationPoint( x, xsol );
		context->getRelaxationPoint( z, zsol );
		if( cut_type == "dcc" ) connectionCuts();
		else if( cut_type == "cec" ) cycleEliminationCuts();
	}
	else {
		xsol.end();
		zsol.end();
		throw IloCplex::Exception( -1, "Unexpected contextID" );
	}

	xsol.end();
	zsol.end();
}

/*
 * separation of directed connection cut inequalities
 */
void CutCallback::connectionCuts()
{
	

	try {

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// TODO find violated directed connection cut inequalities
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


			// add found violated cut to model
//			IloRange r( env, ... );
//			switch( context->getId() ) {
//				case IloCplex::Callback::Context::Id::Candidate:
//					context->rejectCandidate( r );
//					break;
//				case IloCplex::Callback::Context::Id::Relaxation:
//					context->addUserCut( r, IloCplex::UseCutForce, IloFalse );
//					break;
//				default:
//					r.end();
//					throw IloCplex::Exception( -1, "Unexpected contextID" );
//			}
//			r.end();

	}
	catch( IloException &e ) {
		cerr << "CutCallback: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "CutCallback: unknown exception.\n";
		exit( -1 );
	}
}

vector<int> indices;

void strongConnect()
{
	cout << "TODO";
}

void tarjan(const vector<int> &vertices, const vector<int> &edges)
{
	int index = 0;
	stack<int> stack;

	indices.clear();
	for (u_int i = 0; i < vertices.size(); i++)
	{
		indices.push_back(-1);
	}
	
	for (u_int i = 0; i < vertices.size(); i++)
	{
		if (indices[i] == -1)
		{
			strongConnect();
		}
	}
}

void CutCallback::testeig(const vector<int> &vertices, const vector<int> &edges,
					  vector<int> &flags)
{
	
	int real_root = 0;
	for (u_int i = 0; i < edges.size(); i++)
	{
		if (instance.edges[edges[i] / 2].v1 == 0) {
			real_root = instance.edges[edges[i]].v2;
			break;
		}
	}

	flags[real_root] = 1;
	dfs(real_root, vertices, edges, flags);

}

void CutCallback::dfs(const int v,
					  const vector<int> &vertices,
					  const vector<int> &edges,
					  vector<int> &flags)
{
	for (u_int i = 0; i < edges.size(); i++)
	{
		int v1 = dEdges[edges[i]].v1;
		int v2 = dEdges[edges[i]].v2;

		if (v1 == v)
		{
			int neighbor = v1 == v ? v2 : v1;

			flags[neighbor] += 1;

			if (flags[neighbor] == 1) {
				dfs(neighbor, vertices, edges, flags);
			}
		}
	}
	

/*
	for (u_int ei: instance.incidentEdges.at(v))
	{
		int v1 = instance.edges[ei].v1;
		int v2 = instance.edges[ei].v2;
		int neighbor = v1 == v ? v2 : v1;
		
		flags[neighbor] += 1;

		if (flags[neighbor] == 1) {
			dfs(neighbor, vertices, edges, flags);
		}
	}*/
}


/*
 * separation of cycle elimination cut inequalities
 */
void CutCallback::cycleEliminationCuts()
{
	try {
		if (context->inCandidate() == true) {
			const int value0 = xsol[0];
			cout << "test" << value0 << endl;

			vector<int> z_u;
			vector<int> x_u;
			vector<int> flags;
			
			for (u_int i = 0; i < zsol.getSize(); i++)
			{
				flags.push_back(0);
				
				if (zsol[i] == 1) {
					z_u.push_back(i);
				}
			}

			for (u_int i = 0; i < xsol.getSize(); i++)
			{
				if (xsol[i] == 1) {
					x_u.push_back(i);
				}
			}
			
			testeig(z_u, x_u, flags);

			for (u_int i = 0; i < flags.size(); i++)
			{
				if (flags[i] > 1)
				{
					cout << "rej";
				}
			}
			

		}

		//vector<IloNum> xu;
		//vector<IloNum> zu;

		



		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// TODO find violated cycle elimination cut inequalities
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//				// add found violated cut to model
//				IloRange r( env, ... );
//				switch( context->getId() ) {
//					case IloCplex::Callback::Context::Id::Candidate:
//						context->rejectCandidate( r );
//						break;
//					case IloCplex::Callback::Context::Id::Relaxation:
//						context->addUserCut( r, IloCplex::UseCutForce, IloFalse );
//						break;
//					default:
//						r.end();
//						throw IloCplex::Exception( -1, "Unexpected contextID" );
//				}
//				r.end();

	}
	catch( IloException &e ) {
		cerr << "CutCallback: exception " << e.getMessage();
		exit( -1 );
	}
	catch( ... ) {
		cerr << "CutCallback: unknown exception.\n";
		exit( -1 );
	}
}

/*
 * Dijkstra's algorithm to find a shortest path
 * Note: slow implementation with vectors in time O(n^2), instead of heaps etc.
 */
CutCallback::SPResultT CutCallback::shortestPath( u_int source, u_int target )
{
	vector<SPNodeT> nodes( instance.n_nodes );
	vector<bool> finished( instance.n_nodes, false ); // indicates finished nodes

	// initialization
	for( u_int v = 0; v < instance.n_nodes; v++ ) {
		nodes[v].pred = -1;
		nodes[v].pred_arc_id = -1;
		if( v == source ) nodes[v].weight = 0;
		else nodes[v].weight = numeric_limits<double>::max();
	}

	while( true ) {

		// find unfinished node with minimum weight to examine next
		// (should usually be done with heap or similar data structures)
		double wmin = numeric_limits<double>::max();
		u_int v = 0;
		for( u_int u = 0; u < instance.n_nodes; u++ ) {
			if( !finished[u] && nodes[u].weight < wmin ) {
				wmin = nodes[u].weight;
				v = u;
			}
		}

		// if all reachable nodes are finished or target node is reached -> stop
		if( wmin == numeric_limits<double>::max() || v == target ) break;

		// this node is finished now
		finished[v] = true;

		// update all adjacent nodes on outgoing arcs
		for( u_int e : instance.incidentEdges[v] ) {
			u_int a; // according arc id
			u_int u; // adjacent node
			if( instance.edges[e].v1 == v ) {
				a = e;
				u = instance.edges[e].v2;
			}
			else {
				a = e + instance.n_edges;
				u = instance.edges[e].v1;
			}
			// only examine adjacent node if unfinished
			if( !finished[u] ) {
				// check if weight at node u can be decreased
				if( nodes[u].weight > nodes[v].weight + arc_weights[a] ) {
					nodes[u].weight = nodes[v].weight + arc_weights[a];
					nodes[u].pred = v;
					nodes[u].pred_arc_id = a;
				}
			}
		}
	}

	SPResultT sp;
	sp.weight = 0;
	int v = target;
	while( v != (int) source && v != -1 ) {
		int a = nodes[v].pred_arc_id;
		if( a < 0 ) break;
		sp.weight += arc_weights[a];
		sp.path.push_back( a );
		v = nodes[v].pred;
	}
	return sp;
}

