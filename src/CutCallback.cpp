#include "CutCallback.h"

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <cstdlib>

std::vector<int> myvector;

list<pair<u_int, u_int> > arcs;
double* weights;
Maxflow* algorithm;

CutCallback::CutCallback( string _cut_type, double _eps, Instance &_instance, IloBoolVarArray &_x, IloBoolVarArray &_z, vector<Instance::Edge> &_dEdges ) :
	cut_type( _cut_type ), eps( _eps ), instance( _instance ), context( NULL ), x( _x ), z( _z ), dEdges( _dEdges)
{
	arc_weights.resize( 2 * instance.n_edges );

	// initialise for random edge selection in cec
	for (u_int i=0; i < dEdges.size(); i++) myvector.push_back(i);

	// setup data structure for maxflow alg
	for (u_int i = 0; i < dEdges.size(); i++)
	{
		int v1 = dEdges[i].v1;
		int v2 = dEdges[i].v2;

		arcs.push_back( pair<u_int, u_int> (v1, v2));
	}

	weights = new double[dEdges.size()];

	int n = instance.n_nodes;
	int m = dEdges.size();

	// init algorithm
	Maxflow alg( n, m, arcs );
	algorithm = &alg;
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
	if (context->inCandidate() == true) {

		int n = instance.n_nodes;
		int m = dEdges.size();

		double* weights = new double[dEdges.size()];

		for (u_int i = 0; i < dEdges.size(); i++)
		{
			weights[i] = xsol[i];
		}

		// init algorithm
		Maxflow algorithm( n, m, arcs );

		int min_cut_cap =  numeric_limits<int>::max();
		int* min_cut;
		int min_node = -1;

		for (int i = 1; i < n; i++)
		{
			if (zsol[i] == 1) {

				// use algorithm
				algorithm.update( 0, i, weights );

				int* cut = new int[n];
				double f = algorithm.min_cut( 100, cut );

				if (f < min_cut_cap) {
					min_cut_cap = f;
					min_cut = cut;
					min_node = i;
					
					// already found?
					if (min_cut_cap < 1) {
						break;
					}
				}

			}
		}

		if (min_cut_cap < 1) {
			// violated inequality found


			IloExpr expr(env);

			expr += 1;
			expr -= z[min_node];

			for (u_int i = 0; i < dEdges.size(); i++)
			{
				int v1 = dEdges[i].v1;
				int v2 = dEdges[i].v2;

				if (min_cut[v1] != 2 && min_cut[v2] == 2)
				{
					expr += x[i];
				}
			}
			
			IloRange r = IloRange(env, 1, expr, std::numeric_limits<double>::infinity());
			context->rejectCandidate( r );
			expr.end();
			r.end();

		}

	}

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

/*
 * separation of cycle elimination cut inequalities
 */
void CutCallback::cycleEliminationCuts()
{

	try {
		if (context->inCandidate() == true) {

			// substitute arc weights
			for (u_int i = 0; i < xsol.getSize(); i++)
			{
				arc_weights[i] = 1 - xsol[i];
			}

			
  			std::random_shuffle ( myvector.begin(), myvector.end() );

			// search for a cycle
			for (u_int i = 0; i < dEdges.size(); i++)
			{
				int v1 = dEdges[myvector[i]].v1;
				int v2 = dEdges[myvector[i]].v2;

				SPResultT result = shortestPath(v2, v1);

				if (result.weight + arc_weights[myvector[i]] < 1.0) {
					// cycle detected!

					IloExpr expr(env);

					expr += x[myvector[i]];
					for (auto const& j : result.path)
					{
						expr += x[j];
					}
					IloRange r = IloRange(env, 0, expr, result.path.size());

					switch( context->getId() ) {
						case IloCplex::Callback::Context::Id::Candidate:
							context->rejectCandidate( r );
							break;
						case IloCplex::Callback::Context::Id::Relaxation:
							context->addUserCut( r, IloCplex::UseCutForce, IloFalse );
							break;
						default:
							r.end();
							throw IloCplex::Exception( -1, "Unexpected contextID" );
					}
					expr.end();
					r.end();

					break;
				}
			}
		}
		

		// switch( context->getId() ) {
		// 	case IloCplex::Callback::Context::Id::Candidate:
		// 		context->rejectCandidate( r );
		// 		break;
		// 	case IloCplex::Callback::Context::Id::Relaxation:
		// 		context->addUserCut( r, IloCplex::UseCutForce, IloFalse );
		// 		break;
		// 	default:
		// 		r.end();
		// 		throw IloCplex::Exception( -1, "Unexpected contextID" );
		// }
		// r.end();

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

