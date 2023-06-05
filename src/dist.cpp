#include <Rcpp.h>
using namespace Rcpp;

#include "transverse.h"
#include "utils.h"
#include "term.h"


// [[Rcpp::export]]
IntegerMatrix cpp_distance_new(S4 dag, IntegerVector nodes, bool use_max_dist = true) {

	int n = dag.slot("n_terms");
	int m = nodes.size();
	IntegerMatrix d(m, m);
	d = d - 1;

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	IntegerVector depth(n);
	for(int i = 0; i < m; i ++) {
		depth = cpp_dag_depth_bfs(dag, nodes[i], use_max_dist, l_all_ancestor); // distance to all offspring of nodes[i]
		for(int j = 0; j < n; j ++) {
			if(depth[j] >= 0 && nodes_ind[j] >= 0) {
				d(i, nodes_ind[j]) = depth[j];
				d(nodes_ind[j], i) = d(i, nodes_ind[j]);
			}
		}
	}

	return d;
}


// [[Rcpp::export]]
double cpp_find_path_length_single(S4 dag, int from, int to, NumericVector weight, int type) { // type 1: max, -1: min
	
	if(from == to) {
		return 0;
	}

	IntegerVector tpl_sorted = dag.slot("tpl_sorted");
	IntegerVector tpl_pos = dag.slot("tpl_pos"); 
	List lt_children = dag.slot("lt_children");
	List lt_children_relations = dag.slot("lt_children_relations");
	int n = lt_children.size();

	if(from > n || from < 1) {
		stop("'from' node is not in the DAG.");
	}
	if(to > n || to < 1) {
		stop("'to' node is not in the DAG.");
	}

	bool relation_is_defined = lt_children_relations.size() > 0 && weight.size() > 0;

	int i_from = from - 1;
	int i_to = to - 1;

	double default_dist;
	if(type == 1) {
		default_dist = R_NegInf;
	} else {
		default_dist = R_PosInf;
	}

	if(tpl_pos[i_from] > tpl_pos[i_to]) {
		return default_dist;
	}

	int i_pos_from = tpl_pos[i_from] - 1;
	int i_pos_to = tpl_pos[i_to] - 1;

	int len = i_pos_to - i_pos_from + 1;
	
	NumericVector dist(len, default_dist);
	dist[0] = 0; // from -> from

	for(int i_pos = i_pos_from; i_pos <= i_pos_to; i_pos ++) {

		int i_node = tpl_sorted[i_pos] - 1;  // the current node on the sorted list

		IntegerVector children = lt_children[i_node];
		if(relation_is_defined) {
			IntegerVector relations = lt_children_relations[i_node];
			for(int j = 0; j < children.size(); j ++) {
				int i_child = children[j] - 1;
				int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
				if(i_pos_child <= i_pos_to) {
					if(type == 1) { // max
						if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + weight[ relations[j] -1 ]) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + weight[ relations[j] -1 ];
						}
					} else { // min
						if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + weight[ relations[j] -1 ]) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + weight[ relations[j] -1 ];
						}
					}
					
				}
			}
		} else {
			for(int j = 0; j < children.size(); j ++) {
				int i_child = children[j] - 1;
				int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
				if(i_pos_child <= i_pos_to) {
					
					if(type == 1) { // max
						if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + 1) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
						}
					} else {
						if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + 1) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
						}
					}
					
				}
			}
		}
	}
	
	return dist[len-1];
}


// [[Rcpp::export]]
NumericMatrix cpp_find_path_length(S4 dag, IntegerVector nodes, NumericVector weight, int type) { // type 1: max, -1: min

	int n = nodes.size();
	NumericMatrix spl(n, n);
	spl.fill_diag(0);

	if(n == 0) {
		return spl;
	}

	if(n == 1) {
		spl(0, 0) = 0;
		return spl;
	}

	for(int i  = 0; i < n-1; i ++) {
		for(int j = i+1; j < n; j ++) {
			spl(i, j) = cpp_find_path_length_single(dag, nodes[i], nodes[j], weight, type);
			spl(j, i) = cpp_find_path_length_single(dag, nodes[j], nodes[i], weight, type);
		}
	}

	return spl;
}

// [[Rcpp::export]]
NumericMatrix cpp_find_path_length_2(S4 dag, IntegerVector from_nodes, IntegerVector to_nodes, NumericVector weight, int type) { // type 1: max, -1: min

	int n1 = from_nodes.size();
	int n2 = to_nodes.size();
	NumericMatrix spl(n1, n2);

	if(n1 == 0 || n2 == 0) {
		return spl;
	}

	for(int i  = 0; i < n1; i ++) {
		for(int j = 0; j < n2; j ++) {
			spl(i, j) = cpp_find_path_length_single(dag, from_nodes[i], to_nodes[j], weight, type);
		}
	}

	return spl;
}

// [[Rcpp::export]]
IntegerVector cpp_find_path_single(S4 dag, int from, int to, NumericVector weight, int type) { // type 1: max, -1: min
	
	if(from == to) {
		IntegerVector path(1, from);
		return path;
	}

	IntegerVector tpl_sorted = dag.slot("tpl_sorted");
	IntegerVector tpl_pos = dag.slot("tpl_pos"); 
	List lt_children = dag.slot("lt_children");
	List lt_children_relations = dag.slot("lt_children_relations");
	int n = lt_children.size();

	if(from > n || from < 1) {
		stop("'from' node is not in the DAG.");
	}
	if(to > n || to < 1) {
		stop("'to' node is not in the DAG.");
	}

	bool relation_is_defined = lt_children_relations.size() > 0 && weight.size() > 0;

	int i_from = from - 1;
	int i_to = to - 1;

	if(tpl_pos[i_from] > tpl_pos[i_to]) {
		IntegerVector path(0);
		return path;
	}

	int i_pos_from = tpl_pos[i_from] - 1;
	int i_pos_to = tpl_pos[i_to] - 1;

	int len = i_pos_to - i_pos_from + 1; // include both from and to
	double default_dist;
	if(type == 1) {
		default_dist = R_NegInf;
	} else {
		default_dist = R_PosInf;
	}
	NumericVector dist(len, default_dist);
	dist[0] = 0; // from -> from

	IntegerVector predecessor(len, -1);

	for(int i_pos = i_pos_from; i_pos <= i_pos_to; i_pos ++) {

		int i_node = tpl_sorted[i_pos] - 1;  // the current node on the sorted list

		IntegerVector children = lt_children[i_node];
		if(relation_is_defined) {
			IntegerVector relations = lt_children_relations[i_node];
			for(int j = 0; j < children.size(); j ++) {
				int i_child = children[j] - 1;
				int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
				if(i_pos_child <= i_pos_to) {
					if(type == 1) {
						if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + weight[ relations[j] -1 ]) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + weight[ relations[j] -1 ];
							predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
						}
					} else {
						if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + weight[ relations[j] -1 ]) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + weight[ relations[j] -1 ];
							predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
						}
					}
				}
			}
		} else {
			for(int j = 0; j < children.size(); j ++) {
				int i_child = children[j] - 1;
				int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
				if(i_pos_child <= i_pos_to) {
					
					if(type == 1) {
						if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + 1) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
							predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
						}
					} else {
						if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + 1) {
							dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
							predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
						}
					}
				}
			}
		}
	}

	IntegerVector path(len, -1);
	path[len-1] = to;

	int i = len - 1;
	int i_previous;
	while(1) {
		i_previous = predecessor[i];
		path[i-1] = tpl_sorted[i_previous + i_pos_from];
		i = i_previous;

		if(i_previous < 0) {
			break;
		}
	}

	path = path[path > 0];
	
	return path;
}

// taking no weight, shortest distance
// [[Rcpp::export]]
int cpp_distance_single(S4 dag, int from, int to) { // type 1: max, -1: min

	if(from == to) {
		return 0;
	}

	NumericVector weight;
	double dist = cpp_find_path_length_single(dag, from, to, weight, -1);

	int dist_int = int(dist);
	return dist_int;
}

// [[Rcpp::export]]
NumericMatrix cpp_distance(S4 dag, IntegerVector nodes) {
	NumericVector weight;
	return cpp_find_path_length(dag, nodes, weight, -1);
}

// [[Rcpp::export]]
NumericMatrix cpp_distance2(S4 dag, IntegerVector from_nodes, IntegerVector to_nodes) {
	NumericVector weight;
	return cpp_find_path_length_2(dag, from_nodes, to_nodes, weight, -1);
}

// treat the DAG as undirected, defined as min_{c in common_ancestor} d(c-a) + d(c-a)
// [[Rcpp::export]]
int cpp_distance_undirected_single(S4 dag, int v1, int v2) {
	NumericVector weight;

	if(v1 == v2) {
		return 0;
	}

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor1(n);
	LogicalVector l_ancestor2(n);

	_find_ancestor(lt_parents, v1-1, l_ancestor1);
	_find_ancestor(lt_parents, v2-1, l_ancestor2);

	int dist = R_PosInf;
	int dist2;
	for(int i = 0; i < n; i ++) {
		if(l_ancestor1[i] && l_ancestor2[i]) {
			dist2 = cpp_distance_single(dag, i+1, v1) + cpp_distance_single(dag, i+1, v2);
			if(dist2 < dist) {
				dist = dist2;
			}
		}
	}

	return dist;
}

// [[Rcpp::export]]
IntegerMatrix cpp_distance_undirected(S4 dag, IntegerVector nodes) {
	NumericVector weight;

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	int m = nodes.size();
	IntegerMatrix dist(m, m);

	if(m <= 1) {
		return dist;
	}

	LogicalVector l_ancestor1(n);
	LogicalVector l_ancestor2(n);

	for(int i = 0; i < m-1; i ++) {
		_find_ancestor(lt_parents, nodes[i], l_ancestor1);
		for(int j = i+1; j < m; j ++) {
			_find_ancestor(lt_parents, nodes[j], l_ancestor2);

			int d = R_PosInf;
			int d2;
			for(int k = 0; k < n; k ++) {
				if(l_ancestor1[k] && l_ancestor2[k]) {
					d2 = cpp_distance_single(dag, k+1, nodes[i]) + cpp_distance_single(dag, k+1, nodes[j]);
					if(d2 < d) {
						d = d2;
					}
				}
			}

			dist(i, j) = d;
			dist(j, i) = d;
		}
	}

	return dist;
}

// [[Rcpp::export]]
int cpp_longest_distance_single(S4 dag, int from, int to) { // type 1: max, -1: min

	if(from == to) {
		return 0;
	}

	NumericVector weight;
	double dist = cpp_find_path_length_single(dag, from, to, weight, 1);

	int dist_int = int(dist);
	return dist_int;
}

// [[Rcpp::export]]
NumericMatrix cpp_longest_distance(S4 dag, IntegerVector nodes) {
	NumericVector weight;
	return cpp_find_path_length(dag, nodes, weight, 1);
}


// [[Rcpp::export]]
IntegerMatrix cpp_longest_distance_undirected(S4 dag, IntegerVector nodes) {
	NumericVector weight;

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	int m = nodes.size();
	IntegerMatrix dist(m, m);

	if(m <= 1) {
		return dist;
	}

	LogicalVector l_ancestor1(n);
	LogicalVector l_ancestor2(n);

	for(int i = 0; i < m-1; i ++) {
		_find_ancestor(lt_parents, nodes[i]-1, l_ancestor1);
		for(int j = i+1; j < m; j ++) {
			_find_ancestor(lt_parents, nodes[j]-1, l_ancestor2);

			int d = R_PosInf;
			int d2;
			for(int k = 0; k < n; k ++) {
				if(l_ancestor1[k] && l_ancestor2[k]) {
					d2 = cpp_longest_distance_single(dag, k+1, nodes[i]) + cpp_longest_distance_single(dag, k+1, nodes[j]);
					if(d2 < d) {
						d = d2;
					}
				}
			}

			dist(i, j) = d;
			dist(j, i) = d;
		}
	}

	return dist;
}

// [[Rcpp::export]]
IntegerMatrix cpp_longest_distance_undirected_via_LCA(S4 dag, IntegerVector nodes, int which) { // which 1:a+b, 2:a, 3:b
	NumericVector weight;

	IntegerVector depth = _dag_depth(dag);

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	int m = nodes.size();
	IntegerMatrix dist(m, m);

	if(m <= 1) {
		return dist;
	}

	LogicalVector l_ancestor1(n);
	LogicalVector l_ancestor2(n);

	for(int i = 0; i < m-1; i ++) {
		_find_ancestor(lt_parents, nodes[i]-1, l_ancestor1);
		for(int j = i+1; j < m; j ++) {
			_find_ancestor(lt_parents, nodes[j]-1, l_ancestor2);

			int d = 0;
			int v = 0;
			for(int k = 0; k < n; k ++) {
				if(l_ancestor1[k] && l_ancestor2[k]) {
					if(depth[k] > d) {
						d = depth[k];
						if(which == 1) {
							v = cpp_find_path_length_single(dag, k+1, nodes[i], weight, 1) + 
							    cpp_find_path_length_single(dag, k+1, nodes[j], weight, 1);
						} else if(which == 2) {
							v = cpp_find_path_length_single(dag, k+1, nodes[i], weight, 1);
						} else {
							v = cpp_find_path_length_single(dag, k+1, nodes[j], weight, 1);
						}
					}
				}
			}

			dist(i, j) = v;
			dist(j, i) = v;
		}
	}

	return dist;
}

