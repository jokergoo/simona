#include <Rcpp.h>
using namespace Rcpp;

#include <limits.h>

#include "dist.h"
#include "transverse.h"
#include "utils.h"


// ----------------------------
//   depth and height


// [[Rcpp::export]]
IntegerVector cpp_dag_depth_bfs(S4 dag, int from_node = 0, bool use_max_dist = true, LogicalVector l_background = LogicalVector(0)) {
	List lt_children = dag.slot("lt_children");

	if(from_node == 0) {
		from_node = dag.slot("root");
	}
	int i_root = from_node - 1;

	int n = lt_children.size();
	IntegerVector d(n, -1);

	bool has_background = false;
	if(l_background.size() > 0) {
		has_background = true;
	}

	if(has_background) {
		if(!l_background[i_root]) {
			return d;
		}
	}

	LogicalVector current_nodes(n, false);
	current_nodes[i_root] = true;
	d[i_root] = 0;

	int n_current_nodes = 1;
	while(n_current_nodes) {
		for(int i = 0; i < n; i ++) {
			if(current_nodes[i]) {
				int cr = i;

				current_nodes[i] = false;
				
				IntegerVector children = lt_children[cr];

				if(children.size()) {
					for(int j = 0; j < children.size(); j ++) {
						int i_child = children[j]-1;

						if( (has_background && l_background[i_child]) || !has_background ) {
							if(d[i_child] == -1) {
								d[i_child] = d[cr] + 1;
							} else {
								if(use_max_dist) {
									if(d[cr] + 1 > d[i_child]) {
										d[i_child] = d[cr] + 1;
									}
								} else {
									if(d[cr] + 1 < d[i_child]) {
										d[i_child] = d[cr] + 1;
									}
								}
							}
							current_nodes[i_child] = true;
						}
					}
				}
			}
		}
		n_current_nodes = sum(current_nodes);
	}

	return d;
}

int _dag_depth_calc_d(int node, IntegerVector d, List lt_parents) {

	int i_node = node - 1;
	IntegerVector parents = lt_parents[i_node];
	
	if(parents.size() == 0) {
		d[i_node] = 0;
		return(0);
	}
	if(d[i_node] > -1) {
		return(d[i_node]);
	}

	IntegerVector v(parents.size());
	for(int i = 0; i < parents.size(); i ++) {
		v[i] = _dag_depth_calc_d(parents[i], d, lt_parents);
	}
	d[i_node] = max(v) + 1;

	return(d[i_node]);
}


// [[Rcpp::export]]
IntegerVector cpp_dag_depth_recursive(S4 dag) {
	List lt_parents = dag.slot("lt_parents");
	int nv = lt_parents.size();
	IntegerVector d(nv, -1);

	for(int node = 1; node <= nv; node ++) {
		_dag_depth_calc_d(node, d, lt_parents);
	}

	return(d);
}

////////////////////////////////////////////////////

// [[Rcpp::export]]
IntegerVector cpp_dag_height_bfs(S4 dag) {
	List lt_parents = dag.slot("lt_parents");
	IntegerVector leaves = dag.slot("leaves");

	int n = lt_parents.size();
	IntegerVector d(n, 0);

	LogicalVector current_nodes(n, false);
	for(int i = 0; i < leaves.size(); i ++) {
		int i_leaf = leaves[i] - 1;
		current_nodes[i_leaf] = true;
	}

	int n_current_nodes = sum(current_nodes);
	while(n_current_nodes) {
		for(int i = 0; i < n; i ++) {
			if(current_nodes[i]) {
				int cr = i;

				current_nodes[i] = false;
				
				IntegerVector parents = lt_parents[cr];

				if(parents.size()) {
					for(int j = 0; j < parents.size(); j ++) {
						int i_parent = parents[j]-1;

						if(d[cr] + 1 > d[i_parent]) {
							d[i_parent] = d[cr] + 1;
						}

						current_nodes[i_parent] = true;
					}
				}
			}
		}
		n_current_nodes = sum(current_nodes);
	}
	return d;
}

int _dag_height_calc_d(int node, IntegerVector d, List lt_children) {
	int i_node = node - 1;
	IntegerVector children = lt_children[i_node];
	
	if(children.size() == 0) {
		d[i_node] = 0;
		return(0);
	}
	if(d[i_node] > -1) {
		return(d[i_node]);
	}
	
	IntegerVector v(children.size());
	for(int i = 0; i < children.size(); i ++) {
		v[i] = _dag_height_calc_d(children[i], d, lt_children);
	}
	d[i_node] = max(v) + 1;

	return(d[i_node]);
}



// [[Rcpp::export]]
IntegerVector cpp_dag_height_recursive(S4 dag) {
	List lt_children = dag.slot("lt_children");

	int nv = lt_children.size();
	IntegerVector d(nv, -1);

	for(int node = 1; node <= nv; node ++) {
		_dag_height_calc_d(node, d, lt_children);
	}

	return(d);
}


// -------------------------
// number of connected leaves
// number of ancestor nodes
// number of offspring nodes


// [[Rcpp::export]]
IntegerVector cpp_n_ancestor(S4 dag, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");

	int n = lt_parents.size();
	IntegerVector num(n, 0);

	LogicalVector l_ancestor(n, false);
	for(int i = 0; i < n; i ++) {
		IntegerVector parents = lt_parents[i];
		if(parents.size() > 0) {
			
			_find_ancestor(lt_parents, i, l_ancestor, include_self);
			num[i] = sum(l_ancestor);

			reset_logical_vector_to_false(l_ancestor);
		}
	}
	return num;
}

// [[Rcpp::export]]
IntegerVector cpp_n_offspring(S4 dag, bool include_self = false) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < n; i ++) {
		IntegerVector children = lt_children[i];
		if(children.size() > 0) {
			
			_find_offspring(lt_children, i, l_offspring, include_self);
			num[i] = sum(l_offspring);

			reset_logical_vector_to_false(l_offspring);
		}
	}
	return num;
}

// [[Rcpp::export]]
IntegerVector cpp_n_leaves(S4 dag) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	LogicalVector l_leaf(n, false);
	for(int i = 0; i < n; i ++) {
		IntegerVector children = lt_children[i];
		if(children.size() > 0) {
			
			_find_connected_leaves(lt_children, i, l_leaf);
			num[i] = sum(l_leaf);

			reset_logical_vector_to_false(l_leaf);
		}
	}
	return num;
}


// -----------------------------
// cyclic node
void _go_child(List lt_children, int node, IntegerVector path, CharacterVector terms) {
	int i_node = node - 1;
	if(path.size() > 0) {
		for(int i = 0; i < path.size(); i ++) {
			if(path[i] == node) {
				String message("find a cyclic node:\n  [");
				for(int j = i; j < path.size(); j ++) {
					if(j == path.size()-1) {
						message = message + terms[ path[j]-1 ] + " " + terms[node-1] + "]";
					} else {
						message = message + terms[ path[j]-1 ] + " ";
					}
				}
				stop(message);
			}
		}
	}

	IntegerVector children = lt_children[i_node];
	for(int i = 0; i < children.size(); i ++) {
		IntegerVector path2 = clone(path);
		path2.push_back(node);
		_go_child(lt_children, children[i], path2, terms);
	}
}

// [[Rcpp::export]]
void cpp_check_cyclic_node(S4 dag) {
	List lt_children = dag.slot("lt_children");
	int root = dag.slot("root");
	CharacterVector terms = dag.slot("terms");

	IntegerVector path(0);
	_go_child(lt_children, root, path, terms);
}


// -------------------------
// IC


// [[Rcpp::export]]
NumericVector cpp_ic_meng(S4 dag, int correct) {
	List lt_children = dag.slot("lt_children");
	IntegerVector depth = _dag_depth(dag);
	int n_terms = dag.slot("n_terms");

	int max_depth = max(depth);
	int n = lt_children.size();

	NumericVector ic(n);
	for(int i = 0; i < n; i ++) {
		if(depth[i] == 0 || (correct == 0 && depth[i] == 1)) {
			ic[i] = 0;
		} else {
			LogicalVector l_offspring(n);
			_find_offspring(lt_children, i, l_offspring);

			double x = 0.0;
			for(int j = 0; j < n; j ++) {
				if(l_offspring[j]) {
					x = x + 1.0/depth[j];
				}
			}

			if(correct) {
				ic[i] = log(depth[i]+1)/log(max_depth+1)*(1 - log(x + 1)/log(n_terms));
			} else {
				ic[i] = log(depth[i])/log(max_depth)*(1 - log(x + 1)/log(n_terms));
			}
		}
	}

	return ic;
}


// [[Rcpp::export]]
IntegerVector cpp_n_annotations(S4 dag) {

	List lt_children = dag.slot("lt_children");
	List annotation = dag.slot("annotation");
	List lt_annotation = annotation["list"];
	CharacterVector anno_names = annotation["names"];
	int n_all_anno = anno_names.size();

	int n = lt_children.size();
	IntegerVector n_anno(n, 0);

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < n; i ++) {
		_find_offspring(lt_children, i, l_offspring, true);  //include self

		LogicalVector l_anno(n_all_anno, false);
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) {
				IntegerVector anno = lt_annotation[j];
				for(int k = 0; k < anno.size(); k ++) {
					l_anno[anno[k]-1] = true;
				}
			}
		}
		n_anno[i] = sum(l_anno);

		reset_logical_vector_to_false(l_offspring);
	}

	return n_anno;
}

// // [[Rcpp::export]]
// NumericVector cpp_ic_wang(S4 dag, NumericVector contribution) {

// 	List lt_parents = dag.slot("lt_parents");

// 	int n = lt_parents.size();
// 	NumericVector t(n);
// 	NumericVector weight = -log(contribution);

// 	for(int i = 0; i < n; i ++) {
// 		Rcout << i << "/" << n << "\n";
// 		LogicalVector ancestor(n, false);
// 		_find_ancestor(lt_parents, i+1, ancestor);

// 		for(int j = 0; j < n; j ++) {
// 			if(ancestor[j]) {
// 				t[i] += cpp_find_path_length_single(dag, j+1, i+1, weight, -1);
// 			}
// 		}

// 		reset_logical_vector_to_false(ancestor);
// 	}

// 	return 1/exp(t);
// }


