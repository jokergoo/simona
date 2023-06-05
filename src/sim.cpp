#include <Rcpp.h>
using namespace Rcpp;
#include <limits.h>


#include "transverse.h"
#include "dist.h"
#include "utils.h"
#include "term.h"

// ---------------------------------
// common ancestor 


void _assign_ancestor_max_value(S4 dag, NumericVector v, NumericMatrix& score, int i_ancestor, int id1, int id2) {
	score(id1, id2) = v[i_ancestor] > score(id1, id2)? v[i_ancestor] : score(id1, id2);
	score(id2, id1) = score(id1, id2);
}

void _assign_ancestor_id(S4 dag, NumericVector v, NumericMatrix& score, IntegerMatrix& id, IntegerMatrix& dd, LogicalVector l_all_ancestor, int i_ancestor, int id1, int id2, int offspring_1, int offspring_2) {
	if(id(id1, id2) == 0 || v[i_ancestor] > score(id1, id2)) {  // if id(id1, id2) is 0, it means it is not visited yet
		score(id1, id2) = v[i_ancestor];
		id(id1, id2) = i_ancestor + 1;
	} else if( std::abs(v[i_ancestor] - score(id1, id2)) < 1e-10 ) {  // if the v is the same, we compare the sumsq of distances to the max ancestor
		IntegerVector depth(v.size());
		if(dd(id1, id2) <= 0) {
			depth = cpp_dag_depth_bfs(dag, id(id1, id2), true, l_all_ancestor);
			dd(id1, id2) = depth[offspring_1 - 1] + depth[offspring_2 - 1];
		}
		depth = cpp_dag_depth_bfs(dag, i_ancestor+1, true, l_all_ancestor);
		int dist = depth[offspring_1-1] + depth[offspring_2 - 1];

		if(dist > dd(id1, id2)) {
			score(id1, id2) = v[i_ancestor];
			id(id1, id2) = i_ancestor + 1;
			dd(id1, id2) = dist;
		}
	}
	
	id(id2, id1) = id(id1, id2);
}

void _assign_ancestor_longest_distance(S4 dag, IntegerMatrix& dd, IntegerVector depth, int id1, int id2, int offspring_1, int offspring_2) {

	int dist = depth[offspring_1 - 1] + depth[offspring_2 - 1];

	if(dist > dd(id1, id2)) {
		dd(id1, id2) = dist;
		dd(id2, id1) = dd(id1, id2);
	}
}

void _assign_ancestor_shortest_distance(S4 dag, IntegerMatrix& dd, IntegerVector depth, int id1, int id2, int offspring_1, int offspring_2) {

	int dist = depth[offspring_1-1] + depth[offspring_2 - 1];

	if(dist < dd(id1, id2)) {
		dd(id1, id2) = dist;
		dd(id2, id1) = dd(id1, id2);
	}
}

// [[Rcpp::export]]
NumericMatrix cpp_max_ancestor_v(S4 dag, IntegerVector nodes, NumericVector v) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	NumericMatrix score(m, m);  // value of v

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	IntegerVector depth(n);
	double dist = 0;
	for(int k = 0; k < all_ancestor.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestor[k]-1, l_offspring, l_all_ancestor, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							_assign_ancestor_max_value(dag, v, score, k, id1, id2);
						}
					}
				}
			}
		}
	}

	return score;
}

// [[Rcpp::export]]
IntegerMatrix cpp_max_ancestor_id(S4 dag, IntegerVector nodes, NumericVector v) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	NumericMatrix score(m, m);  // value of v
	IntegerMatrix id(m, m);  // id, positive
	IntegerMatrix dd(m, m);  // sumsq dist to ancestor 

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		id(i, i) = nodes[i];
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	IntegerVector depth(n);
	double dist = 0;
	for(int k = 0; k < all_ancestor.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestor[k]-1, l_offspring, l_all_ancestor, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							_assign_ancestor_id(dag, v, score, id, dd, l_all_ancestor, k, id1, id2, offspring[i], offspring[j]);
						}
					}
				}
			}
		}
	}

	return id;
}

// [[Rcpp::export]]
IntegerMatrix cpp_distances(S4 dag, IntegerVector nodes, int type = 3) { // 3: longest distance or 4: shortest distance
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd(m, m);  // sumsq dist to ancestor 
	if(type == 4) {
		dd.fill(INT_MAX);
	}

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		dd(i, i) = 0;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	double dist = 0;
	for(int k = 0; k < all_ancestor.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestor[k]-1, l_offspring, l_all_ancestor, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			IntegerVector depth = cpp_dag_depth_bfs(dag, k+1, true, l_all_ancestor);

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							if(type == 3) {
								_assign_ancestor_longest_distance(dag, dd, depth, id1, id2, offspring[i], offspring[j]);
							} else if(type == 4) {
								_assign_ancestor_shortest_distance(dag, dd, depth, id1, id2, offspring[i], offspring[j]);
							}
						}
					}
				}
			}
		}
	}

	return dd;
}

// a leaf's most informative leaf is itself
// [[Rcpp::export]]
IntegerVector cpp_max_leaves_id(S4 dag, IntegerVector nodes, NumericVector v) {
	
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	int m = nodes.size();
	IntegerVector cl(m);
	LogicalVector is_leaf(n);
	for(int i = 0; i < m; i ++) {
		cl[i] = nodes[i];
		_find_connected_leaves(lt_children, nodes[i]-1, is_leaf);

		double max_v = 0;
		for(int j = 0; j < n; j ++) {
			if(is_leaf[j]) {
				if(v[j] > max_v) {
					max_v = v[j];
					cl[i] = j+1;
				}
			}
		}

		reset_logical_vector_to_false(is_leaf);
	}

	return cl;
}


// // aggregate to all ancestors
// // [[Rcpp::export]]
// NumericMatrix cpp_common_ancestor_aggregate_aic(S4 dag, IntegerVector nodes, NumericVector ic) {
// 	List lt_parents = dag.slot("lt_parents");

// 	int n = lt_parents.size();
// 	int m = nodes.size();
// 	NumericMatrix sim(m, m);
// 	sim.fill_diag(1.0);

// 	if(m <= 1) {
// 		return sim;
// 	}

// 	NumericVector v = 1.0/(1 + exp(-1.0/ic));

// 	LogicalVector ancestor1(n);
// 	LogicalVector ancestor2(n);

// 	for(int i = 0 ; i < m - 1; i ++) {

// 		_find_ancestor(lt_parents, nodes[i], ancestor1);
// 		ancestor1[nodes[i]-1] = true; // the node itself
// 		double sv1 = 0;
// 		for(int k = 0; k < n; k ++) {
// 			if(ancestor1[k]) {
// 				sv1 += v[k];
// 			}
// 		}

// 		for(int j = i + 1; j < m; j ++) {

// 			_find_ancestor(lt_parents, nodes[j], ancestor2);
// 			ancestor2[nodes[j]-1] = true; // the node itself
// 			double sv2 = 0;
// 			for(int k = 0; k < n; k ++) {
// 				if(ancestor2[k]) {
// 					sv2 += v[k];
// 				}
// 			}

// 			double sv3 = 0;
// 			for(int k = 0; k < n; k ++) {
// 				if(ancestor1[k] && ancestor2[k]) {
// 					sv3 += v[k];
// 				}
// 			}
// 			sim(i, j) = 2*sv3/sv1/sv2;
// 			sim(j, i) = sim(i, j);

// 			reset_logical_vector_to_false(ancestor2);
			
// 		}
// 		reset_logical_vector_to_false(ancestor1);
// 	}

// 	return sim;
// }

// // [[Rcpp::export]]
// NumericMatrix cpp_eps_XGraSM(S4 dag, IntegerVector nodes, NumericVector ic) {
// 	List lt_parents = dag.slot("lt_parents");

// 	int n = lt_parents.size();
// 	int m = nodes.size();
// 	NumericMatrix eps(m, m);
// 	eps.fill_diag(1);

// 	if(m <= 1) {
// 		return eps;
// 	}

// 	LogicalVector ancestor1(n);
// 	LogicalVector ancestor2(n);

// 	for(int i = 0 ; i < m - 1; i ++) {

// 		_find_ancestor(lt_parents, nodes[i], ancestor1);
// 		ancestor1[nodes[i]-1] = true; // the node itself

// 		for(int j = i + 1; j < m; j ++) {

// 			_find_ancestor(lt_parents, nodes[j], ancestor2);
// 			ancestor2[nodes[j]-1] = true; // the node itself

// 			int n_ca = 0;
// 			double ic_sum = 0;
// 			double ic_max = 0;

// 			for(int k = 0; k < n; k ++) {
// 				if(ancestor1[k] && ancestor2[k]) {
// 					if(ic[k] > ic_max) {
// 						ic_max = ic[k];
// 					}

// 					if(ic[k] > 0) {
// 						n_ca ++;
// 						ic_sum += ic[k];
// 					}
// 				}
// 			}

// 			eps(i, j) = ic_sum/n_ca/ic_max;
// 			eps(j, i) = eps(i, j);
// 		}
// 	}

// 	return eps;
// }

// // [[Rcpp::export]]
// NumericMatrix cpp_eps_EISI(S4 dag, IntegerVector nodes, NumericVector ic) {
// 	List lt_parents = dag.slot("lt_parents");
// 	List lt_children = dag.slot("lt_children");

// 	int n = lt_parents.size();
// 	int m = nodes.size();
// 	NumericMatrix eps(m, m);
// 	eps.fill_diag(1);

// 	if(m <= 1) {
// 		return eps;
// 	}

// 	LogicalVector ancestor1(n);
// 	LogicalVector ancestor2(n);

// 	for(int i = 0 ; i < m - 1; i ++) {

// 		_find_ancestor(lt_parents, nodes[i], ancestor1);
// 		ancestor1[nodes[i]-1] = true; // the node itself

// 		for(int j = i + 1; j < m; j ++) {

// 			_find_ancestor(lt_parents, nodes[j], ancestor2);
// 			ancestor2[nodes[j]-1] = true; // the node itself

// 			int n_ca = 0;
// 			double ic_sum = 0;
// 			double ic_max = 0;

// 			for(int k = 0; k < n; k ++) {
// 				if(ancestor1[k] && ancestor2[k]) {

// 					IntegerVector children = lt_children[k];
// 					for(int p = 0; p < children.size(); p ++) {
// 						if( (ancestor1[ children[p] - 1 ] && !ancestor2[ children[p] - 1 ]) || 
// 						    (!ancestor1[ children[p] - 1 ] && ancestor2[ children[p] - 1 ]) ) {
// 							n_ca ++;
// 							ic_sum += ic[k];
// 						}
// 					}

// 					if(ic[k] > ic_max) {
// 						ic_max = ic[k];
// 					}
// 				}
// 			}

// 			eps(i, j) = ic_sum/n_ca/ic_max;
// 			eps(j, i) = eps(i, j);
// 		}
// 	}

// 	return eps;
// }

// // [[Rcpp::export]]
// NumericMatrix cpp_sim_wang(S4 dag, IntegerVector nodes, NumericVector contribution) {
// 	List lt_parents = dag.slot("lt_parents");

// 	int n = lt_parents.size();
// 	int m = nodes.size();
// 	NumericMatrix sim(m, m);
// 	sim.fill_diag(1.0);

// 	if(m <= 1) {
// 		return sim;
// 	}

// 	NumericVector weight = -log(contribution);

// 	LogicalVector ancestor1(n, false);
// 	LogicalVector ancestor2(n, false);
		
// 	for(int i = 0; i < m - 1; i ++) {
// 		_find_ancestor(lt_parents, nodes[i], ancestor1);
// 		for(int j = i+1; j < m; j ++) {
// 			_find_ancestor(lt_parents, nodes[i], ancestor2);

// 			double x1;
// 			double x2;
// 			double x3;
// 			for(int k = 0; k < n; k ++) {
// 				if(ancestor1[k] && ancestor2[k]) {
// 					x1 += cpp_find_path_length_single(dag, k+1, i+1, weight, -1) +
// 					      cpp_find_path_length_single(dag, k+1, j+1, weight, -1);
// 				} else if(ancestor1[k]) {
// 					x2 += cpp_find_path_length_single(dag, k+1, i+1, weight, -1);
// 				} else if(ancestor2[k]) {
// 					x3 += cpp_find_path_length_single(dag, k+1, j+1, weight, -1);
// 				}
// 			}

// 			sim(i, j) = x1/(x1 + x2 + x3);
// 			sim(j, i) = sim(i, j);

// 			reset_logical_vector_to_false(ancestor2);
// 		}

// 		reset_logical_vector_to_false(ancestor1);
// 	}
// }

// // [[Rcpp::export]]
// NumericMatrix cpp_sim_shen(S4 dag, IntegerVector nodes, NumericVector ic) {
// 	int m = nodes.size();
// 	NumericMatrix sim(m, m);

// 	for(int i = 0; i < m; i ++) {
// 		sim(i, i) = 1 - atan(1/ic[i])/3.1415926*2;
// 	}

// 	if(m <= 1) {
// 		return sim;
// 	}

// 	NumericVector weight;

// 	IntegerMatrix mica_nodes = cpp_max_ancestor_id(dag, nodes, ic);

// 	for(int i = 0; i < m - 1; i ++) {
// 		for(int j = i+1; j < m; j ++) {
// 			IntegerVector path1 = cpp_find_path_single(dag, mica_nodes(i, j), j+1, weight, -1);
// 			IntegerVector path2 = cpp_find_path_single(dag, mica_nodes(i, j), j+1, weight, -1);

// 			double v = 0;
// 			for(int k = 0; k < path1.size(); k ++) {
// 				v += 1/ic[ path1[k]-1 ];
// 			}
// 			if(path2.size() > 1) {
// 				for(int k = 1; k < path2.size(); k ++) {
// 					v += 1/ic[ path2[k]-1 ];
// 				}
// 			}

// 			sim(i, j) = 1 - atan(v)/3.1415926*2;
// 			sim(j, i) = sim(i, j);
// 		}
// 	}

// 	return sim;
// }

// // [[Rcpp::export]]
// NumericMatrix cpp_sim_SSDD(S4 dag, IntegerVector nodes, NumericVector t) {
// 	int m = nodes.size();
// 	NumericMatrix sim(m, m);

// 	for(int i = 0; i < m; i ++) {
// 		sim(i, i) = 1 - atan(t[i])/3.1415926*2;
// 	}

// 	if(m <= 1) {
// 		return sim;
// 	}

// 	NumericVector weight;
// 	IntegerVector depth = _dag_depth(dag);
// 	NumericVector depth2 = as<NumericVector>(depth);

// 	IntegerMatrix lca_nodes = cpp_max_ancestor_id(dag, nodes, depth2);

// 	for(int i = 0; i < m - 1; i ++) {
// 		for(int j = i+1; j < m; j ++) {
// 			IntegerVector path1 = cpp_find_path_single(dag, lca_nodes(i, j), j+1, weight, -1);
// 			IntegerVector path2 = cpp_find_path_single(dag, lca_nodes(i, j), j+1, weight, -1);

// 			double v = 0;
// 			for(int k = 0; k < path1.size(); k ++) {
// 				v += t[ path1[k]-1 ];
// 			}
// 			if(path2.size() > 1) {
// 				for(int k = 1; k < path2.size(); k ++) {
// 					v += t[ path2[k]-1 ];
// 				}
// 			}

// 			sim(i, j) = 1 - atan(v)/3.1415926*2;
// 		}
// 	}

// 	return sim;
// }
