#include <Rcpp.h>
using namespace Rcpp;
#include <limits.h>


#include "transverse.h"
#include "utils.h"

// ---------------------------------
// All the functions have a similar algorithm. The aim is to pick the common ancestor which maximize a certain metric,
// let's denote is as s_c(a, b), where c is the common ancestor of a and b. We want to pick the c which has the maximal s.
//
//   for a group of terms, first everything is restricted in the union of all ancestors of these temrs (include themselves).
//   Then for each of the ancestor k, we look for all its offsprings. For each pair of the offspring a and b, we can have
//   a score filled for k, s_k(a, b). When we go through other ancestor k', we can update s_k to s_k'(a, b), if the latter value
//   is higher.
//
// - for LCA, the score is the depth
// - for MICA, the score is the IC
// - for distance, the score is dc_a + dc_b, which is the sum of distance from c to a and distance from c to b
//
void _assign_ancestor_max_value(S4 dag, NumericVector v, NumericMatrix& score, int i_ancestor, int id1, int id2) {
	score(id1, id2) = v[i_ancestor] > score(id1, id2)? v[i_ancestor] : score(id1, id2);
	score(id2, id1) = score(id1, id2);
}

void _assign_ancestor_id(S4 dag, NumericVector v, NumericMatrix& score, IntegerMatrix& id, IntegerMatrix& dd, IntegerVector depth, int i_ancestor, int id1, int id2, int offspring_1, int offspring_2) {
	int dist = depth[offspring_1-1] + depth[offspring_2 - 1];
	if(id(id1, id2) == 0 || v[i_ancestor] > score(id1, id2)) {  // if id(id1, id2) is 0, it means it is not visited yet
		score(id1, id2) = v[i_ancestor];
		id(id1, id2) = i_ancestor + 1;
		dd(id1, id2) = dist;
	} else if( std::abs(v[i_ancestor] - score(id1, id2)) < 1e-10 ) {  // if the v is the same, we compare the sumsq of distances to the max ancestor

		if(dist < dd(id1, id2)) {
			score(id1, id2) = v[i_ancestor];
			id(id1, id2) = i_ancestor + 1;
			dd(id1, id2) = dist;
		}
	}
	
	id(id2, id1) = id(id1, id2);
}

void _assign_ancestor_longest_distance(S4 dag, IntegerMatrix& dd, IntegerVector depth, int id1, int id2, int offspring_1, int offspring_2) {

	int dist = depth[offspring_1 - 1] + depth[offspring_2 - 1];

	if(dd(id1, id2) == -1 || dist > dd(id1, id2)) {
		dd(id1, id2) = dist;
		dd(id2, id1) = dd(id1, id2);
	}
}

void _assign_ancestor_shortest_distance(S4 dag, IntegerMatrix& dd, IntegerVector shortest_dist_to_offspring, int id1, int id2, int offspring_1, int offspring_2) {

	int dist = shortest_dist_to_offspring[offspring_1 - 1] + shortest_dist_to_offspring[offspring_2 - 1];

	if(dd(id1, id2) == -1 || dist < dd(id1, id2)) {
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
		score(i, i) = v[nodes[i]-1];
	}

	if(m <= 1) {
		return score;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

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
							_assign_ancestor_max_value(dag, v, score, all_ancestors[k]-1, id1, id2);
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
	IntegerMatrix dd(m, m); 

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		id(i, i) = nodes[i];
	}

	if(m <= 1) {
		return id;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {

			IntegerVector depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							_assign_ancestor_id(dag, v, score, id, dd, depth, all_ancestors[k]-1, id1, id2, offspring[i], offspring[j]);
						}
					}
				}
			}
		}
	}

	return id;
}

const int USE_LONGEST_DISTANCE = 1;
const int USE_SHORTEST_DISTANCE = 2;

IntegerMatrix cpp_distances(S4 dag, IntegerVector nodes, int type = 1) { // 1: longest distance or 0: shortest distance
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd(m, m);  // sumsq dist to ancestor 
	dd.fill(-1);
	
	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		dd(i, i) = 0;
	}

	if(m <= 1) {
		return dd;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			IntegerVector depth(n);
			IntegerVector shortest_dist_to_offspring(n);
			if(type == USE_LONGEST_DISTANCE) {
				depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);
			} else if(type == USE_SHORTEST_DISTANCE) {
				shortest_dist_to_offspring = cpp_dag_shortest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);
			}

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							if(type == USE_LONGEST_DISTANCE) {
								_assign_ancestor_longest_distance(dag, dd, depth, id1, id2, offspring[i], offspring[j]);
							} else if(type == USE_SHORTEST_DISTANCE) {
								_assign_ancestor_shortest_distance(dag, dd, shortest_dist_to_offspring, id1, id2, offspring[i], offspring[j]);
							}
						}
					}
				}
			}
		}
	}

	return dd;
}

IntegerMatrix cpp_longest_distances_via_CA(S4 dag, IntegerVector nodes) {
	return cpp_distances(dag, nodes, USE_LONGEST_DISTANCE);
}

// [[Rcpp::export]]
IntegerMatrix cpp_shortest_distances_via_CA(S4 dag, IntegerVector nodes) {
	return cpp_distances(dag, nodes, USE_SHORTEST_DISTANCE);
}


// [[Rcpp::export]]
IntegerMatrix cpp_longest_distances_via_LCA(S4 dag, IntegerVector nodes) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd(m, m);
	dd.fill(-1);
	IntegerMatrix LCA_depth(m, m);
	LCA_depth.fill(-1);
	
	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		dd(i, i) = 0;
	}

	if(m <= 1) {
		return dd;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

	IntegerVector global_depth = cpp_dag_depth(dag);
	
	double dist = 0;
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			IntegerVector depth(n);
			depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
						
						if(id2 >= 0) {
							if(LCA_depth(id1, id2) < global_depth[ all_ancestors[k]-1 ]) {
								LCA_depth(id1, id2) = global_depth[ all_ancestors[k]-1 ];
								dist = depth[offspring[i] - 1] + depth[offspring[j] - 1];
								dd(id1, id2) = dist;
								dd(id2, id1) = dd(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[ all_ancestors[k]-1 ]) {
								dist = depth[offspring[i] - 1] + depth[offspring[j] - 1];
								if(dd(id1, id2) == -1 || dist > dd(id1, id2)) {
									dd(id1, id2) = dist;
									dd(id2, id1) = dd(id1, id2);
								}
							}

						}
					}
				}
			}
		}
	}

	return dd;
}

// [[Rcpp::export]]
List cpp_longest_distances_from_LCA(S4 dag, IntegerVector nodes) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd1(m, m);
	dd1.fill(-1);
	IntegerMatrix dd2(m, m);
	dd2.fill(-1);
	IntegerMatrix LCA_depth(m, m);
	LCA_depth.fill(-1);
	
	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		dd1(i, i) = 0;
		dd2(i, i) = 0;
	}

	if(m <= 1) {
		List lt = List::create(Named("left") = dd1 , Named("right") = dd2);
		return lt;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

	IntegerVector global_depth = cpp_dag_depth(dag);
	
	double dist = 0;
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			IntegerVector depth(n);
			depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
						
						if(id2 >= 0) {
							if(LCA_depth(id1, id2) < global_depth[ all_ancestors[k]-1 ]) {
								LCA_depth(id1, id2) = global_depth[ all_ancestors[k]-1 ];

								dd1(id1, id2) = depth[offspring[i] - 1];
								dd1(id2, id1) = dd1(id1, id2);
								dd2(id1, id2) = depth[offspring[j] - 1];
								dd2(id2, id1) = dd2(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[ all_ancestors[k]-1 ]) {
								dist = depth[offspring[i] - 1] + depth[offspring[j] - 1];
								if(dd1(id1, id2) == -1 || dist > dd1(id1, id2) + dd2(id1, id2)) {
									dd1(id1, id2) = depth[offspring[i] - 1];
									dd1(id2, id1) = dd1(id1, id2);
									dd2(id1, id2) = depth[offspring[j] - 1];
									dd2(id2, id1) = dd2(id1, id2);
								}
							}

						}
					}
				}
			}
		}
	}

	List lt = List::create(Named("left") = dd1 , Named("right") = dd2);

	return lt;
}

IntegerMatrix cpp_distances_directed(S4 dag, IntegerVector nodes, int type = 1) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd(m, m);  // sumsq dist to ancestor 
	dd.fill(-1);

	if(m <= 1) {
		return dd;
	}
	
	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	for(int k = 0; k < m; k ++) {
		_find_offspring_within_background(lt_children, nodes[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		IntegerVector depth(n);
		if(type == USE_LONGEST_DISTANCE) {
			depth = cpp_dag_longest_dist_to_offspring(dag, nodes[k], l_all_ancestors);
		} else if(type == USE_SHORTEST_DISTANCE) {
			depth = cpp_dag_shortest_dist_to_offspring(dag, nodes[k], l_all_ancestors);
		}

		int id_p = nodes_ind[ nodes[k]-1 ];

		for(int i = 0; i < noff; i ++) {
			int id = nodes_ind[ offspring[i]-1 ];
			if(id >= 0) {
				dd(id_p, id) = depth[offspring[i]-1];
			}
		}
	}

	return dd;
}

// [[Rcpp::export]]
IntegerMatrix cpp_longest_distances_directed(S4 dag, IntegerVector nodes) {
	return cpp_distances_directed(dag, nodes, USE_LONGEST_DISTANCE);
}

// [[Rcpp::export]]
IntegerMatrix cpp_shortest_distances_directed(S4 dag, IntegerVector nodes) {
	return cpp_distances_directed(dag, nodes, USE_SHORTEST_DISTANCE);
}


// [[Rcpp::export]]
IntegerMatrix cpp_nearest_common_ancestor(S4 dag, IntegerVector nodes) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix dd(m, m); 
	dd.fill(INT_MAX - 1);
	IntegerMatrix id(m, m);

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		dd(i, i) = 0;
		id(i, i) = nodes[i];
	}

	if(m <= 1) {
		return id;
	}

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

	IntegerVector global_depth = _dag_depth(dag);
	
	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 1) {
			IntegerVector shortest_dist_to_offspring = cpp_dag_shortest_dist_to_offspring(dag, all_ancestors[k], l_all_ancestors);

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							int dist = shortest_dist_to_offspring[offspring[i] - 1] + shortest_dist_to_offspring[offspring[j] - 1];
							if(dd(id1, id2) > dist) {
								dd(id1, id2) = dist;
								id(id1, id2) = all_ancestors[k];
								id(id2, id1) = id(id1, id2);
							} else if(dd(id1, id2) == dist) {
								if(global_depth[all_ancestors[k] - 1] > global_depth[ id(id1, id2)-1 ]) {
									id(id1, id2) = all_ancestors[k];
									id(id2, id1) = id(id1, id2);
								}
							}
						}
					}
				}
			}
		}
	}

	return id;
}
