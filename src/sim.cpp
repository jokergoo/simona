#include <Rcpp.h>
using namespace Rcpp;


#include "transverse.h"
#include "utils.h"
#include "term.h"
#include "dist.h"
#include "common_ancestor.h"


// [[Rcpp::export]]
NumericMatrix cpp_sim_aic(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	NumericMatrix score(m, m);  // value of v

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
		score(i, i) = 1;
	}

	if(m <= 1) {
		return score;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);

	NumericVector sw = 1/(1+exp(-1/ic));
	NumericVector sv(n);

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
			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					sv[offspring[i] - 1] += sw[offspring[i] - 1];
				}
			}
		}
	}

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
							score(id1, id2) += 2*sw[all_ancestor[k]-1]/(sv[offspring[i]-1] + sv[offspring[j]-1]);
							score(id2, id1) = score(id1, id2);
						}
					}
				}
			}
		}
	}

	return score;
}




// [[Rcpp::export]]
NumericVector cpp_sim_wang(S4 dag, IntegerVector nodes, NumericVector contribution) {

	List lt_children = dag.slot("lt_children");
	List lt_children_relations = dag.slot("lt_children_relations");

	int n = lt_children.size();
	int m = nodes.size();
	NumericMatrix sim(m, m);
	sim.fill_diag(1);
	NumericVector ic(m);

	if(m <= 1) {
		return sim;
	}

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}
	
	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);

	for(int k = 0; k < all_ancestor.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestor[k]-1, l_offspring, l_all_ancestor, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff > 0) {
			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					LogicalVector l_ancestor(n);
					ic[id1] += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestor[k]-1, offspring[i]-1, l_all_ancestor);
				}
			}
		}

		if(noff > 1) {
			
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							sim(id1, id2) += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestor[k]-1, offspring[i]-1, l_all_ancestor) +
							                 _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestor[k]-1, offspring[j]-1, l_all_ancestor);
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {
			sim(i, j) = sim(i, j)/(ic[i] + ic[j]);
			sim(j, i) = sim(i, j);
		}
	}
	return sim;
}


void _assign_ancestor_max_wang_edge(S4 dag, NumericVector v, NumericMatrix& score, int i_ancestor, int id1, int id2) {
	score(id1, id2) = v[i_ancestor] > score(id1, id2)? v[i_ancestor] : score(id1, id2);
	score(id2, id1) = score(id1, id2);
}

// [[Rcpp::export]]
NumericMatrix cpp_sim_wang_edge(S4 dag, IntegerVector nodes) {

	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerVector global_depth = _dag_depth(dag);

	NumericMatrix sim(m, m);
	sim.fill_diag(1);
	NumericMatrix LCA_depth(m, m);
	LCA_depth.fill(-1);

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	if(m <= 1) {
		return sim;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	double sim_new = 0;
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
			IntegerVector depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestor[k], l_all_ancestor);
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							if(LCA_depth(id1, id2) < global_depth[all_ancestor[k]-1]) {
								LCA_depth(id1, id2) = global_depth[all_ancestor[k]-1];
								sim(id1, id2) = pow(global_depth[all_ancestor[k]-1], 2)/(global_depth[all_ancestor[k]-1] + depth[offspring[i]-1])/(global_depth[all_ancestor[k]-1] + depth[offspring[j]-1]);
								sim(id2, id1) = sim(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[all_ancestor[k]-1]) {
								if(sim(id1, id2) >= 0) {
									sim_new = pow(global_depth[all_ancestor[k]-1], 2)/(global_depth[all_ancestor[k]-1] + depth[offspring[i]-1])/(global_depth[all_ancestor[k]-1] + depth[offspring[j]-1]);
									if(sim_new < sim(id1, id2)) {
										sim(id1, id2) = sim_new;
										sim(id2, id1) = sim(id1, id2);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return sim;
}

// [[Rcpp::export]]
NumericMatrix cpp_sim_zhong(S4 dag, IntegerVector nodes, bool depth_via_LCA) {

	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	IntegerVector global_depth = _dag_depth(dag);

	NumericMatrix sim(m, m);
	sim.fill_diag(1);
	NumericMatrix LCA_depth(m, m);
	LCA_depth.fill(-1);

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	if(m <= 1) {
		return sim;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
	double sim_new = 0;
	int sigma_c;
	int sigma_a;
	int sigma_b;
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
			IntegerVector depth = cpp_dag_longest_dist_to_offspring(dag, all_ancestor[k], l_all_ancestor);
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							if(LCA_depth(id1, id2) < global_depth[all_ancestor[k]-1]) {
								LCA_depth(id1, id2) = global_depth[all_ancestor[k]-1];

								sigma_c = global_depth[all_ancestor[k]-1];
								if(depth_via_LCA) {
									sigma_a = sigma_c + depth[offspring[i] - 1];
									sigma_b = sigma_c + depth[offspring[j] - 1];
								} else {
									sigma_a = global_depth[offspring[i] - 1];
									sigma_b = global_depth[offspring[j] - 1];
								}
								sim(id1, id2) = 1 - 1.0/pow(2, sigma_c) - 0.5*(1.0/pow(2, sigma_a) + 1.0/pow(2, sigma_b));
								sim(id2, id1) = sim(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[all_ancestor[k]-1]) {
								if(sim(id1, id2) >= 0) {
									sigma_c = global_depth[all_ancestor[k]-1];
									if(depth_via_LCA) {
										sigma_a = sigma_c + depth[offspring[i] - 1];
										sigma_b = sigma_c + depth[offspring[j] - 1];
									} else {
										sigma_a = global_depth[offspring[i] - 1];
										sigma_b = global_depth[offspring[j] - 1];
									}
									sim_new = 1 - 1.0/pow(2, sigma_c) - 0.5*(1.0/pow(2, sigma_a) + 1.0/pow(2, sigma_b));
									if(sim_new < sim(id1, id2)) {
										sim(id1, id2) = sim_new;
										sim(id2, id1) = sim(id1, id2);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return sim;
}

const int PI = 3.1415926;


// [[Rcpp::export]]
NumericMatrix cpp_sim_shen(S4 dag, IntegerVector nodes, NumericVector ic) {
	int m = nodes.size();
	NumericMatrix sim(m, m);

	for(int i = 0; i < m; i ++) {
		sim(i, i) = 1 - atan(1/ic[ nodes[i]-1 ])/PI*2;
	}

	if(m <= 1) {
		return sim;
	}

	NumericVector weight;

	IntegerMatrix mica_nodes = cpp_max_ancestor_id(dag, nodes, ic);

	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {
			IntegerVector path1 = cpp_tpl_shortest_path(dag, mica_nodes(i, j), i+1);
			IntegerVector path2 = cpp_tpl_shortest_path(dag, mica_nodes(i, j), j+1);

			double v = 0;
			for(int k = 0; k < path1.size(); k ++) {
				v += 1/ic[ path1[k]-1 ];
			}
			if(path2.size() > 1) {
				for(int k = 1; k < path2.size(); k ++) {
					v += 1/ic[ path2[k]-1 ];
				}
			}

			sim(i, j) = 1 - atan(v)/PI*2;
			sim(j, i) = sim(i, j);
		}
	}

	return sim;
}


// [[Rcpp::export]]
NumericMatrix cpp_sim_SSDD(S4 dag, IntegerVector nodes, NumericVector t) {
	int m = nodes.size();
	NumericMatrix sim(m, m);

	for(int i = 0; i < m; i ++) {
		sim(i, i) = 1 - atan(t[i])/PI*2;
	}

	if(m <= 1) {
		return sim;
	}

	NumericVector weight;
	IntegerVector depth = _dag_depth(dag);
	NumericVector depth2 = as<NumericVector>(depth);

	IntegerMatrix lca_nodes = cpp_max_ancestor_id(dag, nodes, depth2);

	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {
			IntegerVector path1 = cpp_tpl_shortest_path(dag, lca_nodes(i, j), i+1);	
			IntegerVector path2 = cpp_tpl_shortest_path(dag, lca_nodes(i, j), j+1);

			double v = 0;
			for(int k = 0; k < path1.size(); k ++) {
				v += t[ path1[k]-1 ];
			}
			if(path2.size() > 1) {
				for(int k = 1; k < path2.size(); k ++) {
					v += t[ path2[k]-1 ];
				}
			}

			sim(i, j) = 1 - atan(v)/PI*2;
			sim(j, i) = sim(i, j);
		}
	}

	return sim;
}


// [[Rcpp::export]]
NumericMatrix cpp_common_ancestor_mean_IC_XGraSM(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	NumericMatrix score(m, m);
	IntegerVector na(m, m);

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	if(m <= 0) {
		return score;
	}

	IntegerVector all_ancestor = cpp_ancestor_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestor = integer_to_logical_vector(all_ancestor - 1, n);
	
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
			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							score(id1, id2) += ic[all_ancestor[k] - 1];
							if(ic[all_ancestor[k] - 1] > 0) {
								na(id1, id2) += 1;
							}
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < m; i ++) {
		for(int j = i; j < m; j ++) {
			score(i, j) = score(i, j)/na(i, j);
			score(j, i) = score(i, j);
		}
	}

	return score;
}

// [[Rcpp::export]]
NumericMatrix cpp_eps_EISI(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_parents = dag.slot("lt_parents");
	List lt_children = dag.slot("lt_children");

	int n = lt_parents.size();
	int m = nodes.size();
	NumericMatrix eps(m, m);
	eps.fill_diag(1);

	if(m <= 1) {
		return eps;
	}

	LogicalVector ancestor1(n);
	LogicalVector ancestor2(n);

	for(int i = 0 ; i < m - 1; i ++) {

		_find_ancestor(lt_parents, nodes[i], ancestor1);
		ancestor1[nodes[i]-1] = true; // the node itself

		for(int j = i + 1; j < m; j ++) {

			_find_ancestor(lt_parents, nodes[j], ancestor2);
			ancestor2[nodes[j]-1] = true; // the node itself

			int n_ca = 0;
			double ic_sum = 0;
			double ic_max = 0;

			for(int k = 0; k < n; k ++) {
				if(ancestor1[k] && ancestor2[k]) {

					IntegerVector children = lt_children[k];
					for(int p = 0; p < children.size(); p ++) {
						if( (ancestor1[ children[p] - 1 ] && !ancestor2[ children[p] - 1 ]) || 
						    (!ancestor1[ children[p] - 1 ] && ancestor2[ children[p] - 1 ]) ) {
							n_ca ++;
							ic_sum += ic[k];
						}
					}

					if(ic[k] > ic_max) {
						ic_max = ic[k];
					}
				}
			}

			eps(i, j) = ic_sum/n_ca/ic_max;
			eps(j, i) = eps(i, j);
		}
	}

	return eps;
}

