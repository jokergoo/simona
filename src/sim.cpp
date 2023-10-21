#include <Rcpp.h>
using namespace Rcpp;


#include "traverse.h"
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

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

	NumericVector sw = 1/(1+exp(-1/ic));
	NumericVector sv(n);

	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff >= 1) {
			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					sv[offspring[i] - 1] += sw[all_ancestors[k]-1];
				}
			}
		}
	}

	for(int k = 0; k < all_ancestors.size(); k ++) {
		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

		IntegerVector offspring = _which(l_offspring);
		reset_logical_vector_to_false(l_offspring);

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		offspring = offspring + 1;

		if(noff >= 1) {
			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							score(id1, id2) += 2*sw[all_ancestors[k]-1]/(sv[offspring[i]-1] + sv[offspring[j]-1]);
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
NumericVector cpp_sim_wang(S4 dag, IntegerVector nodes, NumericVector contribution, bool correct = false) {

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
	
	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

	double c = 0;
	if(correct) {
		c = max(contribution)/(1 - max(contribution));
	}

	for(int k = 0; k < all_ancestors.size(); k ++) {

		if(k % 100 == 0) {
			message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
			message("going through " + std::to_string(k) + " / " + std::to_string(all_ancestors.size()) + " ancestors ...", false);
		}

		_find_offspring_within_background(lt_children, all_ancestors[k]-1, l_offspring, l_all_ancestors, true);

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
					LogicalVector l_ancestors(n);
					ic[id1] += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[i]-1, l_all_ancestors, correct, c);
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
							sim(id1, id2) += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[i]-1, l_all_ancestors, correct, c) +
							                 _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[j]-1, l_all_ancestors, correct, c);
						}
					}
				}
			}
		}
	}

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(all_ancestors.size()) + " / " + std::to_string(all_ancestors.size()) + " ancestors ... Done.", true);

	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {
			sim(i, j) = sim(i, j)/(ic[i] + ic[j]);
			sim(j, i) = sim(i, j);
		}
	}
	return sim;
}

// [[Rcpp::export]]
NumericMatrix cpp_wang_sv_to_sim(NumericMatrix sv) {
	int na = sv.nrow();
	int n = sv.ncol();

	NumericVector ic(n);

	for(int i = 0; i < n; i ++) {
		for(int j = 0; j < na; j ++) {
			ic[i] += sv(j, i);
		}
	}

	NumericMatrix sim(n, n);
	sim.fill_diag(1);

	if(n <= 1) {
		return sim;
	}

	for(int i = 0; i < n - 1; i ++) {
		for(int j = i + 1; j < n; j ++) {
			for(int k = 0; k < na; k ++) {
				if(std::abs(sv(k, i)) > 1e-10 && std::abs(sv(k, j)) > 1e-10) {
					sim(i, j) += sv(k, i) + sv(k, j);
				}
			}
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

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	double sim_new = 0;
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
							if(LCA_depth(id1, id2) < global_depth[all_ancestors[k]-1]) {
								LCA_depth(id1, id2) = global_depth[all_ancestors[k]-1];
								if(std::abs(global_depth[all_ancestors[k]-1]) < 1e-10) {
									sim(id1, id2) = 0;
								} else {
									sim(id1, id2) = pow(global_depth[all_ancestors[k]-1], 2)/(global_depth[all_ancestors[k]-1] + depth[offspring[i]-1])/(global_depth[all_ancestors[k]-1] + depth[offspring[j]-1]);
								}
								sim(id2, id1) = sim(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[all_ancestors[k]-1]) {
								if(sim(id1, id2) >= 0) {
									if(std::abs(global_depth[all_ancestors[k]-1]) < 1e-10) {
										sim_new = 0;
									} else {
										sim_new = pow(global_depth[all_ancestors[k]-1], 2)/(global_depth[all_ancestors[k]-1] + depth[offspring[i]-1])/(global_depth[all_ancestors[k]-1] + depth[offspring[j]-1]);
									}
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

	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);
	
	double sim_new = 0;
	int sigma_c;
	int sigma_a;
	int sigma_b;
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
							if(LCA_depth(id1, id2) < global_depth[all_ancestors[k]-1]) {
								LCA_depth(id1, id2) = global_depth[all_ancestors[k]-1];

								sigma_c = global_depth[all_ancestors[k]-1];
								sigma_a = depth[offspring[i] - 1]; // dist from lca to a
								sigma_b = depth[offspring[j] - 1]; // dist from lca to b
								
								sim(id1, id2) = 1 - 1.0/pow(2, sigma_c)*(1- 1.0/pow(2, sigma_a+1) - 1.0/pow(2, sigma_b+1));
								sim(id2, id1) = sim(id1, id2);
							} else if(LCA_depth(id1, id2) == global_depth[all_ancestors[k]-1]) {
								if(sim(id1, id2) >= 0) {
									sigma_c = global_depth[all_ancestors[k]-1];
									sigma_a = depth[offspring[i] - 1];
									sigma_b = depth[offspring[j] - 1];
									
									sim_new = 1 - 1.0/pow(2, sigma_c)*(1- 1.0/pow(2, sigma_a+1) - 1.0/pow(2, sigma_b+1));
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
NumericMatrix cpp_common_ancestor_mean_IC_XGraSM(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_children = dag.slot("lt_children");
	
	int n = lt_children.size();
	int m = nodes.size();

	NumericMatrix score(m, m);
	IntegerMatrix na(m, m); // number of ancestors

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	if(m <= 0) {
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
			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					for(int j = i; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
							
						if(id2 >= 0) {
							score(id1, id2) += ic[all_ancestors[k] - 1];
							if(ic[all_ancestors[k] - 1] > 0) {
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
			if(na(i, j) == 0) {
				score(i, j) = 0;
			} else {
				score(i, j) = score(i, j)/na(i, j);
			}
			score(j, i) = score(i, j);
		}
	}

	return score;
}

// [[Rcpp::export]]
NumericMatrix cpp_common_ancestor_mean_IC_EISI(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_parents = dag.slot("lt_parents");
	List lt_children = dag.slot("lt_children");

	int n = lt_parents.size();
	int m = nodes.size();
	NumericMatrix mean_ic(m, m);
	for(int i = 0; i < m; i ++) {
		mean_ic(i, i) = ic[nodes[i] - 1];
	}

	if(m <= 1) {
		return mean_ic;
	}

	LogicalMatrix m_ancestors(m, n);
	LogicalVector la(n);
	for(int i = 0 ; i < m; i ++) {
		_find_ancestors(lt_parents, nodes[i] - 1, la, true);
		m_ancestors(i, _) = la;
		reset_logical_vector_to_false(la);
	}

	for(int i = 0 ; i < m - 1; i ++) {
		for(int j = i + 1; j < m; j ++) {

			int n_ca = 0;
			double ic_sum = 0;

			// here looks for EICA (exclusively inherited and all common ancestors)
			for(int k = 0; k < n; k ++) {
				if(m_ancestors(i, k) && m_ancestors(j, k)) { // if common ancestor

					// check C_h(a) \cap ( A-B union B-A) != empty
					bool flag = false;
					IntegerVector children = lt_children[k];
					for(int p = 0; p < children.size(); p ++) {
						if( (m_ancestors(i, children[p] - 1 ) && !m_ancestors(j, children[p] - 1 )) || 
						    (!m_ancestors(i, children[p] - 1 ) && m_ancestors(j, children[p] - 1 )) ) {
							
							flag = true;
							break;
						}
					}

					if(flag) {
						n_ca ++;
						ic_sum += ic[k];
					}
				}
			}

			mean_ic(i, j) = ic_sum/n_ca;
			mean_ic(j, i) = mean_ic(i, j);
		}
	}

	return mean_ic;
}

// [[Rcpp::export]]
NumericMatrix cpp_sim_ancestor(S4 dag, IntegerVector nodes) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	int m = nodes.size();
	NumericMatrix sim(m, m);
	NumericVector v_occur(m);
	NumericMatrix m_intersect(m, m);

	IntegerVector nodes_ind(n, -1);  // mapping between n and m indices
	for(int i = 0; i < m; i ++) {
		nodes_ind[ nodes[i]-1 ] = i;
	}

	if(m <= 1) {
		if(m == 1) {
			sim(0, 0) = 1;
		}
		return sim;
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

		if(noff >= 1) {

			for(int i = 0; i < noff; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				if(id1 >= 0) {
					v_occur[id1] = v_occur[id1] + 1;
				}
			}

			for(int i = 0; i < noff - 1; i ++) {
				int id1 = nodes_ind[ offspring[i]-1 ];
				
				if(id1 >= 0) {

					for(int j = i+1; j < noff; j ++) {
						int id2 = nodes_ind[ offspring[j]-1 ];
						
						if(id2 >= 0) {	
							m_intersect(id1, id2) = m_intersect(id1, id2) + 1;
							m_intersect(id2, id1) = m_intersect(id1, id2);
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < m; i ++) {
		sim(i, i) = 1;
	}

	for(int i = 0; i < m-1; i ++) {
		for(int j = i + 1; j < m; j ++) {
			sim(i, j) = m_intersect(i, j)/(v_occur[i] + v_occur[j] - m_intersect(i, j));
			sim(j, i) = sim(i, j);
		}
	}


	return sim;

}


