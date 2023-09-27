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
	
	IntegerVector all_ancestors = cpp_ancestors_of_a_group(dag, nodes, 1, true);
	LogicalVector l_offspring(n);
	LogicalVector l_all_ancestors = integer_to_logical_vector(all_ancestors - 1, n);

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
					ic[id1] += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[i]-1, l_all_ancestors);
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
							sim(id1, id2) += _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[i]-1, l_all_ancestors) +
							                 _calc_wang_s(lt_children, lt_children_relations, contribution, all_ancestors[k]-1, offspring[j]-1, l_all_ancestors);
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

const double PI = 3.1415926;


// [[Rcpp::export]]
NumericMatrix cpp_sim_shen(S4 dag, IntegerVector nodes, NumericVector ic) {
	int m = nodes.size();
	NumericMatrix sim(m, m);

	for(int i = 0; i < m; i ++) {
		if(std::abs(ic[ nodes[i]-1 ]) < 1e-10) {
			sim(i, i) = 0;
		} else {
			sim(i, i) = 1 - atan(1/ic[ nodes[i]-1 ])/PI*2;
		}
	}

	if(m <= 1) {
		return sim;
	}

	NumericVector weight;

	IntegerMatrix mica_nodes = cpp_max_ancestor_id(dag, nodes, ic);

	int i_pair = 0;
	int n_pairs = m*(m-1)/2;
	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {

			i_pair ++;
			if(i_pair % 1000 == 0) {
				message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
				message("going through " + std::to_string(i_pair) + " / " + std::to_string(n_pairs) + " pairs ...", false);
			}

			if(std::abs(ic[ mica_nodes(i, j)-1 ]) < 1e-10) {
				sim(i, j) = 0;
				sim(j, i) = 0;
			} else {
				IntegerVector path1 = cpp_tpl_shortest_path(dag, mica_nodes(i, j), nodes[i]);
				IntegerVector path2 = cpp_tpl_shortest_path(dag, mica_nodes(i, j), nodes[j]);

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
	}

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(n_pairs) + " / " + std::to_string(n_pairs) + " pairs ... Done.", true);

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

	int i_pair = 0;
	int n_pairs = m*(m-1)/2;
	for(int i = 0; i < m - 1; i ++) {
		for(int j = i+1; j < m; j ++) {

			i_pair ++;
			if(i_pair % 1000 == 0) {
				message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
				message("going through " + std::to_string(i_pair) + " / " + std::to_string(n_pairs) + " pairs ...", false);
			}

			IntegerVector path1 = cpp_tpl_shortest_path(dag, lca_nodes(i, j), nodes[i]);	
			IntegerVector path2 = cpp_tpl_shortest_path(dag, lca_nodes(i, j), nodes[j]);

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

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(n_pairs) + " / " + std::to_string(n_pairs) + " pairs ... Done.", true);


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
		_find_ancestors(lt_parents, nodes[i], la, true);
		m_ancestors(i, _) = la;
		reset_logical_vector_to_false(la);
	}

	LogicalVector ancestors1(n);
	LogicalVector ancestors2(n);

	for(int i = 0 ; i < m - 1; i ++) {
		for(int j = i + 1; j < m; j ++) {

			int n_ca = 0;
			double ic_sum = 0;

			// here looks for EICA (exclusively inherited and all common ancestors)
			for(int k = 0; k < n; k ++) {
				if(m_ancestors(i, k) && m_ancestors(j, k)) {

					IntegerVector children = lt_children[k];
					for(int p = 0; p < children.size(); p ++) {
						if( (m_ancestors(i, children[p] - 1 ) && !m_ancestors(j, children[p] - 1 )) || 
						    (!m_ancestors(i, children[p] - 1 ) && m_ancestors(j, children[p] - 1 )) ) {
							n_ca ++;
							ic_sum += ic[k];
						}
					}
				}
			}

			mean_ic(i, j) = ic_sum/n_ca;
			mean_ic(j, i) = mean_ic(i, j);
		}
	}

	return mean_ic;
}


bool _is_disjunctive_ancestor_pair(List lt_parents, int i_ancestor1, int i_ancestor2, int i_node) {
	int n = lt_parents.size();
	LogicalVector l_ancestors(n);
	LogicalVector l_background(n);

	// remove ancestor 1
	reset_logical_vector_to_true(l_background);
	l_background[i_ancestor1] = false;

	_find_ancestors_with_background(lt_parents, i_node, l_ancestors, l_background);
	if(!l_ancestors[i_ancestor2]) {
		return true;
	} else {
		reset_logical_vector_to_false(l_ancestors);
		reset_logical_vector_to_true(l_background);
		l_background[i_ancestor2] = false;

		_find_ancestors_with_background(lt_parents, i_node, l_ancestors, l_background);
		if(!l_ancestors[i_ancestor1]) {
			return true;
		}
	}

	return false;
}

// given a node, returns all disjunctive ancestor pairs
LogicalVector _disjunctive_common_ancestors_single(List lt_parents, int i_node1, int i_node2, NumericVector ic) {
	int n = lt_parents.size();

	LogicalVector l_ancestors1(n);
	LogicalVector l_ancestors2(n);
	LogicalVector l_background(n);

	_find_ancestors(lt_parents, i_node1, l_ancestors1, true);
	_find_ancestors(lt_parents, i_node2, l_ancestors2, true);

	LogicalVector l_common_ancestors = l_ancestors1 & l_ancestors2;

	LogicalVector l_DCA(n);

	for(int i = 0; i < n; i ++) { 
		if(l_common_ancestors[i]) { // for every a in common ancestor
			for(int j = 0; j < n; j ++) {
				if(l_common_ancestors[j]) { // go every c in common ancestor
					if(ic[i] < ic[j]) {
						if(!_is_disjunctive_ancestor_pair(lt_parents, i, j, i_node1)) {
							break;
						}
						if(!_is_disjunctive_ancestor_pair(lt_parents, i, j, i_node2)) {
							break;
						}
						l_DCA[i] = true;
					}
				}
			}
		}
	}

	return l_DCA;
}

// [[Rcpp::export]]
NumericMatrix cpp_common_ancestor_mean_IC_GraSM(S4 dag, IntegerVector nodes, NumericVector ic) {
	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();
	int m = nodes.size();

	NumericMatrix mean_ic(m, m);
	IntegerVector l_DCA(n);

	int i_pair = 0;
	int n_pairs = m*(m-1)/2 + m-1;
	for(int i = 0; i < m; i ++) {
		for(int j = i; j < m; j ++) {

			i_pair ++;
			if(i_pair % 1000 == 0) {
				message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
				message("going through " + std::to_string(i_pair) + " / " + std::to_string(n_pairs) + " pairs ...", false);
			}

			l_DCA = _disjunctive_common_ancestors_single(lt_parents, nodes[i] - 1, nodes[j] - 1, ic);
			int n_DCA = sum(l_DCA);
			if(n_DCA) {
				double ss = 0;
				for(int k = 0; k < l_DCA.size(); k ++) {
					if(l_DCA[k]) {
						ss += ic[k];
					}
				}

				mean_ic(i, j) = ss/n_DCA;
				mean_ic(j, i) = mean_ic(i, j);
			} else {
				mean_ic(i, j) = 0;
				mean_ic(j, i) = 0;
			}
		}
	}

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(n_pairs) + " / " + std::to_string(n_pairs) + " pairs ... Done.", true);

	return mean_ic;
}
