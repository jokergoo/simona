#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix cpp_common_ancestor_int(List lt_offspring, IntegerVector v, List lt_parents) {
	int n = lt_offspring.size();
	IntegerMatrix m(n, n);

	for(int k = 0; k < n; k ++) {  
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}

		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) { // the first element is the node itself, so ignore it
				for(int j = i+1; j < noff; j ++) {
					int id1 = offspring[i]-1;
					int id2 = offspring[j]-1;
					
					m(id1, id2) = v[k] > m(id1, id2)? v[k] : m(id1, id2);
					m(id2, id1) = m(id1, id2);
					
				}
			}
		}

		m(k, k) = v[k];
	}

	return m;
}

// [[Rcpp::export]]
NumericMatrix cpp_common_ancestor_double(List lt_offspring, NumericVector v, List lt_parents) {
	int n = lt_offspring.size();
	NumericMatrix m(n, n);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					int id1 = offspring[i]-1;
					int id2 = offspring[j]-1;
					
					m(id1, id2) = v[k] > m(id1, id2)? v[k] : m(id1, id2);
					m(id2, id1) = m(id1, id2);
				}
			}
		}

		m(k, k) = v[k];
	}

	return m;
}

// [[Rcpp::export]]
IntegerMatrix cpp_common_ancestor_ID(List lt_offspring, NumericVector v, List lt_parents, List dist_offspring) {
	int n = lt_offspring.size();
	IntegerMatrix m(n, n);
	IntegerMatrix dsumsq(n, n);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		IntegerVector dist = dist_offspring[k];

		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					int id1 = offspring[i]-1;
					int id2 = offspring[j]-1;
					if(m[id1, id2] == 0) {  // first time to go over `offspring`
						m(id1, id2) = k + 1;
						m(id2, id1) = m(id1, id2);

						dsumsq(id1, id2) = pow(dist[id1], 2) + pow(dist[id2], 2);
						dsumsq(id2, id1) = dsumsq(id1, id2);
					} else if(v[k] > v[ m[id1, id2]-1 ]) {
						m(id1, id2) = k + 1;
						m(id2, id1) = m(id1, id2);

						dsumsq(id1, id2) = pow(dist[id1], 2) + pow(dist[id2], 2);
						dsumsq(id2, id1) = dsumsq(id1, id2);
					
					} else if(v[k] == v[ m[id1, id2]-1 ]) {
						float current_dsumsq = pow(dist[id1], 2) + pow(dist[id2], 2);
						if(current_dsumsq < dsumsq(id1, id2) && dsumsq(id1, id2) > 0) {
							m(id1, id2) = k + 1;
							m(id2, id1) = m(id1, id2);

							dsumsq(id1, id2) = current_dsumsq;
							dsumsq(id2, id1) = dsumsq(id1, id2);
						}
					}
				}
			}
		}

		m(k, k) = k + 1;

	}

	return m;
}



NumericMatrix cpp_eps_EICA(List lt_offspring, List lt_ancestor, List, lt_children, NumericVector ic) {
	int n = lt_offspring.size();
	NumericMatrix s(n, n);
	IntegerMatrix n_ca(n, n);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					int id1 = offspring[i]-1;
					int id2 = offspring[j]-1;

					it = intersect(lt_children[k], setdiff(union_(lt_ancestor[id1], lt_ancestor[id2]), intersect(lt_ancestor[id1], lt_ancestor[id2])));
					
					if(it.size() > 0) {
						s(id1, id2) = s(id1, id2) + ic[k];
						s(id2, id1) = s(id1, id2);

						if(ic[k] > 0) {
							n_ca(id1, id2) += 1;
							n_ca(id2, id1) = n_ca(id1, id2);
						}
					}
				}
			}
		}
	}
	return s/n_ca;
}


// [[Rcpp::export]]
NumericMatrix cpp_eps_XGraSM(List lt_offspring, NumericVector ic) {
	int n = lt_offspring.size();
	NumericMatrix s(n, n);
	IntegerMatrix n_ca(n, n);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();

		if(noff == 0) {
			continue;
		}
		
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					int id1 = offspring[i]-1;
					int id2 = offspring[j]-1;
					
					s(id1, id2) = s(id1, id2) + ic[k];
					s(id2, id1) = s(id1, id2);

					if(ic[k] > 0) {
						n_ca(id1, id2) += 1;
						n_ca(id2, id1) = n_ca(id1, id2);
					}
				}
			}
		}
	}
	return s/n_ca;
}



// [[Rcpp::export]]
NumericMatrix cpp_ancestor_aggregate_wang(NumericMatrix sv, List lt_offspring) {
	int n = lt_offspring.size();
	NumericMatrix m(n, n);
	m.fill_diag(1);

	NumericVector ic = colSums(sv);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					m(i, j) += (sv[k, i] + sv[k, j])/(ic[i] + ic[j]);
					m(j, i) = m(i, j);
				}
			}
		}
	}

	return m;
}

// [[Rcpp::export]]
NumericMatrix cpp_ancestor_aggregate_aic(List lt_offspring, NumericVector sv, NumericVector sw) {
	int n = lt_offspring.size();
	NumericMatrix m(n, n);
	m.fill_diag(1);

	for(int k = 0; k < n; k ++) {
		IntegerVector offspring = lt_offspring[k];
		int noff = offspring.size();
		if(noff > 1) {
			for(int i = 0; i < noff - 1; i ++) {
				for(int j = i+1; j < noff; j ++) {
					m(i, j) += 2*sw[k]/(sv[i] + sv[j]);
					m(j, i) = m(i, j);
				}
			}
		}
	}

	return m;
}



// [[Rcpp::export]]
NumericMatrix cpp_sim_shen(IntegerMatrix mica, IntegerVector d) {
	int n = mica.size();
	NumericMatrix sim(n, n);

	for(int i = 0; i < n-1; i ++) {
		for(int j = i+1; j < n; j ++) {
			int mica_term = mica[i, j];
			sim(i, j) = 1 - atan(d[mica_term, i] + d[mica_term, j])/3.1415926*2;
		}
	}
	return sim;
}

// [[Rcpp::export]]
NumericMatrix cpp_sim_pekar(IntegerMatrix lca_term, IntegerMatrix d, IntegerVector depth, int factor) {
	int n = lca_term.size();
	NumericMatrix sim(n, n);
	sum.fill_diag(1);

	if(n == 1) {
		return m;
	}
	
	for(int i = 0; i < n-1; i ++) {
		for(int j = i+1; j < n; j ++) {
			int c = lca_term[i, j];
			sim(i, j) = (factor*depth[c] + 0.0) / (factor*depth[c] + d[c, i] + d[c, j] + 0.0);
			sim(j, i) = sim(i, j);
		}
	}
	return sim;
}

// [[Rcpp::export]]
NumericMatrix cpp_sim_wang_edge(IntegerMatrix lca_term, IntegerMatrix d, IntegerVector depth) {
	int n = lca_term.size();
	NumericMatrix sim(n, n);
	sim.fill_diag(1);

	if(n == 1) {
		return m;
	}

	for(int i = 0; i < n-1; i ++) {
		for(int j = i+1; j < n; j ++) {
			int c = lca_term[i, j];
			sim(i, j) = (depth[c] + 0.0)*(depth[c] + 0.0) / (depth[c] + d[c, i] + 0.0)/ (depth[c] + d[c, j] + 0.0);
			sim(j, i) = sim(i, j);
		}
	}
	return sim;
}


