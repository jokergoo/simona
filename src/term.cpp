#include <Rcpp.h>
using namespace Rcpp;

int cpp_dag_depth_calc_d(int node, IntegerVector d, List lt_parents) {

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
		v[i] = cpp_dag_depth_calc_d(parents[i], d, lt_parents);
	}
	d[i_node] = max(v) + 1;

	return(d[i_node]);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_depth(List lt_parents, int nv) {
	IntegerVector d(nv, -1);

	for(int node = 1; node <= nv; node ++) {
		cpp_dag_depth_calc_d(node, d, lt_parents);
	}

	return(d);
}

int cpp_dag_height_calc_d(int node, IntegerVector d, List lt_children) {
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
		v[i] = cpp_dag_height_calc_d(children[i], d, lt_children);
	}
	d[i_node] = max(v) + 1;

	return(d[i_node]);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_height(List lt_children, int nv) {
	IntegerVector d(nv, -1);

	for(int node = 1; node <= nv; node ++) {
		cpp_dag_height_calc_d(node, d, lt_children);
	}

	return(d);
}
