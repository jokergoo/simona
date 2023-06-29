#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"
#include "transverse.h"


NumericVector _get_breaks(double left, double right, int n, NumericVector weight) {
	NumericVector breaks(n+1);
	double rg = right - left;
	breaks[0] = left;
	breaks[n] = right;

	if(n <= 1) {
		return breaks;
	}

	NumericVector bin(n);
	double weight_sum = sum(weight);

	for(int i = 0; i < n; i ++) {
		bin[i] = weight[i]/weight_sum*rg;
	}

	for(int i = 1; i <= n; i ++) {
		breaks[i] = breaks[i-1] + bin[i-1];
	}

	return breaks;
}

// [[Rcpp::export]]
DataFrame cpp_term_pos_on_circle(S4 dag, double start = 0, double end = 360) {

	int root = dag.slot("root");
	List lt_children = dag.slot("lt_children");
	IntegerVector depth = _dag_depth(dag);
	int max_depth = max(depth);
	IntegerVector n_offspring = cpp_n_offspring(dag);
	int n = lt_children.size();

	NumericVector theta(n);  // in degree
	NumericVector rho(n);    // radius
	NumericVector width(n);   // width of each term

	// root
	int i_root = root - 1;
	LogicalVector l_current_parent(n);
	LogicalVector l_current_parent2(n);

	NumericVector parent_range_left(n);
	NumericVector parent_range_right(n);

	l_current_parent[i_root] = true;
	parent_range_left[i_root] = 0;
	parent_range_right[i_root] = 360;

	int current_depth = 0;

	// polar coordinates of root
	theta[i_root] = (end + start)*0.5;
	rho[i_root] = 0;
	width[i_root] = end - start;

	if(max_depth == 1) {
		DataFrame df = DataFrame::create(Named("theta") = theta, Named("rho") = rho, Named("level1_group") = -1, Named("width") = width);

		return df;
	}

	// level 1
	IntegerVector level1 = lt_children[i_root];
	IntegerVector level1_group(n);
	for(int i = 0; i < level1.size(); i ++) {
		level1_group[ level1[i]-1 ] = level1[i];
	}

	while(1) {
		current_depth = current_depth + 1;
		reset_logical_vector_to_false(l_current_parent2);

		for(int i = 0; i < n; i ++) {
			if(l_current_parent[i]) { // check its children
				IntegerVector children = lt_children[i];
				IntegerVector children2;

				for(int j = 0; j < children.size(); j ++) {
					if(depth[ children[j]-1 ] == current_depth) {
						children2.push_back(children[j]);
					}
				}

				// calculate the circular coordinate of #children2 nodes
				int n_children2 = children2.size();
				NumericVector weight(n_children2);
				for(int j = 0; j < n_children2; j ++) {
					if(n_offspring[children2[j]-1] == 0) {
						weight[j] = 1;
					} else {
						weight[j] = n_offspring[children2[j]-1];
					}
				}
				NumericVector breaks = _get_breaks(parent_range_left[i], parent_range_right[i], n_children2, weight);
				for(int j = 0; j < n_children2; j ++) {
					theta[ children2[j]-1 ] = (breaks[j] + breaks[j+1])/2;
					rho[ children2[j]-1 ] = current_depth;
					width[ children2[j]-1 ] = breaks[j+1] - breaks[j];

					parent_range_left[ children2[j]-1 ] = breaks[j];
					parent_range_right[ children2[j]-1 ] = breaks[j+1];

					if(current_depth > 1) {
						level1_group[ children2[j]-1 ] = level1_group[i];
					}
				}

				for(int j = 0; j < n_children2; j ++) {
					l_current_parent2[ children2[j]-1 ] = true;
				}
			}
		}
		
		for(int i = 0; i < n; i ++) {
			l_current_parent[i] = l_current_parent2[i];
		}

		if(current_depth == max_depth) {
			break;
		}
	}

	DataFrame df = DataFrame::create(Named("theta") = theta, Named("rho") = rho, Named("level1_group") = level1_group, Named("width") = width);

	return df;
}



