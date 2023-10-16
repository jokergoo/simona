#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"
#include "traverse.h"


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
DataFrame cpp_node_pos_in_tree(S4 tree, IntegerVector bin_size, double start = 1, double end = 360) {

	int root = tree.slot("root");
	List lt_children = tree.slot("lt_children");
	IntegerVector depth = _dag_depth(tree);
	int max_depth = max(depth);
	int n = lt_children.size();

	NumericVector x(n);  // horizontal direction if the tree is top-down
	NumericVector h(n);    // vertical direction
	NumericVector width(n);   // width of each term

	// root
	int i_root = root - 1;
	LogicalVector l_current_parent(n);
	LogicalVector l_current_parent2(n);

	NumericVector parent_range_left(n);
	NumericVector parent_range_right(n);

	l_current_parent[i_root] = true;
	parent_range_left[i_root] = start;
	parent_range_right[i_root] = end;

	int current_depth = 0;

	x[i_root] = (end + start)*0.5;
	h[i_root] = 0;
	width[i_root] = end - start;

	if(max_depth == 0) {
		DataFrame df = DataFrame::create(Named("x") = x, Named("h") = h, Named("level1_group") = -1, Named("width") = width);

		return df;
	}

	String msg;
	int i_visited = 1;
	while(1) {
		current_depth = current_depth + 1;
		reset_logical_vector_to_false(l_current_parent2);

		for(int i = 0; i < n; i ++) {
			if(l_current_parent[i]) { // check its children
				IntegerVector children = lt_children[i];
				LogicalVector l_children(children.size());

				for(int j = 0; j < children.size(); j ++) {
					if(depth[ children[j]-1 ] == current_depth) {
						l_children[j] = true;
					}
				}

				IntegerVector children2 = children[l_children];

				// calculate the circular coordinate of #children2 nodes
				int n_children2 = children2.size();
				NumericVector weight(n_children2);
				for(int j = 0; j < n_children2; j ++) {
					if(bin_size[children2[j]-1] == 0) {
						weight[j] = 1;
					} else {
						weight[j] = bin_size[children2[j]-1];
					}
				}
				NumericVector breaks = _get_breaks(parent_range_left[i], parent_range_right[i], n_children2, weight);
				for(int j = 0; j < n_children2; j ++) {
					x[ children2[j]-1 ] = (breaks[j] + breaks[j+1])/2;
					h[ children2[j]-1 ] = current_depth;
					width[ children2[j]-1 ] = breaks[j+1] - breaks[j];

					parent_range_left[ children2[j]-1 ] = breaks[j];
					parent_range_right[ children2[j]-1 ] = breaks[j+1];

					i_visited ++;

					if(i_visited % 1000 == 0) {
						message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
						message("going through " + std::to_string(i_visited) + " / " + std::to_string(n) + " nodes ...", false);
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

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(n) + " / " + std::to_string(n) + " nodes ... Done.", true);

	DataFrame df = DataFrame::create(Named("x") = x, Named("h") = h, Named("width") = width);

	return df;
}

// [[Rcpp::export]]
IntegerVector cpp_calc_n_neighbours(NumericVector x, double range) {

	 // x is sorted and in [0, 360], range is small
	int n = x.size();

	IntegerVector k(n);

	int prev_i;
	int next_i;
	double diff;
	bool flag = false;
	for(int i = 0; i < n; i ++) {
		k[i] = 1;

		// forward
		flag = false;
		prev_i = i - 1;
		while(true) {
			if(prev_i < 0) {
				prev_i = n - 1;
				flag = true;
			} else {
				prev_i --;
			}

			if(flag) {
				diff = x[i] - x[prev_i] + 360;
			} else {
				diff = x[i] - x[prev_i];
			}
			if(diff < range) {
				k[i] ++;
			} else {
				break;
			}
		}

		// afterword
		flag = false;
		next_i = i + 1;
		while(true) {
			if(next_i > n - 1) {
				next_i = 0;
				flag = true;
			} else {
				next_i ++;
			}

			if(flag) {
				diff = x[next_i] - x[i] + 360;
			} else {
				diff = x[next_i] - x[i];
			}
			if(diff < range) {
				k[i] ++;
			} else {
				break;
			}
		}
	}

	return k;
}

