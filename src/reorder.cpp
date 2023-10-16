
#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"
#include "traverse.h"


// we have a complete DAG and a reduced tree, the force comes from the additional links in DAG while not in tree,
// the force on a node is defined as the total force from additional links applied to the sub-tree rooted by the node.
//
// we consider node t's parent and its offspring's parents that contribute to the "additional links" to the sub-tree
// [[Rcpp::export]]
List cpp_get_force_counterpart(
	List lt_children_dag,
	List lt_parents_dag,
	List lt_children_tree,
	List lt_parents_tree,
	int root) {

	int n = lt_children_dag.size();

	LogicalVector l_offspring(n);
	LogicalVector l_counterpart(n);

	List lt_counterpart(n);

	int i_root = root - 1;

	for(int i = 0; i < n; i ++) { // for each node
		if(i == i_root) {
			lt_counterpart[i] = IntegerVector(0);
			continue;
		}

		// its parents, children, offspring, and offspring's parents
		_find_offspring(lt_children_tree, i, l_offspring, true);
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) { // for each offspring
				IntegerVector parents = lt_parents_dag[j];
				if(parents.size()) {
					l_counterpart[ parents - 1 ] = true;
				}
			}
		}

		IntegerVector children = lt_children_dag[i];
		if(children.size()) {
			l_counterpart[children - 1] = true;
		}

		// remove those overlap with the tree
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) {
				l_counterpart[j] = false;
			}
		}
		IntegerVector parents_tree;
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) { // for each offspring
				parents_tree = lt_parents_tree[j];
				for(int k = 0; k < parents_tree.size(); k ++) {
					l_counterpart[ parents_tree[k]-1 ] = false;
				}
			}
		}

		l_counterpart[i] = false;  // the node itself is removed
		l_counterpart[i_root] = false;  // root is removed

		IntegerVector counterpart = _which(l_counterpart) + 1; // be consistent to other lt_*, the index starts from 1
		lt_counterpart[i] = counterpart;

		reset_logical_vector_to_false(l_offspring);
		reset_logical_vector_to_false(l_counterpart);
	}

	return lt_counterpart;
}


// the force is defined as  sum(x_diff*depth), x_diff is directional, depth is on the tree
// on circle, x_diff have a different way to calculating, taking consider of the circle
// [[Rcpp::export]]
NumericVector cpp_get_force(List lt_counterpart, NumericVector x, IntegerVector depth, bool on_circle = false) {
	int n = lt_counterpart.size();
	NumericVector force(n);

	double f = 0;
	double diff;

	int j;

	for(int i = 0; i < n; i ++) {
		IntegerVector counterpart = lt_counterpart[i];

		if(counterpart.size()) {
			if(on_circle) {
				for(int k = 0; k < counterpart.size(); k ++) {
					j = counterpart[k] - 1;
					diff = x[j] - x[i];
					diff = diff - int(diff)/360*360;
					if(diff < 0) {
						diff = 360 + diff;
					}

					if(diff > 180) {
						diff = 180 - diff;
					}

					f += diff * depth[i];
					
				}
			} else {
				for(int k = 0; k < counterpart.size(); k ++) {
					j = counterpart[k] - 1;
					// if(depth[j] > depth[i]) { // deeper ones pull the upper ones
						f += (x[j] - x[i]) * depth[i];
					// }
				}
			}
		}

		force[i] = f;

		f = 0;
	}

	return force;
}



// x: value vector, only look at the sign, positve values moved to the left, negative values movec to the right
// sorted_od: the order of abs(x) in decreasing order, ie, order(-abs(x))
// k: how many top value to move
// [[Rcpp::export]]
IntegerVector move_index(NumericVector x, IntegerVector sorted_od, int k, bool decreasing = true) {
	int n = x.size();
	IntegerVector new_od(n, -1);

	LogicalVector l_unused(n, true);
	int i_left = 0;
	int i_right = n-1;
	for(int i = 0; i < k; i ++) {
		if(decreasing) {
			if(x[sorted_od[i]] > 0) {
				new_od[i_left] = sorted_od[i];
				i_left ++;
			} else if(x[sorted_od[i]] < 0) {
				new_od[i_right] = sorted_od[i];
				i_right --;
			}
		} else {
			if(x[sorted_od[i]] < 0) {
				new_od[i_left] = sorted_od[i];
				i_left ++;
			} else if(x[sorted_od[i]] > 0) {
				new_od[i_right] = sorted_od[i];
				i_right --;
			}
		}
		if(x[sorted_od[i]] != 0) {
			l_unused[sorted_od[i]] = false;
		}
	}

	IntegerVector unused = _which(l_unused);
	for(int i = 0; i < unused.size(); i ++) {
		new_od[i + i_left] = unused[i];
	}

	return new_od;
}


// children: a vector of child IDs
// prev_od, new_od: two orderings for `children`
// width: the global vector of width where the width for children can be obtained by `width[children-1]
// [[Rcpp::export]]
NumericVector calc_x_offset(IntegerVector children, IntegerVector prev_od, IntegerVector new_od, NumericVector width) {

	int n = children.size();
	
	IntegerVector children1 = children[prev_od];
	IntegerVector children2 = children[new_od];

	NumericVector w1 = width[children1 - 1];
	NumericVector w2 = width[children2 - 1];

	IntegerVector x1 = IntegerVector(n + 1);
	IntegerVector x2 = IntegerVector(n + 1);

	for(int i = 0; i < n; i ++) {
		x1[i+1] += w1[i] + x1[i];
		x2[i+1] += w2[i] + x2[i];
	}

	IntegerVector match_ind1 = _order(prev_od);
	IntegerVector match_ind2 = _order(new_od);
	
	NumericVector mid1(n);
	NumericVector mid2(n);
	for(int i = 0; i < n; i ++) {
		mid1[i] = (x1[match_ind1[i]+1] + x1[match_ind1[i]])/2;
		mid2[i] = (x2[match_ind2[i]+1] + x2[match_ind2[i]])/2;
	}

	NumericVector offset(n);
	for(int i = 0; i < n; i ++) {
		offset[i] = mid2[i] - mid1[i];
	}

	return offset; // the same order as `children`
}


// two objects to modify
// - children a reordered version of children
// - new_x an updated version of `x`.
// [[Rcpp::export]]
IntegerVector reorder_children(IntegerVector children, IntegerVector n_cp, NumericVector force, 
	NumericVector width, IntegerVector depth, NumericVector new_x, List lt_children) {
	// reorder children and calculate sum of forces

	int n = lt_children.size();
	int nc = children.size();
	NumericVector f = force[children - 1];

	bool f_all_zero = true;
	for(int i = 0; i < nc; i ++) {
		if(std::abs(f[i]) > 100) {
			f_all_zero = false;
			break;
		}
	}
	if(f_all_zero) {
		return children;
	}

	IntegerVector sorted_od = _order(-abs(f));

	IntegerVector prev_od = seq(0, nc - 1);
	IntegerVector new_od(nc);

	NumericVector x_diff(nc);
	NumericVector force_diff(nc);

	double total_force_diff = 0;

	bool flag = false;
	for(int i = 0; i < nc; i ++) {

		if(n_cp[ children[sorted_od[i]] - 1 ] == 0) {
			continue;
		}

		flag = true;

		new_od = move_index(f, sorted_od, i+1, false); // move top i+1 items to the two sides
		x_diff = calc_x_offset(children, prev_od, new_od, width);

		for(int j = 0; j < nc; j ++) {
			force_diff[j] = x_diff[j]*n_cp[children[j]-1]*depth[children[j]-1];
		}

		total_force_diff = 0;
		for(int j = 0; j < nc; j ++) {
			total_force_diff += std::abs(force[j] - force_diff[j]) - std::abs(force[j]);
		}

		prev_od = new_od;
	}

	if(!flag) {
		return children;
	}

	prev_od = seq(0, nc - 1);
	new_od = move_index(f, sorted_od, nc, false); // move top i+1 items to the two sides

	x_diff = calc_x_offset(children, prev_od, new_od, width);

	for(int j = 0; j < nc; j ++) {
		LogicalVector l_offspring(n);
		_find_offspring(lt_children, children[j]-1, l_offspring, true);

		for(int k = 0; k < n; k ++) {
			if(l_offspring[k]) {
				new_x[k] += x_diff[j];
			}
		}
	}

	IntegerVector new_children(nc);
	for(int j = 0; j < nc; j ++) {
		new_children[j] = children[ new_od[j] ];
	}
	return new_children;

}


// the aim of this function is to reorder children of a node, or get a new set of x-positions
// that correspond to the reordered children (note the order of `x` is not changed, always from the first node to the last node).
// `lt_children` can be reordered or not. If not, it can be reordered later by new values in `x`.
// `width` is the width of each node, so it is used in the reordering procedure
// [[Rcpp::export]]
NumericVector cpp_reorder_tree_x(S4 tree, List lt_counterpart, NumericVector x, NumericVector width, int times = 1) {
	int n = tree.slot("n_terms");
	List lt_children = tree.slot("lt_children");
	IntegerVector depth = _dag_depth(tree);

	NumericVector force = cpp_get_force(lt_counterpart, x, depth);

	IntegerVector n_cp(n);
	for(int i = 0; i < n; i ++) {
		IntegerVector cp = lt_counterpart[i];
		n_cp[i] = cp.size();
	}

	double total_force = sum(abs(force));
	double total_force2 = 0;

	NumericVector new_x = clone(x);
	NumericVector new_x2;

	for(int k = 0; k < times; k ++) {

		new_x2 = clone(new_x);
		for(int i = 0; i < n; i ++) {
			IntegerVector children = lt_children[i];
			IntegerVector children2(children.size());

			if(children.size() > 1) {
				// adjusted x-positions are updated in `new_x`.
				children2 = reorder_children(children, n_cp, force, width, depth, new_x, lt_children);
				lt_children[i] = children2;

			} 
		}

		force = cpp_get_force(lt_counterpart, new_x, depth);
		total_force2 = sum(abs(force));
		Rcout << k << ": " << total_force << "/" << total_force2 << "\n";
		if(total_force2 >= total_force) {
			return new_x2;
		}
		total_force = total_force2;
	}

	return new_x;
}


// [[Rcpp::export]]
int n_links_from_two_groups_of_nodes(S4 dag, IntegerVector nodes1, IntegerVector nodes2) {

	// number of cross-cluster links
	List lt_children = dag.slot("lt_children");
	List lt_parents = dag.slot("lt_parents");
	int n = lt_children.size();

	IntegerVector i_nodes1 = nodes1 - 1;
	IntegerVector i_nodes2 = nodes2 - 1;

	int k = 0;

	LogicalVector l2 = integer_to_logical_vector(i_nodes2, n);

	for(int i = 0; i < i_nodes1.size(); i ++) {
		IntegerVector children = lt_children[ i_nodes1[i] ];
		if(children.size()) {
			for(int j = 0; j < children.size(); j ++) {
				if(l2[ children[j]-1 ]) {
					k ++;
				}
			}
		}

		IntegerVector parents = lt_parents[ i_nodes1[i] ];
		if(parents.size()) {
			for(int j = 0; j < parents.size(); j ++) {
				if(l2[ parents[j]-1 ]) {
					k ++;
				}
			}
		}
	}

	return k;
}
