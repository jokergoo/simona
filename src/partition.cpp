#include <Rcpp.h>
using namespace Rcpp;

#include "traverse.h"

// [[Rcpp::export]]
IntegerVector cpp_partition_by_size(S4 tree, int size) {
	List lt_children = tree.slot("lt_children");
	int n = tree.slot("n_terms");
	int root = tree.slot("root");

	IntegerVector n_offspring = cpp_n_offspring_on_tree(tree, true);
	IntegerVector pa(n, -1);

	// breadth-first search
	IntegerVector current_nodes = {root};
	while(current_nodes.size()) {
		IntegerVector current_nodes2;
		for(int i = 0; i < current_nodes.size(); i ++) {
			int i_node = current_nodes[i] - 1;
			IntegerVector children = lt_children[i_node];
				
			if(children.size() == 0) { // leaf
				pa[i_node] = current_nodes[i];
			} else {

				if(n_offspring[i_node] <= size) {
					Rcout << "getting offsprings for " << current_nodes[i] << "\n";
					IntegerVector offspring = cpp_offspring(tree, current_nodes[i], true);
					pa[offspring-1] = current_nodes[i];
					continue;
				}

				// check i_node's children
				bool all_small_children = true;
				for(int j = 0; j < children.size(); j ++) {
					if(n_offspring[children[j]-1] >= size) {
						all_small_children = false;
					}
				}
				if(all_small_children) {
					IntegerVector offspring = cpp_offspring(tree, current_nodes[i], true);
					pa[offspring-1] = current_nodes[i];
				} else {
					for(int j = 0; j < children.size(); j ++) {
						if(n_offspring[children[j]-1] <= size) {
							IntegerVector offspring = cpp_offspring(tree, children[j], true);
							pa[offspring-1] = children[j];
						} else {
							current_nodes2.push_back(children[j]);
						}
					}
				}
			}
		}
		current_nodes = current_nodes2;
	}

	return pa;
}
