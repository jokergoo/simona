
#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"

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
