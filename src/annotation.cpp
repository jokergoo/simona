#include <Rcpp.h>
using namespace Rcpp;

#include "transverse.h"
#include "utils.h"

// [[Rcpp::export]]
IntegerVector cpp_n_annotations(S4 dag) {

	List lt_children = dag.slot("lt_children");
	List annotation = dag.slot("annotation");
	List lt_annotation = annotation["list"];
	CharacterVector anno_names = annotation["names"];
	int n_all_anno = anno_names.size();

	int n = lt_children.size();
	IntegerVector n_anno(n, 0);

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < n; i ++) {
		_find_offspring(lt_children, i, l_offspring, true);  //include self

		LogicalVector l_anno(n_all_anno, false);
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) {
				IntegerVector anno = lt_annotation[j];
				for(int k = 0; k < anno.size(); k ++) {
					l_anno[anno[k]-1] = true;
				}
			}
		}
		n_anno[i] = sum(l_anno);

		reset_logical_vector_to_false(l_offspring);
	}

	return n_anno;
}


// [[Rcpp::export]]
IntegerMatrix cpp_get_term_annotations(S4 dag, IntegerVector nodes) {
	List lt_children = dag.slot("lt_children");
	List annotation = dag.slot("annotation");
	List lt_annotation = annotation["list"];
	CharacterVector anno_names = annotation["names"];
	int n_all_anno = anno_names.size();
	int n = lt_children.size();
	int m = nodes.size();

	IntegerMatrix mat(m, n_all_anno);

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < m; i ++) {
		_find_offspring(lt_children, nodes[i]-1, l_offspring, true);  //include self

		LogicalVector l_anno(n_all_anno, false);
		for(int j = 0; j < n; j ++) {
			if(l_offspring[j]) {
				IntegerVector anno = lt_annotation[j];
				for(int k = 0; k < anno.size(); k ++) {
					l_anno[anno[k]-1] = true;
					mat(i, anno[k]-1) = 1;
				}
			}
		}

		reset_logical_vector_to_false(l_offspring);
	}

	return mat;
}


IntegerVector cpp_get_annotated_terms(S4 dag, int anno_id) {
	List lt_children = dag.slot("lt_children");
	List lt_parents = dag.slot("lt_parents");

	List annotation = dag.slot("annotation");
	List lt_annotation = annotation["list"];

	int n = lt_parents.size();

	IntegerVector x(n);
	LogicalVector l_ancestors(n, false);
	for(int i = 0; i < n; i ++) {
		IntegerVector anno = lt_annotation[i];
		if(anno.size() || x[i] == 0) {
			for(int ia = 0; ia < anno.size(); ia ++) {
				if(anno[ia] == anno_id) {
					_find_ancestors(lt_parents, i, l_ancestors, true);	
					for(int j = 0; j < n; j ++) {
						if(l_ancestors[j]) {
							x[j] = 1;
						}
					}
					break;
				}
			}
		}

		reset_logical_vector_to_false(l_ancestors);
	}

	return x;
}

// [[Rcpp::export]]
IntegerMatrix cpp_get_annotated_terms(S4 dag, IntegerVector anno_id) {
	int m = anno_id.size();
	int n = dag.slot("n_terms");
	IntegerMatrix x(m, n);

	for(int i = 0; i < m; i ++) {
		x(i, _) = cpp_get_annotated_terms(dag, anno_id[i]);
	}

	return x;

}
