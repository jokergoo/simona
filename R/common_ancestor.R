


max_ancestor_v = function(dag, terms, value) {

	if(length(value) != dag@n_terms) {
		stop("Length of `value` should be the same as number of total terms in the DAG.")
	}

	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}

	an = cpp_max_ancestor_v(dag, id, value)
	dimnames(an) = list(dag@terms[id], dag@terms[id])
	an
}


max_ancestor_id = function(dag, terms, value, in_labels = FALSE) {

	if(length(value) != dag@n_terms) {
		stop("Length of `value` should be the same as number of total terms in the DAG.")
	}

	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}

	an = cpp_max_ancestor_id(dag, id, value)
	dimnames(an) = list(dag@terms[id], dag@terms[id])

	if(in_labels) {
		structure(dag@terms[an], dim = dim(an), dimnames = dimnames(an))
	} else {
		an
	}
}


#' Various types of common ancestors
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param IC_method An IC method. Valid values are in [ALL_IC_METHODS].
#' @param in_labels Whether the terms are in their names or as the integer indices?
#' 
#' @details
#' There are the following three types of common ancestors:
#' 
#' - MICA: most informative common ancestor: The common ancestor with the highest IC value.
#' - LCA: lowest common ancestor: The common ancestor with the largest height. If there are multiple ancestors having
#'        the max height, the ancestor with the smallest distance to the two terms is taken.
#' - NCA: nearest common ancestor: The common ancestor with the smallest distance to the two terms. If there are multiple
#'        ancestors with the same smallest distance, the ancestor with the largest height is taken.
#' 
#' @return 
#' - `MICA_term()` returns a character matrix of the MICA terms. 
#' - `MICA_IC()` returns a numeric matrix of the IC of the MICA terms.
#' - `LCA_term()` returns a character matrix of the LCA term.
#' - `LCA_depth()` reutrns an integer matrix of the depth of the LCA terms.
#' - `NCA_term()` returns a character matrix of the nearest common ancestors. The shortest distance can be calculated by `shortest_distances_via_CA()`.
#' @rdname common_ancestor
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' MICA_term(dag, letters[1:6], "IC_universal")
#' MICA_IC(dag, letters[1:6], "IC_universal")
#' LCA_term(dag, letters[1:6])
#' LCA_depth(dag, letters[1:6])
#' NCA_term(dag, letters[1:6])
MICA_term = function(dag, terms, IC_method, in_labels = TRUE) {
	ic = term_IC(dag, IC_method)
	max_ancestor_id(dag, terms, ic, in_labels = in_labels)
}

#' @rdname common_ancestor
#' @export
MICA_IC = function(dag, terms, IC_method) {
	ic = term_IC(dag, IC_method)
	max_ancestor_v(dag, terms, ic)
}


#' @rdname common_ancestor
#' @export
LCA_term = function(dag, terms, in_labels = TRUE) {
	depth = dag_depth(dag)
	
	max_ancestor_id(dag, terms, depth, in_labels = in_labels)
}

#' @rdname common_ancestor
#' @export
LCA_depth = function(dag, terms) {
	depth = dag_depth(dag)
	
	max_ancestor_v(dag, terms, depth)
}

#' @rdname common_ancestor
#' @export
NCA_term = function(dag, terms, in_labels = TRUE) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}

	an = cpp_nearest_common_ancestor(dag, id)
	dimnames(an) = list(dag@terms[id], dag@terms[id])

	if(in_labels) {
		structure(dag@terms[an], dim = dim(an), dimnames = dimnames(an))
	} else {
		an
	}

}



