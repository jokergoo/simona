


#' Various types of common ancestors
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param IC_method An IC method. Valid values are in [`all_term_IC_methods()`].
#' @param in_labels Whether the terms are represented in their names or as integer indices?
#' @param distance If there are multiple LCA or MICA of two terms, whether to take the one with
#'   the longest distance of shortest distance to the two terms. Possible values are "longest" and "shortest".
#' @param verbose Whether to print messages.
#' 
#' @details
#' There are the following three types of common ancestors:
#' 
#' - MICA (most informative common ancestor): The common ancestor with the highest IC value.
#' - LCA (lowest common ancestor): The common ancestor with the largest depth (The depth of a term is the maximal distance from the root term). If there are multiple ancestors having
#'        the same max depth, the ancestor with the smallest distance to the two terms is used.
#' - NCA (nearest common ancestor): The common ancestor with the smallest distance to the two terms. If there are multiple
#'        ancestors with the same smallest distance, the ancestor with the largest depth is used.
#' 
#' `max_ancestor_v()` and `max_ancestor_id()` are more general functions which return common ancestors with
#' the highest value in `value`.
#' 
#' Given a path connecting two terms and their MICA/LCA, `max_ancestor_path_sum()` calculates the sum of terms along the path. The values
#' to be added in specified in `add_v` argument.
#' 
#' @return 
#' - `MICA_term()` returns an integer or a character matrix of the MICA terms depending on the value of `in_labels`. 
#' - `MICA_IC()` returns a numeric matrix of the IC of the MICA terms.
#' - `LCA_term()` returns an integer or a character matrix of the LCA term depending on the value of `in_labels`.
#' - `LCA_depth()` returns an integer matrix of the depth of the LCA terms.
#' - `NCA_term()` returns an integer or a character matrix of the NCA term depending on the value of `in_labels`. The shortest distance from NCA terms can be calculated by [`shortest_distances_via_NCA()`].
#' - `max_ancestor_v()` returns a numeric matrix.
#' - `max_ancestor_id()` returns an integer or a character matrix.
#' 
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
MICA_term = function(dag, terms, IC_method, in_labels = TRUE, distance = "longest", verbose = simona_opt$verbose) {
	ic = term_IC(dag, IC_method, verbose = verbose)
	max_ancestor_id(dag, terms, ic, in_labels = in_labels, distance = distance, verbose = verbose)
}

#' @rdname common_ancestor
#' @export
MICA_IC = function(dag, terms, IC_method, verbose = simona_opt$verbose) {
	ic = term_IC(dag, IC_method, verbose = verbose)
	max_ancestor_v(dag, terms, ic, verbose = verbose)
}


#' @rdname common_ancestor
#' @export
LCA_term = function(dag, terms, in_labels = TRUE, distance = "longest", verbose = simona_opt$verbose) {
	depth = dag_depth(dag)
	
	max_ancestor_id(dag, terms, depth, in_labels = in_labels, distance = distance, verbose = verbose)
}

#' @rdname common_ancestor
#' @export
LCA_depth = function(dag, terms, verbose = simona_opt$verbose) {
	depth = dag_depth(dag)
	
	max_ancestor_v(dag, terms, depth, verbose = verbose)
}

#' @rdname common_ancestor
#' @export
NCA_term = function(dag, terms, in_labels = TRUE, verbose = simona_opt$verbose) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}

	an = exec_under_message_condition({
		cpp_nearest_common_ancestor(dag, id)
	}, verbose = verbose)
	dimnames(an) = list(dag@terms[id], dag@terms[id])

	if(in_labels) {
		structure(dag@terms[an], dim = dim(an), dimnames = dimnames(an))
	} else {
		an
	}

}


#' @param value A numeric vector. The elements should corrrespond to terms in `dag_all_terms()` (should have the same length as the number of terms in the DAG).
#' 
#' @rdname common_ancestor
#' @export
max_ancestor_v = function(dag, terms, value, verbose = simona_opt$verbose) {

	if(length(value) != dag@n_terms) {
		stop("Length of `value` should be the same as number of total terms in the DAG.")
	}

	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	an = exec_under_message_condition({
		cpp_max_ancestor_v(dag, id, value)
	}, verbose = verbose)
	dimnames(an) = list(dag@terms[id], dag@terms[id])
	an
}


#' @rdname common_ancestor
#' @export
max_ancestor_id = function(dag, terms, value, in_labels = FALSE, distance = "longest", verbose = simona_opt$verbose) {

	if(length(value) != dag@n_terms) {
		stop("Length of `value` should be the same as number of total terms in the DAG.")
	}

	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	an = exec_under_message_condition({
		cpp_max_ancestor_id(dag, id, value, distance == "longest")
	}, verbose = verbose)
	dimnames(an) = list(dag@terms[id], dag@terms[id])

	if(in_labels) {
		structure(dag@terms[an], dim = dim(an), dimnames = dimnames(an))
	} else {
		an
	}
}


#' @param add_v Values to be added along the path to the MICA or LCA. The same format as `value`.
#' 
#' @export
#' @rdname common_ancestor
max_ancestor_path_sum = function(dag, terms, value, add_v, distance = "longest", verbose = simona_opt$verbose) {
	if(length(value) != dag@n_terms) {
		stop("Length of `value` should be the same as number of total terms in the DAG.")
	}
	if(length(add_v) != dag@n_terms) {
		stop("Length of `add_v` should be the same as number of total terms in the DAG.")
	}

	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}

	sv = exec_under_message_condition({
		cpp_max_ancestor_path_sum_value(dag, id, value, add_v, distance == "longest")
	}, verbose = verbose)

	dimnames(sv) = list(dag@terms[id], dag@terms[id])
	sv
}

