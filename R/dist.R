
#' Distance on the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param verbose Whether to print messages.
#' 
#' @details
#' Denote two terms as `a` and `b`, a common ancestor as `c`, and the distance function `d()` calculates the longest
#' distance or the shortest distance depending on the function.
#' 
#' - `shortest_distances_via_NCA()`: It calculates the smallest `d(c, a) + d(c, b)` where `d()` calculates the shortest distance between two terms. In this case,
#'     `c` is the NCA (nearest common ancestor) of `a` and `b`.
#' - `longest_distances_via_LCA()`: It calculates the largest `d(c, a) + d(c, b)` where `d()` calculates the longest distance between two terms *via the LCA (lowest common ancestor) term*. In this case,
#'     `c` is the LCA of `a` and `b`.
#' - `shortest_distances_directed()`: It calculates `d(a, b)` where `d()` calculates the shortest distance between two terms. The distance is only calculated when `a` is an ancestor of `b`, otherwise the distance value is -1.
#' - `longest_distances_directed()`: It calculates `d(a, b)` where `d()` calculates the longest distance between two terms. The distance is only calculated when `a` is an ancestor of `b`, otherwise the distance value is -1.
#' @rdname distance
#' @export
#' @returns A numeric distance matrix.
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' shortest_distances_via_NCA(dag, letters[1:6])
#' longest_distances_via_LCA(dag, letters[1:6])
#' shortest_distances_directed(dag, letters[1:6])
#' longest_distances_directed(dag, letters[1:6])
shortest_distances_via_NCA = function(dag, terms, verbose = simona_opt$verbose) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	d = exec_under_message_condition({
		cpp_shortest_distances_via_NCA(dag, id)
	}, verbose = verbose)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
longest_distances_via_LCA = function(dag, terms, verbose = simona_opt$verbose) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	d = exec_under_message_condition({
		cpp_max_ancestor_path_sum_value(dag, id, dag_depth(dag), rep(1, dag@n_terms)) - 1
	}, verbose = verbose)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
shortest_distances_directed = function(dag, terms, verbose = simona_opt$verbose) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	d = exec_under_message_condition({
		cpp_shortest_distances_directed(dag, id)
	}, verbose = verbose)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
longest_distances_directed = function(dag, terms, verbose = simona_opt$verbose) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	if(any(duplicated(id))) {
		stop("`term` should not be duplicated.")
	}
	d = exec_under_message_condition({
		cpp_longest_distances_directed(dag, id)
	}, verbose = verbose)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}


longest_distances_from_LCA = function(dag, id, verbose = simona_opt$verbose) {
	exec_under_message_condition({
		cpp_longest_distances_from_LCA(dag, id)
	}, verbose = verbose)
}
