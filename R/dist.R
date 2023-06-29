
#' Distance on the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' 
#' @details
#' Denote two terms as `a` and `b`, a common ancestor as `c`, and the distance function `d()` calculates the longest
#' distance or the shortest distance depending on the function.
#' 
#' - `shortest_distances_via_NCA()`: It calculates the smallest `d(c, a) + d(c, b)` where `d()` calculates the shortest distance between two terms. In this case,
#'     `c` is the NCA (nearest common ancestor) of `a` and `b`.
#' - `longest_distances_via_LCA()`: It calculates the largest `d(c, a) + d(c, b)` where `d()` calculates the longest distance between two terms via the LCA term. In this case,
#'     `c` is the LCA of `a` and `b`.
#' - `shortest_distances_directed()`: It calculates `d(a, b)`. The distance is only calculated when `a` is an ancestor of `b`.
#' - `longest_distances_directed()`: It calculates `d(a, b)`. The distance is only calculated when `a` is an ancestor of `b`.
#' @rdname distance
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' shortest_distances_via_NCA(dag, letters[1:6])
#' longest_distances_via_LCA(dag, letters[1:6])
#' shortest_distances_directed(dag, letters[1:6])
#' longest_distances_directed(dag, letters[1:6])
shortest_distances_via_NCA = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	d = cpp_shortest_distances_via_NCA(dag, id)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
longest_distances_via_LCA = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	d = cpp_longest_distances_via_LCA(dag, id)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
shortest_distances_directed = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	d = cpp_shortest_distances_directed(dag, id)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname distance
#' @export
longest_distances_directed = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	d = cpp_longest_distances_directed(dag, id)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}
