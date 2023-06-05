

.path_validate_input = function(dag, from, to, weight) {
	if(length(weight)) {
		if(length(dag@lt_children_relations) > 1) {
			relation_levels = attr(dag@lt_children_relations, "levels")
			if(length(weight) != length(relation_levels)) {
				stop("Length of `weight` should be the same as the levels of relations.")
			}
			if(!is.null(names(weight))) {
				if(length(setdiff(relation_levels, names(weight)))) {
					stop("If `weight` has names, they should cover all relation levels.")
				}

				weight = unname(weight[relation_levels])
			}
		}
	}
	if(!is.numeric(from)) {
		term2ind = structure(seq_along(dag@terms), names = dag@terms)
		from = unname(term2ind[from])
		if(any(is.na(from))) {
			stop("All terms in `from` should exist in the DAG.")
		}
	}

	if(!is.null(to)) {
		if(!is.numeric(to)) {
			to = unname(term2ind[to])
			if(any(is.na(to))) {
				stop("All terms in `to` should exist in the DAG.")
			}
		}
	}

	return(list(from = from, to = to, weight = weight))
}

.path_length = function(dag, from, to = NULL, weight = numeric(0), type = -1L) {

	lt = .path_validate_input(dag, from, to, weight)
	from = lt$from
	to = lt$to
	weight = lt$weight

	if(is.null(to)) {
		cpp_find_path_length(dag, from, weight, type)
	} else {
		cpp_find_path_length_2(dag, from, to, weight, type)
	}
}

#' Shortest path in the directed DAG
#' 
#' @param dag A `ontology_DAG` object.
#' @param from A vector of the "from" terms. For `shortest_path()` and `longest_path()`, it should be a scalar.
#' @param to A vector of the "to" terms. For `shortest_path()` and `longest_path()`, it should be a scalar.
#' @param weight Internally used.
#' 
#' @return
#' `shortest_path_length()` and `longest_path_length()` return a numeric matrix
#' `shortest_path()` returns a vector of term indicies or term names.
#' @export
#' @example
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' shortest_path_length(dag, c("a", "b", "c", "d", "e", "f"))
#' shortest_path_length(dag, c("a", "b", "c"), c("d", "e", "f"))
#' longest_path_length(dag, "a", "c")
#' shortest_path(dag, "a", "c")
#' shortest_path(dag, "a", "c", in_labels = TRUE)
#' longest_path(dag, "a", "c", in_labels = TRUE)
shortest_path_length = function(dag, from, to = NULL, weight = numeric(0)) {
	.path_length(dag, from, to, weight, -1L)
}

#' @rdname shortest_path_length
#' @export
longest_path_length = function(dag, from, to = NULL, weight = numeric(0)) {
	.path_length(dag, from, to, weight, 1L)
}

#' @param in_labels Whether the path is represented in their integer indices or in their names?
#' @rdname shortest_path_length
#' @export
shortest_path = function(dag, from, to, weight = numeric(0), in_labels = FALSE) {

	if(length(from) != 1) {
		stop("`from` should be a scalar.")
	}
	if(length(to) != 1) {
		stop("`to` should be a scalar.")
	}

	lt = .path_validate_input(dag, from, to, weight)
	from = lt$from
	to = lt$to
	weight = lt$weight

	path = cpp_find_path_single(dag, from, to, weight, -1L)
	if(in_labels) {
		dag@terms[path]
	} else {
		path
	}
}

#' @rdname shortest_path
#' @export
longest_path = function(dag, from, to, weight = numeric(0), in_labels = FALSE) {

	if(length(from) != 1) {
		stop("`from` should be a scalar.")
	}
	if(length(to) != 1) {
		stop("`to` should be a scalar.")
	}

	lt = .path_validate_input(dag, from, to, weight)
	from = lt$from
	to = lt$to
	weight = lt$weight

	path = cpp_find_path_single(dag, from, to, weight, 1L)
	if(in_labels) {
		dag@terms[path]
	} else {
		path
	}
}


dag_distance_undirected = function(dag, terms) {
	id = term_to_node_id(dag, terms)
	dist = cpp_distance_undirected(dag, id)
	dimnames(dist) = list(dag@terms[id], dag@terms[id])

	dist
}

dag_longest_distance_undirected_via_LCA = function(dag, terms) {
	id = term_to_node_id(dag, terms)
	dist = cpp_longest_distance_undirected_via_LCA(dag, id, 1L)
	dimnames(dist) = list(dag@terms[id], dag@terms[id])

	dist
}
