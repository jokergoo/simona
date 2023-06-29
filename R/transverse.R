

#' Parent/child/ancestor/offspring terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param term The value can be a vector of multiple term names. If it is a vector, it returns
#'             union of the upstream/downstream terms of the selected set of terms.
#' @param in_labels Whether the terms are represented in their names or as the integer indices?
#' @param include_self For `dag_offspring()` and `dag_ancestors()`, this controls whether to also contain the query term itself/themselves.
#' 
#' @return An integer vector or a character vector depending on the value of `in_labels`.
#' @rdname dag_terms
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_parents(dag, "b")
#' dag_parents(dag, "c", in_labels = FALSE)
#' dag_children(dag, "b")
#' dag_ancestors(dag, "e")
#' dag_ancestors(dag, "b")
dag_parents = function(dag, term, in_labels = TRUE) {

	i = term_to_node_id(dag, term)

	parents = unique(unlist(dag@lt_parents[i]))
	if(in_labels) {
		dag@terms[parents]
	} else {
		parents
	}
}


#' @rdname dag_terms
#' @export
dag_children = function(dag, term, in_labels = TRUE) {

	i = term_to_node_id(dag, term)

	children = unique(unlist(dag@lt_children[i]))
	if(in_labels) {
		dag@terms[children]
	} else {
		children
	}
}

#' @rdname dag_terms
#' @export
dag_ancestors = function(dag, term, in_labels = TRUE, include_self = FALSE) {
	if(length(term) == 1) {
		i = term_to_node_id(dag, term)

		ancestors = cpp_ancestors(dag, i, include_self)
		if(in_labels) {
			dag@terms[ancestors]
		} else {
			ancestors
		}
	} else {
		i = term_to_node_id(dag, term)

		ancestors = cpp_ancestors_of_a_group(dag, i, include_self)
		if(in_labels) {
			dag@terms[ancestors]
		} else {
			ancestors
		}
	}
}

#' @rdname dag_terms
#' @export
dag_offspring = function(dag, term, in_labels = TRUE, include_self = FALSE) {
	if(length(term) == 1) {
		i = term_to_node_id(dag, term)

		offspring = cpp_offspring(dag, i, include_self)
		if(in_labels) {
			dag@terms[offspring]
		} else {
			offspring
		}
	} else {
		i = term_to_node_id(dag, term)

		offspring = cpp_offspring_of_a_group(dag, i, include_self)
		if(in_labels) {
			dag@terms[offspring]
		} else {
			offspring
		}
	}
}


dag_all_offspring = function(dag, in_labels = TRUE, include_self = FALSE) {
	m = cpp_all_offspring(dag, include_self)
	lt = apply(m, 1, which)
	names(lt) = dag@terms

	if(in_labels) {
		lt = lapply(lt, function(x) dag@terms[x])
	}

	lt
}

#' Number of parent/child/ancestor/offspring/leaf terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names. If the value is `NULL`, it returns for all terms in the DAG.
#' @param use_cache Internally used.
#' @param include_self For `n_offspring()` and `n_ancestors()`, this controls whether to also contain the query term itself.
#' 
#' @return An integer vector.
#' @rdname n_terms
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' n_parents(dag)
#' n_children(dag)
#' n_offspring(dag)
#' n_ancestors(dag)
#' n_leaves(dag)
n_offspring = function(dag, terms = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_offspring) || !use_cache) {
		dag@term_env$n_offspring = cpp_n_offspring(dag, FALSE)
	}
	n = dag@term_env$n_offspring
	names(n) = dag@terms

	if(include_self) {
		n = n + 1
	}

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_ancestors = function(dag, terms = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_ancestors) || !use_cache) {
		dag@term_env$n_ancestors = cpp_n_ancestors(dag, FALSE)
	}
	n = dag@term_env$n_ancestors
	names(n) = dag@terms

	if(include_self) {
		n = n + 1
	}

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @details Leaf nodes have value of zero, so you can identify leaf terms based on the value.
#' @export
n_leaves = function(dag, terms = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$n_leaves) || !use_cache) {
		dag@term_env$n_leaves = cpp_n_leaves(dag)
	}
	n = dag@term_env$n_leaves
	names(n) = dag@terms

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_parents = function(dag, terms = NULL, use_cache = TRUE) {
	n = dag@term_env$n_parents
	names(n) = dag@terms

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_children = function(dag, terms = NULL, use_cache = TRUE) {
	n = dag@term_env$n_children
	names(n) = dag@terms

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}
}

dag_ancestors_of_two_groups = function(dag, group1, group2, type = "union", in_labels = FALSE, include_self = FALSE) {

	group1 = term_to_node_id(dag, group1)
	group2 = term_to_node_id(dag, group2)

	all_types = c("union" = 1, "intersect" = 2, "group 1" = 3, "group 2" = 4)
	type = all_types[tolower(type)]
	if(is.na(type)) {
		stop("Values of `type` can only be one of 'union', 'intersect', 'group 1' and 'group 2'.")
	}
	ancestors = cpp_ancestors_of_two_groups(dag, group1, group2, type, include_self)
	if(in_labels) {
		dag@terms[ancestors]
	} else {
		ancestors
	}
}



#' Depth and height in the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names. If it is set, the returned vector will be subsetted to the terms that have been set here.
#' @param use_cache Internally used.
#' @details
#' The depth of a term in the DAG is defined as the maximal distance from the root. The height
#' of a term in the DAG is the maximal finite distance to the leaf terms
#' 
#' `dag_shortest_dist_from_root()` and `dag_shortest_dist_to_leaves()` calculate the minimal distance from the root or to the leaves.
#' 
#' @return An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_depth(dag)
#' dag_height(dag)
#' dag_shortest_dist_from_root(dag)
#' dag_shortest_dist_to_leaves(dag)
dag_depth = function(dag, terms = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_depth) || !use_cache) {
		d = cpp_dag_depth(dag)
		dag@term_env$dag_depth = d
	}
	d = structure(dag@term_env$dag_depth, names = dag@terms)

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		d[i]
	} else {
		d
	}
}

#' @rdname dag_depth
#' @export
dag_height = function(dag, terms = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_height) | !use_cache) {
		d = cpp_dag_height(dag)
		dag@term_env$dag_height = d
	}
	d = structure(dag@term_env$dag_height, names = dag@terms)

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		d[i]
	} else {
		d
	}
}

#' @rdname dag_depth
#' @export
dag_shortest_dist_from_root = function(dag, terms = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_root) || !use_cache) {
		d = cpp_dag_dist_from_root(dag)
		dag@term_env$dag_dist_to_root = d
	}
	d = structure(dag@term_env$dag_dist_to_root, names = dag@terms)

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		d[i]
	} else {
		d
	}
}

#' @rdname dag_depth
#' @export
dag_shortest_dist_to_leaves = function(dag, terms = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_leaves) || !use_cache) {
		d = cpp_dag_dist_to_leaves(dag)
		dag@term_env$dag_dist_to_leaves = d
	}
	d = structure(dag@term_env$dag_dist_to_leaves, names = dag@terms)

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		d[i]
	} else {
		d
	}
}



.dag_singular_dist = function(dag, from, terms = NULL, background = NULL, dist_type = "longest", towards = "ancestors") {

	if(dist_type == "longest" && towards == "ancestors") {
		f = cpp_dag_longest_dist_from_ancestors
	} else if(dist_type == "shortest" && towards == "ancestors") {
		f = cpp_dag_shortest_dist_from_ancestors
	} else if(dist_type == "longest" && towards == "offspring") {
		f = cpp_dag_longest_dist_to_offspring
	} else if(dist_type == "shortest" && towards == "offspring") {
		f = cpp_dag_shortest_dist_to_offspring
	}
	
	from = term_to_node_id(dag, from)
	if(is.null(background)) {
		d = f(dag, from)
	} else {
		background = term_to_node_id(dag, background, strict = FALSE)
		l_background = rep(FALSE, dag@n_terms)
		l_background[background] = TRUE
		d = f(dag, from, l_background)
	}
	
	names(d) = dag@terms

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		d[i]
	} else {
		d
	}
}

#' Distance from terms in the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param from A single term or a vector of terms. For all downstream terms, the maximal or the minimal distance to the terms specified in this variable will be taken.
#' @param to Same format as the `from` argument.
#' @param terms A vector of term names. If it is set, the returned vector will be subsetted to the terms that have been set here.
#' @param background A vector of terms. Then the lookup will only be applied in this set of terms.
#' 
#' @details
#' As a special case, when `from` is the root term, `dag_longest_dist_to_offspring()` is the same as `dag_depth()`,
#' and when `to` are all leaf terms, `dag_shortest_dist_to_offspring()` is the same as `dag_height()`.
#' 
#' @export
dag_longest_dist_to_offspring = function(dag, from, terms = NULL, background = NULL) {
	.dag_singular_dist(dag, from, terms, background, dist_type = "longest", towards = "offspring")
}

#' @rdname dag_longest_dist_to_offspring
#' @export
dag_shortest_dist_to_offspring = function(dag, from, terms = NULL, background = NULL) {
	.dag_singular_dist(dag, from, terms, background, dist_type = "shortest", towards = "offspring")
}

#' @rdname dag_longest_dist_to_offspring
#' @export
dag_longest_dist_from_ancestors = function(dag, to, terms = NULL, background = NULL) {
	.dag_singular_dist(dag, to, terms, background, dist_type = "longest", towards = "ancestors")
}

#' @rdname dag_longest_dist_to_offspring
#' @export
dag_shortest_dist_from_ancestors = function(dag, to, terms = NULL, background = NULL) {
	.dag_singular_dist(dag, to, terms, background, dist_type = "shortest", towards = "ancestors")
}

