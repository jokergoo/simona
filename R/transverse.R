

#' Terms' parents/children/ancestor/offspring terms
#' 
#' @param dag A `ontology_DAG` object.
#' @param term For `dag_parents()` and `dag_children()`, the value should be a single term name.
#'             For `dag_ancestor()` and `dag_offspring()`, the value can be a vector of term names.
#' @param in_labels Whether to return the numeric indices or in their term names.
#' 
#' @return An integer vector or a character vector depending on the value of `in_labels`.
#' @export
#' @example
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_parents(dag, "b")
#' dag_parents(dag, "c", in_labels = FALSE)
#' dag_children(dag, "b")
#' dag_ancestor(dag, "e")
#' dag_ancestor(dag, "b")
dag_parents = function(dag, term, in_labels = TRUE) {
	if(length(term) != 1) {
		stop("`term` should be a scalar.")
	}

	i = term_to_node_id(dag, term)

	parents = dag@lt_parents[[i]]
	if(in_labels) {
		dag@terms[parents]
	} else {
		parents
	}
}


#' @rdname dag_parents
#' @export
dag_children = function(dag, term, in_labels = TRUE) {
	if(length(term) != 1) {
		stop("`term` should be a scalar.")
	}

	i = term_to_node_id(dag, term)

	children = dag@lt_children[[i]]
	if(in_labels) {
		dag@terms[children]
	} else {
		children
	}
}

#' @rdname dag_parents
#' @export
dag_ancestor = function(dag, term, in_labels = TRUE, include_self = FALSE) {
	if(length(term) == 1) {
		i = term_to_node_id(dag, term)

		ancestor = cpp_ancestor(dag, i, include_self)
		if(in_labels) {
			dag@terms[ancestor]
		} else {
			ancestor
		}
	} else {
		i = term_to_node_id(dag, term)

		ancestor = cpp_ancestor_of_a_group(dag, i, include_self)
		if(in_labels) {
			dag@terms[ancestor]
		} else {
			ancestor
		}
	}
}

#' @rdname dag_parents
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


#' Number of parent/child/ancestor/offspring/leaf terms
#' 
#' @param dag A `ontology_DAG` object.
#' @param term A vector of term names. If the value is `NULL`, it returns for all terms.
#' @param use_cache Internally used.
#' 
#' @return An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' n_parents(dag)
#' n_children(dag)
#' n_offspring(dag)
#' n_ancestor(dag)
#' n_leaves(dag)
n_offspring = function(dag, term = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_offspring) || !use_cache) {
		dag@term_env$n_offspring = cpp_n_offspring(dag, include_self)
	}
	n = dag@term_env$n_offspring

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_offspring
#' @export
n_ancestor = function(dag, term = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_ancestor) || !use_cache) {
		dag@term_env$n_ancestor = cpp_n_ancestor(dag, include_self)
	}
	n = dag@term_env$n_ancestor

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_offspring
#' @details Leaf nodes have zero `n_leaves` value.
#' @export
n_leaves = function(dag, term = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$n_leaves) || !use_cache) {
		dag@term_env$n_leaves = cpp_n_leaves(dag)
	}
	n = dag@term_env$n_leaves

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_offspring
#' @export
n_parents = function(dag, term = NULL, use_cache = TRUE) {
	n = dag@term_env$n_parents

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_offspring
#' @export
n_children = function(dag, term = NULL, use_cache = TRUE) {
	n = dag@term_env$n_children

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}


dag_ancestor_of_two_groups = function(dag, group1, group2, type = "union", in_labels = FALSE, include_self = FALSE) {

	group1 = term_to_node_id(dag, group1)
	group2 = term_to_node_id(dag, group2)

	all_types = c("union" = 1, "intersect" = 2, "group 1" = 3, "group 2" = 4)
	type = all_types[tolower(type)]
	if(is.na(type)) {
		stop("Values of `type` can only be one of 'union', 'intersect', 'group 1' and 'group 2'.")
	}
	ancestor = cpp_ancestor_of_two_groups(dag, group1, group2, type, include_self)
	if(in_labels) {
		dag@terms[ancestor]
	} else {
		ancestor
	}
}



#' Depth and height of the DAG
#' 
#' @param dag A `ontology_DAG` object.
#' @param use_cache Internally used.
#' @details
#' The depth of a term in the DAG is defined as the maximal distance to the root. The height
#' of a term in the DAG is the maximal finite distance to the leaf nodes.
#' 
#' @return An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_depth(dag)
#' dag_height(dag)
dag_depth = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_depth) || !use_cache) {
		d = cpp_dag_depth(dag)
		dag@term_env$dag_depth = d
	}
	dag@term_env$dag_depth
}

#' @rdname dag_depth
#' @export
dag_height = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_height) | !use_cache) {
		d = cpp_dag_height(dag)
		dag@term_env$dag_height = d
	}
	dag@term_env$dag_height
}

#' @rdname dag_depth
#' @export
dag_dist_to_root = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_root) || !use_cache) {
		d = cpp_dag_dist_to_root(dag)
		dag@term_env$dag_dist_to_root = d
	}
	dag@term_env$dag_dist_to_root
}

#' @rdname dag_depth
#' @export
dag_dist_to_leaves = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_leaves) || !use_cache) {
		d = cpp_dag_dist_to_leaves(dag)
		dag@term_env$dag_dist_to_leaves = d
	}
	dag@term_env$dag_dist_to_leaves
}
