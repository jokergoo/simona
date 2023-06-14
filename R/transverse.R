

#' Parent/child/ancestor/offspring terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param term For `dag_parents()` and `dag_children()`, the value should be a single term name.
#'             For `dag_ancestors()` and `dag_offspring()`, the value can be a vector of multiple term names. If it is a vector, it returns
#'             union of the ancestor or offspring terms of the selected set of terms.
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


#' @rdname dag_terms
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


#' Number of parent/child/ancestor/offspring/leaf terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param term A vector of term names. If the value is `NULL`, it returns for all terms in the DAG.
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
n_offspring = function(dag, term = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_offspring) || !use_cache) {
		dag@term_env$n_offspring = cpp_n_offspring(dag, include_self)
	}
	n = dag@term_env$n_offspring
	names(n) = dag@terms

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_ancestors = function(dag, term = NULL, use_cache = TRUE, include_self = FALSE) {
	if(is.null(dag@term_env$n_ancestors) || !use_cache) {
		dag@term_env$n_ancestors = cpp_n_ancestors(dag, include_self)
	}
	n = dag@term_env$n_ancestors
	names(n) = dag@terms

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @details Leaf nodes have zero `n_leaves` value, so you can identify leaf terms based on the value.
#' @export
n_leaves = function(dag, term = NULL, use_cache = TRUE) {
	if(is.null(dag@term_env$n_leaves) || !use_cache) {
		dag@term_env$n_leaves = cpp_n_leaves(dag)
	}
	n = dag@term_env$n_leaves
	names(n) = dag@terms

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_parents = function(dag, term = NULL, use_cache = TRUE) {
	n = dag@term_env$n_parents
	names(n) = dag@terms

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
		n[i]
	} else {
		n
	}
}

#' @rdname n_terms
#' @export
n_children = function(dag, term = NULL, use_cache = TRUE) {
	n = dag@term_env$n_children
	names(n) = dag@terms

	if(!is.null(term)) {
		i = term_to_node_id(dag, term)
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
#' @param use_cache Internally used.
#' @details
#' The depth of a term in the DAG is defined as the maximal distance from the root. The height
#' of a term in the DAG is the maximal finite distance to the leaf terms
#' 
#' `dag_dist_from_root()` and `dag_dist_from_leaves()` calculate the minimal distance from the root or to the leaves.
#' 
#' @return An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_depth(dag)
#' dag_height(dag)
#' dag_dist_from_root(dag)
#' dag_dist_from_leaves(dag)
dag_depth = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_depth) || !use_cache) {
		d = cpp_dag_depth(dag)
		dag@term_env$dag_depth = d
	}
	structure(dag@term_env$dag_depth, names = dag@terms)
}

#' @rdname dag_depth
#' @export
dag_height = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_height) | !use_cache) {
		d = cpp_dag_height(dag)
		dag@term_env$dag_height = d
	}
	structure(dag@term_env$dag_height, names = dag@terms)
}

#' @rdname dag_depth
#' @export
dag_dist_from_root = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_root) || !use_cache) {
		d = cpp_dag_dist_from_root(dag)
		dag@term_env$dag_dist_to_root = d
	}
	structure(dag@term_env$dag_dist_to_root, names = dag@terms)
}

#' @rdname dag_depth
#' @export
dag_dist_from_leaves = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_dist_to_leaves) || !use_cache) {
		d = cpp_dag_dist_from_leaves(dag)
		dag@term_env$dag_dist_to_leaves = d
	}
	structure(dag@term_env$dag_dist_to_leaves, names = dag@terms)
}
