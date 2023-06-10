


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


#' Most informative common ancestor (MICA) and lowest common ancestor (LCA)
#' 
#' @param dag A `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param IC_method An IC method. Valid values are in [ALL_IC_METHODS].
#' 
#' @return 
#' `MICA_term()` returns a character matrix of the MICA terms. 
#' `MICA_IC()` returns a numeric matrix of the IC of the MICA terms.
#' `LCA_term()` returns a character matrix of the LCA term.
#' `LCA_depth()` reutrns an integer matrix of the depth of the LCA terms.
MICA_term = function(dag, terms, IC_method, in_labels = TRUE) {
	ic = term_IC(dag, IC_method)
	max_ancestor_id(dag, terms, ic, in_labels = in_labels)
}

#' @rdname MICA_term
MICA_IC = function(dag, terms, IC_method) {
	ic = term_IC(dag, IC_method)
	max_ancestor_v(dag, terms, ic)
}


#' @rdname MICA_term
LCA_term = function(dag, terms, in_labels = TRUE) {
	depth = dag_depth(dag)
	
	max_ancestor_id(dag, terms, depth, in_labels = in_labels)
}

#' @rdname MICA_term
LCA_depth = function(dag, terms) {
	depth = dag_depth(dag)
	
	max_ancestor_v(dag, terms, depth)
}

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


#' Distance on the DAG
#' 
#' @param dag A `ontology_DAG` object.
#' @param terms A vector of term names.
#' 
#' @details
#' - `shortest_distances_via_CA()`: it is a sum of two distances: d(c, a) + d(c, b)
#' where c is a common ancestor of a and b and d(c, a) is the shortest distance from c to a.
#' - `longest_distances_via_LCA()`: here c is the lowest common ancestor (with the largest depth) of a and b,
#' and d(c, a) is the longest distance from c to a.
#' - `shortest_distances_directed()`: it calculates d(a, b). The distance is only calculated when a is an ancestor of b.
#' - `longest_distances_directed()`: it calculates d(a, b). The distance is only calculated when a is an ancestor of b.
shortest_distances_via_CA = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	d = cpp_shortest_distances_via_CA(dag, id)

	dimnames(d) = list(dag@terms[id], dag@terms[id])
	d
}

#' @rdname shortest_distances_via_CA
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

#' @rdname shortest_distances_via_CA
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

#' @rdname shortest_distances_via_CA
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



