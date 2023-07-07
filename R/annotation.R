


#' Number of annotated items
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names. If it is set, the returned vector will be subsetted to the terms that have been set here.
#' @param uniquify Whether to uniquify items that are annotated to the term? See **Details**.
#' @param use_cache Internally used.
#' 
#' @details
#' Due to the nature of the DAG, a parent term includes all annotated items of its child terms, and an ancestor term includes
#' all annotated items from its offspring recursively. In current tools, there are two different implementations to deal with
#' such recursive merging.
#' 
#' For a term `t`, denote `S_1`, `S_2`, ... as the sets of annotated items for its child 1, 2, ..., also denote `S_t` as the set
#' of items that are **directly** annotated to `t`. The first method takes the union of annotated items on `t` and all its child terms:
#' 
#' ```
#' n = length(union(S_t, S_1, S_2, ...))
#' ```
#' 
#' And the second method takes the sum of numbers of items on `t` and on all its child terms:
#' 
#' ```
#' n = sum(length(s_t) + length(S_1) + length(S_2) + ...)
#' ```
#' 
#' In `n_annotations()`, when `uniquify = TRUE`, the first method is used; and when `uniquify = FALSE`, the second method is used.
#' 
#' For some annotation sources, it is possible that an item is annotated to multiple terms, thus, the second method which simply
#' adds numbers of all its child terms may not be proper because an item may be counted duplicatedly, thus over-estimating `n`. The two methods
#' are identical only if an item is annotated to a unique term in the DAG.
#' 
#' We suggest to always set `uniquify = TRUE` (the default), and the scenario of `uniquify = FALSE` is only for the testing or benchmarking purpose. 
#' 
#' @returns An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' n_annotations(dag)
n_annotations = function(dag, terms = NULL, uniquify = simone_opt$anno_uniquify, use_cache = simone_opt$use_cache) {
	if(!uniquify && is.null(dag@term_env$n_annotations)) {
		use_cache = FALSE
	} else if(uniquify && is.null(dag@term_env$n_annotations_unique)) {
		use_cache = FALSE
	}
	if(!use_cache) {
		validate_dag_has_annotation(dag)
	
		n = cpp_n_annotations(dag, uniquify)
		
		attr(n, "N") = max(n)

		if(uniquify) {
			dag@term_env$n_annotations_unique = n
		} else {
			dag@term_env$n_annotations = n
		}
	}
	
	if(uniquify) {
		n = structure(dag@term_env$n_annotations_unique, names = dag@terms)
	} else {
		n = structure(dag@term_env$n_annotations, names = dag@terms)
	}	
	
	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms)
		n[i]
	} else {
		n
	}

}

validate_dag_has_annotation = function(dag) {
	if(length(dag@annotation$list) == 0) {
		stop("`annotation` should be set in `create_ontology_DAG()`.")
	}
}


# remove terms with zero annotation and the root term
validate_annotated_terms = function(dag, id) {
	n_anno = n_annotations(dag)[id]

	l1 = n_anno == 0
	if(any(l1)) {
		message("remove ", sum(l1), " terms with no annotation.")
		
		if(length(sum(!l1)) == 0) {
			stop("No term is lef.")
		}
	}
	l2 = id %in% dag@root
	if(any(l2)) {
		message("remove root term.")
		
		if(length(sum(!l1 & !l2)) == 0) {
			stop("No term is lef.")
		}
	}

	!l1 & !l2
}


#' Term-item associations
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param return Whether the returned object is a list or a matrix?
#' 
#' @details
#' If an item is annotated to a term, all this term's ancestor terms are also annotated.
#' 
#' @rdname annotation
#' @return
#' A list or a binary matrix showing annotation relations between terms and items.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' term_annotations(dag, letters[1:6])
#' annotated_terms(dag, c("t1", "t2", "t3"))
term_annotations = function(dag, terms, return = "list") {
	validate_dag_has_annotation(dag)

	id = term_to_node_id(dag, terms, strict = FALSE)

	m = cpp_get_term_annotations(dag, id)
	rownames(m) = dag@terms[id]
	colnames(m) = dag@annotation$names
	
	if(return == "matrix") {
		m
	} else {
		cn = colnames(m)
		apply(m, 1, function(x) {
			cn[x > 0]
		}, simplify = FALSE)
	}
}

#' @param anno A vector of annotated item names.
#' 
#' @rdname annotation
#' @export
annotated_terms = function(dag, anno, return = "list") {
	validate_dag_has_annotation(dag)

	id = which(dag@annotation$names %in% anno)
	if(length(id) == 0) {
		stop("Cannot find `anno` in DAG annotations.")
	}

	m = cpp_get_annotated_terms(dag, id)
	rownames(m) = dag@annotation$names[id]
	colnames(m) = dag@terms

	if(return == "matrix") {
		m
	} else {
		cn = colnames(m)
		apply(m, 1, function(x) {
			cn[x > 0]
		}, simplify = FALSE)
	}
}
