


#' Number of annotated items
#' 
#' @param dag A `ontology_DAG` object.
#' @param use_cache Internally used.
#' 
#' @details
#' `annotation` argument should be set in [create_ontology_DAG()].
#' @returns An integer vector.
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")
#' head(n_annotations(dag))
#' }
n_annotations = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$n_annotations)) {
		use_cache = FALSE
	} 
	if(!use_cache) {
		validate_dag_has_annotation(dag)
		n_all_anno = length(dag@annotation$names)
	
		n = cpp_n_annotations(dag)
		
		attr(n, "N") = n_all_anno

		dag@term_env$n_annotations = n
	}
	
	structure(dag@term_env$n_annotations, names = dag@terms)
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
		message("removed ", sum(l1), " terms with no annotation.")
		
		if(length(sum(!l1)) == 0) {
			stop("No term is lef.")
		}
	}
	l2 = id %in% dag@root
	if(any(l2)) {
		message("removed root term.")
		
		if(length(sum(!l1 & !l2)) == 0) {
			stop("No term is lef.")
		}
	}

	!l1 & !l2
}


#' Get term-item associations
#' 
#' @param dag A `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param return Whether to return as a list or a matrix?
#' 
#' @return A list of annotated items or a list of terms or a binary matrix showing whether the term is annotated with the item.
#' @export
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
#' @rdname term_annotations
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
