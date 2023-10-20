
#' Filter the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names. The sub-DAG will only contain these terms.
#' @param relations A vector of relations. The sub-DAG will only contain these relations. 
#'                  Valid values of "relations" should correspond to the values set in the 
#'                  `relations` argument in the [`create_ontology_DAG()`]. If `relations_DAG` is 
#'                  already provided, offspring relation types will all be selected. Note "is_a"
#'                  is always included.
#' @param root A vector of term names which will be used as roots of the sub-DAG. Only 
#'             these with their offspring terms will be kept. If there are multiple root terms, 
#'             a super root will be automatically added.
#' @param leaves A vector of leaf terms. Only these with their ancestor terms will be kept.
#' @param mcols_filter Filtering on columns in the meta data frame.
#' 
#' @details If the DAG is reduced into several disconnected parts after the filtering, a
#'          super root is automatically added.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_filter(dag, terms = c("b", "d", "f"))
#' dag_filter(dag, root = "b")
#' dag_filter(dag, leaves = c("c", "b"))
#' dag_filter(dag, root = "b", leaves = "e")
#' 
#' \donttest{
#' dag = create_ontology_DAG_from_GO_db()
#' dag_filter(dag, relations = "is_a")
#' }
dag_filter = function(dag, terms = NULL, relations = NULL, root = NULL, leaves = NULL,
	mcols_filter = NULL) {

	lt_children = dag@lt_children
	children = unlist(lt_children)
	parents = rep(seq_along(dag@terms), times = vapply(lt_children, length, FUN.VALUE = integer(1)))
	if(length(dag@lt_children_relations) > 0) {
		relation_levels = attr(dag@lt_children_relations, "levels")
		v_relations = unlist(dag@lt_children_relations)
		v_relations = relation_levels[v_relations]
	} else {
		v_relations = NULL
	}

	if(length(dag@annotation$list) > 0) {
		annotation = lapply(dag@annotation$list, function(x) {
			dag@annotation$names[x]
		})
		names(annotation) = dag@terms
	} else {
		annotation = NULL
	}

	l = rep(TRUE, length(parents))

	if(!is.null(terms)) {
		terms = term_to_node_id(dag, terms)
		l = l & parents %in% terms & children %in% terms
	}

	if(!is.null(relations)) {

		relations = normalize_relation_type(relations)
		if(!"is_a" %in% relations) {
			stop("'is_a' must be included.")
		}

		if(!is.null(dag@relations_DAG)) {
			relations = merge_offspring_relation_types(dag@relations_DAG, relations)
		}

		if(!is.null(v_relations)) {
			l2 = v_relations %in% relations
			if(!any(l2)) {
				stop("Cannot find any relation.")
			}
			l = l & l2
		}
	}

	if(!is.null(root)) {
		offspring = dag_offspring(dag, root, FALSE)
		offspring = c(offspring, which(dag@terms %in% root))
		l = l & (parents %in% offspring & children %in% offspring)
	}

	if(!is.null(leaves)) {
		ancestors = dag_ancestors(dag, leaves, FALSE)
		ancestors = c(ancestors, which(dag@terms %in% leaves))
		l = l & (parents %in% ancestors & children %in% ancestors)
	}

	## filter in mcols
	df = mcols(dag)
	if(!is.null(df)) {
		l2 = eval(substitute(mcols_filter), envir = df)
		if(!is.null(l2)) {
			l2[is.na(l2)] = FALSE
			if(!is.null(l2)) {
				l = l & l2
			}
		}
	}

	parents = parents[l]
	children = children[l]
	if(!is.null(v_relations)) {
		v_relations = v_relations[l]
	}

	parents = dag@terms[parents]
	children = dag@terms[children]

	all_terms = unique(c(parents, children))
	if(!is.null(annotation)) {
		annotation = annotation[all_terms]
	}

	dag2 = create_ontology_DAG(parents = parents, children = children,
		relations = v_relations, relations_DAG = dag@relations_DAG, 
		annotation = annotation, source = dag@source)

	meta = mcols(dag)
	if(!is.null(meta)) {
		mcols(dag2) = meta[dag2@terms, , drop = FALSE]
	}

	dag2
}



#' Create sub-DAGs
#' 
#' @param x An `ontology_DAG` object.
#' @param i A single term name. The value should be a character vector. It corresponds to the roots.
#' @param j A single term name. The value should be a character vector. It corresponds to the leaves.
#' @param ... Ignored.
#' @param drop Ignored.
#' 
#' @details It returns a sub-DAG taking node `i` as the root and `j` as the leaves. If `i` is a vector, a super root will be added.
#' @return An `ontology_DAG` object.
#' 
#' @rdname subset
#' @exportMethod [
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag["b"]
#' dag[["b"]]
#' dag["b", "f"]
#' dag[, "f"]
setMethod("[", 
	signature = c("ontology_DAG", i = "ANY", j = "ANY", drop = "missing"),
	definition = function(x, i, j, ..., drop = FALSE) {

	if(!is.character(i)) {
		stop("Only the character term name should be used as index.")
	}

	if(!is.character(j)) {
		stop("Only the character term name should be used as index.")
	}
	
	dag_filter(x, root = i, leaves = j)
	
})


#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "ANY", j = "ANY", drop = "ANY"),
	definition = function(x, i, j, ..., drop = FALSE) {

	x[i, j]
	
})


#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "ANY", j = "missing", drop = "missing"),
	definition = function(x, i, j, ..., drop = FALSE) {

	if(!is.character(i)) {
		stop("Only the character term name should be used as index.")
	}

	dag_filter(x, root = i)
})

#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "ANY", j = "missing", drop = "ANY"),
	definition = function(x, i, j, ..., drop = FALSE) {

	x[i]
})

#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "missing", j = "ANY", drop = "missing"),
	definition = function(x, i, j, ..., drop = FALSE) {

	if(!is.character(j)) {
		stop("Only the character term name should be used as index.")
	}

	dag_filter(x, leaves = j)

})

#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "missing", j = "ANY", drop = "ANY"),
	definition = function(x, i, j, ..., drop = FALSE) {

	x[, j]

})

#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "missing", j = "missing", drop = "missing"),
	definition = function(x, i, j, ..., drop = FALSE) {

	x
	
})

#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", i = "missing", j = "missing", drop = "ANY"),
	definition = function(x, i, j, ..., drop = FALSE) {

	x
	
})

#' @rdname subset
#' @exportMethod [[
setMethod("[[", 
	signature = c("ontology_DAG", "character", "missing"),
	definition = function(x, i, j, ...) {

	if(!is.character(i)) {
		stop("Only the character term name should be used as index.")
	}
	dag_filter(x, root = i)
})
