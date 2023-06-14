

#' The ontology_DAG class
#' 
#' This class defines the DAG structure of an ontology.
#' 
#' @slot lt_parents A list of length `n`. Each element is an integer index vector of parent terms.
#' @slot lt_children foo
#' @slot lt_children_relations foo
#' @slot source The source of the ontology. A character scalar only used as a mark of the returned object.
#' @slot terms A character vector of length `n`, which contains the term names.
#' @slot n_terms An integer scalar of the total number of terms in the DAG.
#' @slot root An integer index of the root term.
#' @slot leaves An integer vector of the indicies of leaf terms.
#' @slot tpl_sorted An integer vector of reordered term indices which has been topologically sorted in the DAG.
#' @slot tpl_pos The position of the original term in the topologically sorted path (similar as the rank), e.g. the position of term 1 in the sorted path.
#' @slot annotation A list of two elements: `list` and `names`. The `list` element contains a list of length `n` and each element
#'         there is a vector of integer indices of annotated items. The full list of annotated items is in the `names` element.
#' @slot term_env An environment which contains various term-level statistics, mainly for cache purpose.
#' 
#' @export ontology_DAG
#' @exportClass ontology_DAG
ontology_DAG = setClass("ontology_DAG",
	slots = c("lt_parents" = "list",
		      "lt_children" = "list",
		      "lt_children_relations" = "list",
		      "source" = "character",
		      "terms" = "character",
		      "n_terms" = "integer",
		      "root" = "integer",
		      "leaves" = "integer",
		      "tpl_sorted" = "integer",
		      "tpl_pos" = "integer",
		      "annotation" = "list",
		      "term_env" = "environment"
		      )
)

#' Create the ontology_DAG object
#' 
#' @param parents A character vector of parent terms. 
#' @param children A character vector of child terms. 
#' @param relations A character vector of parent-child relations, e.g. "isa" or "part of".
#' @param source Source of the ontology. It is simply used as a mark of the object.
#' @param annotation A list of character items that are annotated to the terms. Names of the list should be the term names. In the DAG, items
#'                   annotated to a term will also be annotated to its parents. Such upstreaming merging
#'                   is applied automatically in the package.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom methods new
#' @import Rcpp
#' @useDynLib ontsim, .registration = TRUE
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' 
#' # with annotations
#' annotation = list(
#'     "a" = 1:3,
#'     "b" = 3:4,
#'     "c" = 5,
#'     "d" = 7,
#'     "e" = 4:7,
#'     "f" = 8
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' 
#' # with relations
#' dag = create_ontology_DAG(parents, children, 
#'     relations = c("isa", "part of", "isa", "part of", "isa", "part of"))
create_ontology_DAG = function(parents, children, relations = NULL, 
	source = "Ontology", annotation = NULL) {

	if(!is.character(parents)) {
		stop("`parents` should be in character mode.")
	}

	if(!is.character(children)) {
		stop("`children` should be in character mode.")
	}

	if(length(parents) == 0 || length(children) == 0) {
		stop("There is no relation.")
	}

	has_relations = !is.null(relations)
	lt_relations = list()

	terms = sort(unique(c(parents, children)))
	n_terms = length(terms)
	terms2ind = structure(seq_along(terms), names = terms)

	lt_parents2 = split(unname(terms2ind[parents]), children)
	lt_children2 = split(unname(terms2ind[children]), parents)
	if(has_relations) {
		relations = as.factor(relations)
		relation_levels = levels(relations)
		relations = as.integer(relations)
		lt_relations2 = split(relations, parents)
	}

	lt_parents = rep(list(integer(0)), n_terms)
	lt_children = lt_parents
	lt_parents[ terms2ind[names(lt_parents2)] ] = lt_parents2
	lt_children[ terms2ind[names(lt_children2)] ] = lt_children2
	if(has_relations) {
		lt_relations = rep(list(integer(0)), n_terms)
		lt_relations[ terms2ind[names(lt_relations2)] ] = lt_relations2
		attr(lt_relations, "levels") = relation_levels
	}

	root = which(sapply(lt_parents, length) == 0)
	leaves = which(sapply(lt_children, length) == 0)

	if(length(root) == 0) {
		stop("Cannot find the root. There might exist a cycle.")
	}
	
	if(length(root) > 1) {

		warning("There are more than one root:\n", "  ", paste(terms[root], collapse = ", "), "\nA super root (_all_) is added.")
		
		super_root = n_terms + 1L
		lt_parents[[super_root]] = integer(0)
		for(r in root) {
			lt_parents[[r]] = super_root
		}

		lt_children[[super_root]] = root
		if(has_relations) lt_relations[[super_root]] = rep(which(relation_levels == "isa"), length(root))
		terms = c(terms, "_all_")
		root = super_root
		n_terms = n_terms + 1
	}
		
	dag = ontology_DAG(
		lt_parents = lt_parents,
		lt_children = lt_children,
		lt_children_relations = lt_relations, 
		source = source,
		terms = terms,
		n_terms = length(terms),
		root = root,
		leaves = leaves,
		term_env = new.env(parent = emptyenv())
	)

	cpp_check_cyclic_node(dag)

	dag@term_env$n_parents = sapply(lt_parents, length)
	dag@term_env$n_children = sapply(lt_children, length)

	tpl_order = order(dag_depth(dag), dag@term_env$n_children, dag@term_env$n_parents, dag@terms)
	dag@tpl_sorted = seq_len(n_terms)[tpl_order]
	dag@tpl_pos = seq_len(n_terms)
	
	# pos has the same order as terms, but the value are the positions on `sorted`, i.e. the rank
	dag@tpl_pos[tpl_order] = seq_len(n_terms)

	if(!is.null(annotation)) {
		annotation = lapply(annotation, as.character)
		all_anno = unique(unlist(annotation))
		n_all_anno = length(all_anno)
		anno2ind = structure(seq_along(all_anno), names = all_anno)
		annotation = lapply(annotation, function(x) unname(anno2ind[as.character(x)]))

		n_terms = dag@n_terms
		annotation2 = rep(list(integer(0)), n_terms)
		term2ind = structure(seq_along(dag@terms), names = dag@terms)

		cn = intersect(dag@terms, names(annotation))
		if(length(cn) == 0) {
			stop("No overlap between annotation names and terms.")
		}
		annotation2[ term2ind[cn] ] = annotation[cn]

		dag@annotation = list(list = annotation2,
			                  names = all_anno)
	} else {
		dag@annotation = list(list = vector("list", 0), names = character(0))
	}

	return(dag)
}

#' Print the ontology_DAG object
#' 
#' @param object An `ontology_DAG` object.
#' @exportMethod show
#' @importFrom stats quantile
setMethod("show",
	signature = "ontology_DAG",
	definition = function(object) {

		n_terms = object@n_terms
		n_relations = sum(sapply(object@lt_children, length))

		cat("An ontology_DAG object:\n")
		cat("  Source:", object@source, "\n")
		cat("  ", n_terms, " terms / ", n_relations, " relations\n", sep = "")
		cat("  Root:", paste(object@terms[object@root], collapse = ", "), "\n")

		if(!is.null(object@term_env$dag_depth)) {
			depth = object@term_env$dag_depth
			cat("  Max depth:", max(depth), "\n")
		}

		depth = dag_depth(object)
		tb1 = table(depth)
		aspect_ratio1 = max(tb1)/quantile(depth, 0.99)
		aspect_ratio1 = round(aspect_ratio1, 2)

		dist = cpp_dag_dist_from_root(object)
		tb2 = table(dist)
		aspect_ratio2 = max(tb2)/quantile(dist, 0.99)
		aspect_ratio2 = round(aspect_ratio2, 2)

		cat("  Aspect ratio: ", aspect_ratio1, ":1 (based on the longest_dist to root)\n", sep = "")
		cat("                ", aspect_ratio2, ":1 (based on the shortest_dist to root)\n", sep = "")

		if(length(object@lt_children_relations)) {
			cat("  Relations:", paste(attr(object@lt_children_relations, "levels"), collapse = ", "), "\n")
		}

		if(length(object@annotation$list)) {
			cat("  Annotations are available.\n")
		}
	}
)

#' Create the ontology_DAG object from the GO.db package
#' 
#' @param namespace One of "BP", "CC" and "MF".
#' @param relations Type of the GO term relations.
#' @param org_db The name of the organism package or the corresponding database object, e.g. `"org.Hs.eg.db"` or 
#'            directly the [org.Hs.eg.db::org.Hs.eg.db] object for human. Then the gene annotation to GO terms will be added
#'            to the object.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db()
#' dag = create_ontology_DAG_from_GO_db("BP", 
#'     relations = "isa", org_db = "org.Hs.eg.db")
#' }
create_ontology_DAG_from_GO_db = function(namespace = "BP", relations = c("isa", "part of"), org_db = NULL) {

	check_pkg("GO.db", bioc = TRUE)

	if(namespace == "BP") {
		df = AnnotationDbi::toTable(GO.db::GOBPCHILDREN)
	} else if(namespace == "CC") {
		df = AnnotationDbi::toTable(GO.db::GOCCCHILDREN)
	} else if(namespace == "MF") {
		df = AnnotationDbi::toTable(GO.db::GOMFCHILDREN)
	} else {
		stop("Value of `namespace` can only be one of 'BP', 'MF' and 'CC'.")
	}

	df = df[df[, 3] %in% relations, , drop = FALSE]

	if(!is.null(org_db)) {
		if(is.character(org_db)) {
			check_pkg(org_db, bioc = TRUE)
			org_db = getFromNamespace(org_db, ns = org_db)
		}
	
		suppressMessages(tb <- AnnotationDbi::select(org_db, keys = AnnotationDbi::keys(org_db), columns = c("GO", "ONTOLOGY")))
		tb = tb[tb$ONTOLOGY == namespace, , drop = FALSE]
		tb = tb[, c(1, which(colnames(tb) == "GO")), drop = FALSE]
		annotation = split(tb[, 1], tb[, 2])
		annotation = lapply(annotation, unique)
		annotation = lapply(annotation, as.character)
	} else {
		annotation = NULL
	}

	create_ontology_DAG(parents = df[, 2], children = df[, 1], relations = df[, 3],
		annotation = annotation, source = paste0("GO ", namespace))
}


#' Create a sub-DAG
#' 
#' @param x An `ontology_DAG` object.
#' @param i A single term name.
#' @param j Ignored.
#' @param ... Ignored.
#' @param drop Ignored.
#' 
#' @details It returns a sub-DAG taking node `i` as the root.
#' @return An `ontology_DAG` object.
#' 
#' @rdname subset
#' @exportMethod [
setMethod("[", 
	signature = c("ontology_DAG", "character", "missing", "ANY"),
	definition = function(x, i, j, ..., drop = FALSE) {

	dag_filter(x, root = i)
})

#' @rdname subset
#' @exportMethod [[
setMethod("[[", 
	signature = c("ontology_DAG", "character", "missing"),
	definition = function(x, i, j, ...) {

	dag_filter(x, root = i)
})

#' Names of all terms
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @return A vector of term names.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_all_terms(dag)
dag_all_terms = function(dag) {
	dag@terms
}

#' Root or leaves of the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param in_labels Whether the terms are in their names or as integer indices?
#' 
#' @return A character or an integer vector.
#' @export
dag_root = function(dag, in_labels = TRUE) {
	if(in_labels) {
		dag@terms[dag@root]
	} else {
		dag@root
	}
}

#' @rdname dag_root
#' @export
dag_leaves = function(dag, in_labels = TRUE) {
	if(in_labels) {
		dag@terms[dag@leaves]
	} else {
		dag@leaves
	}
}

#' Convert to an igraph object
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @details if `relations` is set in `create_ontology_DAG()`, relations are also set as an edge attribute of the `igraph` object.
#' 
#' @return An `igraph` object.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_as_igraph(dag)
dag_as_igraph = function(dag) {

	check_pkg("igraph")

	lt_children = dag@lt_children
	children = unlist(lt_children)
	parents = rep(seq_along(dag@terms), times = sapply(lt_children, length))
	if(length(dag@lt_children_relations) > 0) {
		relation_levels = attr(dag@lt_children_relations, "levels")
		relations = unlist(dag@lt_children_relations)
		relations = relation_levels[relations]
	} else {
		relations = NULL
	}

	children = dag@terms[children]
	parents = dag@terms[parents]

	g = igraph::graph_from_edgelist(cbind(parents, children))
	if(!is.null(relations)) {
		igraph::set_edge_attr(g, "relation", value = relations)
	}

	g
}

#' Filter the DAG
#' 
#' @param dag A `ontology_DAG` object.
#' @param terms A vectof of term names. The sub DAG will only contain these terms.
#' @param relations A vector of relations. The sub DAG will only contain these relations.
#' @param root A vector of root terms. Only they with their offspring terms will be kept.
#' @param leaves A vector of leaf terms. Only they with their ancestor terms will are kept.
#' 
#' @details If the DAG is reduced into several disconnected parts after the term filtering, a
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
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db()
#' dag_filter(dag, relations = "isa")
#' }
dag_filter = function(dag, terms = NULL, relations = NULL, root = NULL, leaves = NULL) {

	lt_children = dag@lt_children
	children = unlist(lt_children)
	parents = rep(seq_along(dag@terms), times = sapply(lt_children, length))
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

	l = rep(TRUE, dag@n_terms)

	if(!is.null(terms)) {
		terms = term_to_node_id(dag, terms)
		l = l & parents %in% terms & children %in% terms
	}

	if(!is.null(relations)) {
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

	parents = parents[l]
	children = children[l]
	if(!is.null(relations)) {
		relations = relations[l]
	}
	if(!is.null(annotation)) {
		annotation = annotation[l]
	}
	parents = dag@terms[parents]
	children = dag@terms[children]

	create_ontology_DAG(parents = parents, children = children,
		relations = relations, annotation = annotation, source = dag@source)
}

