


#' The ontology_DAG class
#' 
#' This class defines the DAG structure of an ontology.
#' 
#' @slot terms A character vector of length `n` of all term names. Other slots that store term-level information use the integer indices of terms.
#' @slot n_terms An integer scalar of the total number of terms in the DAG.
#' @slot lt_parents A list of length `n`. Each element in the list is an integer index vector of the parent terms of the i^th term.
#' @slot lt_children A list of length `n`. Each element in the list is an integer index vector of the child terms of the i^th term.
#' @slot lt_children_relations A list of length `n`. Each element is a vector of the semantic relations between the i^th term and its child terms.
#'       The relations are represented as integers. The character name of the relations is in `attr(dag@lt_children_relations, "levels")`.
#' @slot relations_DAG A simple `ontology_DAG` object but constructed for relation types.
#' @slot source The source of the ontology. A character scalar only used as a mark of the returned object.
#' @slot root An integer scalar of the root term.
#' @slot leaves An integer vector of the indicies of leaf terms.
#' @slot tpl_sorted An integer vector of reordered term indices which has been topologically sorted in the DAG. Terms are sorted first by the depth (maximal
#'       distance from root), then the number of child terms, then the number of parent terms, and last the term names.
#' @slot tpl_pos The position of the original term in the topologically sorted path (similar as the rank), e.g. the value of the first element in the vector
#'       is the position of term 1 in the topologically sorted path.
#' @slot annotation A list of two elements: `list` and `names`. The `dag@annotation$list` element contains a list of length `n` and each element
#'       is a vector of integer indices of annotated items. The full list of annotated items is in `dag@annotation$names`.
#' @slot term_env An environment which contains various term-level statistics. It is mainly for cache purpose.
#' @slot aspect_ratio A numeric vector of length two. The aspect ratio is calculated as `w/h`. For each term, there is a distance to root, 
#'       `h` is the maximal distance of all terms, `w` is the maximal number of items with the same distance. The two values in the `aspect_ratio` slot
#'       use maximal distance to root (the height) and the shortest distance to root as the distance measure.
#' @slot elementMetadata An additional data frame with the same number of rows as the number of terms in DAG. Order of rows should be the same as order of terms in `dag@terms`.
#' 
#' @export ontology_DAG
#' @exportClass ontology_DAG
ontology_DAG = setClass("ontology_DAG",
	slots = c("terms" = "character",
		      "n_terms" = "integer",
		      "lt_parents" = "list",
		      "lt_children" = "list",
		      "lt_children_relations" = "list",
		      "relations_DAG" = "ANY",
		      "source" = "character",
		      "root" = "integer",
		      "leaves" = "integer",
		      "tpl_sorted" = "integer",
		      "tpl_pos" = "integer",
		      "annotation" = "list",
		      "term_env" = "environment",
		      "aspect_ratio" = "numeric",
		      "elementMetadata" = "ANY"
		      )
)

#' Create the ontology_DAG object
#' 
#' @param parents A character vector of parent terms. 
#' @param children A character vector of child terms. 
#' @param relations A character vector of parent-child relations, e.g. "is_a", "part_of", or self-defined semantic relations.
#' @param relations_DAG If the relation types also have hierarchical relations, it can also be constructed by `create_ontology_DAG()` first. See **Examples**.
#'          When the DAG for relation types is provided, in downstream analysis, the ancestor/offspring relationship for relation types will be taken into consideration.
#' @param source Source of the ontology. It is only used as a mark of the object.
#' @param annotation A list of character vectors which contain items annotated to the terms. Names of the list should be the term names. In the DAG, items
#'                   annotated to a term will also be annotated to its parents. Such merging
#'                   is applied automatically in the package.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom methods new
#' @importFrom stats quantile
#' @import Rcpp
#' @useDynLib simone, .registration = TRUE
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' 
#' # with annotations
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' 
#' # with relations
#' dag = create_ontology_DAG(parents, children, 
#'     relations = c("is_a", "part_of", "is_a", "part_of", "is_a", "part_of"))
#' 
#' # with relations_DAG
#' relations_DAG = create_ontology_DAG(c("r2", "r2"), c("r3", "r4"))
#' dag = create_ontology_DAG(parents, children, 
#'     relations = c("r1", "r2", "r1", "r3", "r1", "r4"),
#'     relations_DAG = relations_DAG)
create_ontology_DAG = function(parents, children, relations = NULL, relations_DAG = NULL,
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
		relations[grepl("^is.*a$", relations, ignore.case = TRUE)] = "is_a"
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
		txt = strwrap(paste(terms[root], collapse = ", "), width = 60)
		txt = paste(paste0("  ", txt), collapse = "\n")
		message("There are more than one root:\n", txt, "\n  A super root (_all_) is added.")
		
		super_root = n_terms + 1L
		lt_parents[[super_root]] = integer(0)
		for(r in root) {
			lt_parents[[r]] = super_root
		}

		lt_children[[super_root]] = root
		if(has_relations) lt_relations[[super_root]] = rep(which(relation_levels == "is_a"), length(root))
		terms = c(terms, "_all_")
		root = super_root
		n_terms = n_terms + 1
	}

	if(!is.null(relations_DAG)) {
		if(!inherits(relations_DAG, "ontology_DAG")) {
			stop("`relations_DAG` should be constructed by `create_ontology_DAG()`.")
		}
	}
		
	dag = ontology_DAG(
		terms = terms,
		n_terms = length(terms),
		lt_parents = lt_parents,
		lt_children = lt_children,
		lt_children_relations = lt_relations, 
		relations_DAG = relations_DAG,
		source = source,
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

	depth = dag_depth(dag)
	tb1 = table(depth)
	aspect_ratio1 = max(tb1)/quantile(depth, 0.99)
	aspect_ratio1 = round(aspect_ratio1, 2)

	dist = cpp_dag_dist_from_root(dag)
	tb2 = table(dist)
	aspect_ratio2 = max(tb2)/quantile(dist, 0.99)
	aspect_ratio2 = round(aspect_ratio2, 2)

	dag@aspect_ratio = c(aspect_ratio1, aspect_ratio2)

	return(dag)
}

#' Print the ontology_DAG object
#' 
#' @param object An `ontology_DAG` object.
#' @exportMethod show
setMethod("show",
	signature = "ontology_DAG",
	definition = function(object) {

		n_terms = object@n_terms
		n_relations = sum(sapply(object@lt_children, length))

		cat("An ontology_DAG object:\n")
		cat("  Source:", object@source, "\n")
		if(n_terms == n_relations + 1) {
			cat("  ", n_terms, " terms / ", n_relations, " relations / a tree\n", sep = "")
		} else {
			cat("  ", n_terms, " terms / ", n_relations, " relations\n", sep = "")
		}

		cat("  Root:", paste(object@terms[object@root], collapse = ", "), "\n")
		txt = strwrap(paste(c(object@terms[1:min(4, object@n_terms)], "..."), collapse = ", "), width = 60)
		txt[1] = paste0("  Terms: ", txt[1])
		txt[-1] = paste0("         ", txt[-1])
		cat(txt, sep = "\n")
		
		if(!is.null(object@term_env$dag_depth)) {
			depth = object@term_env$dag_depth
			cat("  Max depth:", max(depth), "\n")
		}

		if(n_terms == n_relations + 1) {
			cat("  Aspect ratio: ", object@aspect_ratio[1], ":1\n", sep = "")
		} else {
			cat("  Aspect ratio: ", object@aspect_ratio[1], ":1 (based on the longest distance to root)\n", sep = "")
			cat("                ", object@aspect_ratio[2], ":1 (based on the shortest distance to root)\n", sep = "")
		}

		if(length(object@lt_children_relations)) {
			txt = strwrap(paste(attr(object@lt_children_relations, "levels"), collapse = ", "), width = 60)
			txt[1] = paste0("  Relations: ", txt[1])
			txt[-1] = paste0("             ", txt[-1])
			cat(txt, sep = "\n")

			if(!is.null(object@relations_DAG)) {
				if(length(intersect(object@relations_DAG@terms, attr(object@lt_children_relations, "levels")))) {
					cat("  Relation types may have hierarchical relations.\n")
				}
			}
		}

		if(length(object@annotation$list)) {
			cat("  Annotations are available.\n")
		}

		if(!is.null(object@elementMetadata)) {
			cat("\n")
			cat("With the following columns in the metadata data frame:\n")
			txt = strwrap(paste(colnames(object@elementMetadata), collapse = ", "), width = 60)
			txt = paste0("  ", txt)
			cat(txt, sep = "\n")
		}
	}
)

#' Create the ontology_DAG object from the GO.db package
#' 
#' @param namespace One of "BP", "CC" and "MF".
#' @param relations Types of the GO term relations. In the **GO.db** package, the GO term relations can be "is_a", "part_of",
#'               "regulates", "negatively regulates", "positively regulates". Note since "regulates" is a parent relation
#'               of "negatively regulates", "positively regulates", if "regulates" is selected, "negatively regulates" and "positively regulates"
#'               are also selected. Note "is_a" is always included.
#' @param org_db The name of the organism package or the corresponding database object, e.g. `"org.Hs.eg.db"` or 
#'            directly the [`org.Hs.eg.db::org.Hs.eg.db`] object for human, then the gene annotation to GO terms will be added
#'            to the object.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db()
#' dag = create_ontology_DAG_from_GO_db("BP", 
#'     relations = NULL, org_db = "org.Hs.eg.db")
#' }
create_ontology_DAG_from_GO_db = function(namespace = "BP", relations = "part of", org_db = NULL) {

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

	l = df[, 3] == "isa"
	df[l, 3] = "is_a"
	df[, 3] = gsub(" ", "_", df[, 3])

	if(length(relations) == 0) {
		relations = character(0)
	}

	if(identical(relations, NA)) {
		relations = character(0)
	}

	relations = c("is_a", relations)
	if("regulates" %in% relations) {
		relations = c(relations, "negatively regulates", "positively regulates")
	}
	relations = gsub(" ", "_", relations)
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

	relations_DAG = create_ontology_DAG(c("regulates", "regulates"), c("negatively regulates", "positively regulates"))

	dag = create_ontology_DAG(parents = df[, 2], children = df[, 1], relations = df[, 3], relations_DAG = relations_DAG,
		annotation = annotation, source = paste0("GO ", namespace, " / GO.db package"))

	go = GO.db::GOTERM[dag@terms]
	meta = data.frame(id = AnnotationDbi::GOID(go),
		        name = AnnotationDbi::Term(go),
		        definition = AnnotationDbi::Definition(go))
	rownames(meta) = dag@terms

	mcols(dag) = meta

	dag
}


#' Create sub-DAGs
#' 
#' @param x An `ontology_DAG` object.
#' @param i A single term name. The value should be a character vector. This corresponds to the roots.
#' @param j A single term name. The value should be a character vector. This corresponds to the leaves.
#' @param ... Ignored.
#' @param drop Ignored.
#' 
#' @details It returns a sub-DAG taking node `i` as the root and `j` as the leaves.
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
	signature = c("ontology_DAG"),
	definition = function(x, i, j, ..., drop = FALSE) {

	if(!missing(i)) {
		if(!is.character(i)) {
			stop("Only the character term name should be used as index.")
		}
	}
	if(!missing(j)) {
		if(!is.character(j)) {
			stop("Only the character term name should be used as index.")
		}
	}

	if(!missing(i) && missing(j)) {  ## dag[i] or dag[i, ]
		dag_filter(x, root = i)
	} else if(!missing(i) && !missing(j)) {  ## dag[i, j]
		dag_filter(x, root = i, leaves = j)
	} else if(missing(i) && !missing(j)) {  ## dag[, j]
		dag_filter(x, leaves = j)
	} else {    ## dag[, ]
		x
	}
	
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
#' dag_n_terms(dag)
dag_all_terms = function(dag) {
	dag@terms
}

#' @rdname dag_all_terms
#' @export
dag_n_terms = function(dag) {
	dag@n_terms
}

#' Root or leaves of the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param in_labels Whether the terms are represented in their names or as integer indices?
#' 
#' @return A character or an integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_root(dag)
#' dag_leaves(dag)
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

#' @param terms A vector of term names.
#' 
#' @rdname dag_root
#' @export
dag_is_leaf = function(dag, terms) {
	id = term_to_node_id(dag, terms)
	l = id %in% dag@leaves
	structure(l, names = dag@terms[id])
}

#' Convert to an igraph object
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @details If `relations` is already set in [`create_ontology_DAG()`], relations are also set as an edge attribute in the [`igraph::igraph`] object.
#' 
#' @return An [`igraph::igraph`] object.
#' @export
#' @import igraph
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

	g = graph_from_edgelist(cbind(parents, children))
	V(g)$name = dag@terms

	if(!is.null(relations)) {
		g = set_edge_attr(g, "relation", value = relations)
	}

	g
}

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

	l = rep(TRUE, length(parents))

	if(!is.null(terms)) {
		terms = term_to_node_id(dag, terms)
		l = l & parents %in% terms & children %in% terms
	}

	if(!is.null(relations)) {

		if("is_a" %in% relations) {
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


#' Get meta columns on DAG
#' 
#' @param x An `ontology_DAG` object.
#' @param use.names Please ignore.
#' @param ... Other argument. For `mcols()`, it can be a vector of column names in the meta data frame.
#' @exportMethod mcols
#' @importMethodsFrom S4Vectors mcols
setMethod("mcols", 
	signature = "ontology_DAG",
	definition = function(x, use.names = TRUE, ...) {
	columns = unlist(list(...))
	if(length(columns)) {
		x@elementMetadata[, columns]
	} else {
		x@elementMetadata
	}
})


#' Set meta columns on DAG
#' 
#' @param x An `ontology_DAG` object.
#' @param ... Other argument. For `mcols()`, it can be a vector of column names in the meta data frame.
#' @param value A data frame or a matrix where rows should correspond to terms in `x@terms`.
#' @exportMethod 'mcols<-'
#' @importMethodsFrom S4Vectors 'mcols<-'
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' mcols(dag) = data.frame(id = letters[1:6], v = 1:6)
#' mcols(dag)
#' mcols(dag, "id")
#' dag
setMethod("mcols<-", 
	signature = "ontology_DAG",
	definition = function(x, ..., value) {

	if(is.null(value)) {
		return(x)
	}
	
	df = as.data.frame(value)
	if(nrow(df) != x@n_terms) {
		stop("value should be a table with the same number of rows as the number of terms in DAG.")
	}

	x@elementMetadata = value
	invisible(x)
})

