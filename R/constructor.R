


#' The ontology_DAG class
#' 
#' This class defines the DAG structure of an ontology.
#' 
#' @slot terms A character vector of length `n` of all term names. Other slots that store term-level information use the integer indices of terms.
#' @slot n_terms An integer scalar of the total number of terms in the DAG.
#' @slot n_relations An integer scalar of the total number of relations in the DAG.
#' @slot lt_parents A list of length `n`. Each element in the list is an integer index vector of the parent terms of the i^th term.
#' @slot lt_children A list of length `n`. Each element in the list is an integer index vector of the child terms of the i^th term.
#' @slot lt_children_relations A list of length `n`. Each element is a vector of the semantic relations between the i^th term and its child terms, e.g. a child "is_a" parent.
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
#' @returns An `ontology_DAG` object.
#' @examples
#' 1
#' # This function should not be used directly.
ontology_DAG = setClass("ontology_DAG",
	slots = c("terms" = "character",
		      "n_terms" = "integer",
		      "n_relations" = "integer",
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
#' @param parents A character vector of parent terms. You can also construct the `ontology_DAG` object by a list of parent-child links. See **Examples**.
#' @param children A character vector of child terms. 
#' @param relations A character vector of parent-child relations, e.g. "is_a", "part_of", or self-defined semantic relations.
#'        If it is set, it should have the same length as `parents` and `children`.
#' @param relations_DAG If the relation types have hierarchical relations, it can also be constructed by `create_ontology_DAG()` first. See **Examples**.
#'          When the DAG for relation types is provided, the ancestor/offspring relationship of relation types will be taken into consideration automatically.
#' @param source Source of the ontology. It is only used as a label of the object.
#' @param annotation A list of character vectors which contain items annotated to the terms. Names of the list should be the term names. In the DAG, items
#'                   annotated to a term will also be annotated to its parents. Such merging
#'                   is applied automatically in the package.
#' @param remove_cyclic_paths Whether to remove cyclic paths If a cyclic path is represented as `[a, b, ..., z, a]`,
#'       the last link (i.e. `z->a`) is simply removed. If the value is set to `FALSE` and if there are cyclic paths, there
#'       will be an error that lists all cyclic paths.
#' @param remove_rings There might be rings that are isolated to the main DAG where there are no roots on the rings, thus they cannot be attached to the main DAG. If the value
#'        of `remove_rings` is set to `TRUE`, such rings are removed.
#' @param verbose Whether to print messages.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom methods new
#' @importFrom stats quantile
#' @import Rcpp
#' @useDynLib simona, .registration = TRUE
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
#' 
#' # with a list of parent-child relations
#' dag = create_ontology_DAG(c("a-b", "a-c", "b-c", "b-d", "c-e", "e-f"))
create_ontology_DAG = function(parents, children, relations = NULL, relations_DAG = NULL,
	source = "Ontology", annotation = NULL, remove_cyclic_paths = FALSE, remove_rings = FALSE,
	verbose = simona_opt$verbose) {

	if(missing(children)) {
		if(any(grepl("\\s+-\\s+", parents))) {
			lt = strsplit(parents, "\\s+-\\s+")
		} else {
			lt = strsplit(parents, "\\s*-\\s*")
		}
		
		parents = sapply(lt, function(x) x[1])
		children = sapply(lt, function(x) x[2])
	}

	if(!is.character(parents)) {
		stop("`parents` should be in character mode.")
	}

	if(!is.character(children)) {
		stop("`children` should be in character mode.")
	}

	l = parents == children
	if(any(l)) {
		if(verbose) {
			message("removed ", sum(l), " relations where parents are the same as children.")
		}
		parents = parents[!l]
		children = children[!l]
	}

	if(length(parents) == 0 || length(children) == 0) {
		stop("There is no relation.")
	}

	has_relations = !is.null(relations)
	lt_relations = list()

	if(has_relations && any(l)) {
		relations = relations[!l]
	}

	terms = sort(unique(c(parents, children)))
	n_terms = length(terms)
	terms2ind = structure(seq_along(terms), names = terms)

	lt_parents2 = split(unname(terms2ind[parents]), children)
	lt_children2 = split(unname(terms2ind[children]), parents)
	if(has_relations) {
		relations = normalize_relation_type(relations)
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

	root = which(vapply(lt_parents, length, FUN.VALUE = integer(1)) == 0)
	leaves = which(vapply(lt_children, length, FUN.VALUE = integer(1)) == 0)

	if(length(root) == 0) {
		stop("Cannot find the root. There might exist a cycle.")
	}
	
	if(length(root) > 1) {
		if(length(root) > 5) {
			txt = strwrap(paste(terms[root][seq_len(5)], collapse = ", "), width = 80)
			txt[length(txt)] = paste0(txt[length(txt)], ",")
			txt = c(txt, qq("  and other @{length(root)-5} term@{ifelse(length(root)-5 == 1, '', 's')} ..."))
		} else {
			txt = strwrap(paste(terms[root], collapse = ", "), width = 80)
		}
		txt = paste(paste0("  ", txt), collapse = "\n")
		if(verbose) {
			message("There are more than one root:\n", txt, "\n  A super root (", SUPER_ROOT, ") is added.")
		}
		
		super_root = n_terms + 1L
		lt_parents[[super_root]] = integer(0)
		for(r in root) {
			lt_parents[[r]] = super_root
		}

		lt_children[[super_root]] = root
		if(has_relations) lt_relations[[super_root]] = rep(which(relation_levels == "is_a"), length(root))
		terms = c(terms, SUPER_ROOT)
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
		n_relations = sum(vapply(lt_children, length, FUN.VALUE = integer(1))),
		lt_parents = lt_parents,
		lt_children = lt_children,
		lt_children_relations = lt_relations, 
		relations_DAG = relations_DAG,
		source = source,
		root = root,
		leaves = leaves,
		term_env = new.env(parent = emptyenv())
	)

	cyclic_paths = cpp_check_cyclic_node(dag)

	if(length(cyclic_paths)) {
		if(remove_cyclic_paths) {
			message_wrap(qq("Removed @{length(cyclic_paths)} cyclic path@{ifelse(length(cyclic_paths) == 1, '', 's')}."))
			removed_p = vapply(cyclic_paths, function(x) x[length(x)-1], FUN.VALUE = integer(1))
			removed_c = vapply(cyclic_paths, function(x) x[length(x)], FUN.VALUE = integer(1))

			for(ir in seq_along(cyclic_paths)) {
				ip = removed_p[ir]
				ic = removed_c[ir]

				lt_parents[[ic]] = setdiff(lt_parents[[ic]], ip)
				ind = lt_children[[ip]] != ic
				lt_children[[ip]] = lt_children[[ip]][ind]
				if(length(lt_relations)) {
					lt_relations[[ip]] = lt_relations[[ip]][ind]
				}
			}
			root = which(vapply(lt_parents, length, FUN.VALUE = integer(1)) == 0)
			leaves = which(vapply(lt_children, length, FUN.VALUE = integer(1)) == 0)

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
		} else {
			if(verbose) {
				print_cyclic_paths(cyclic_paths, terms)
			}
			stop_wrap("Found cyclic nodes (paths are listed above). Set `remove_cyclic_paths = TRUE` to remove them.")
		}
	}

	dag@term_env$n_parents = vapply(lt_parents, length, FUN.VALUE = integer(1))
	dag@term_env$n_children = vapply(lt_children, length, FUN.VALUE = integer(1))

	tpl_order = order(dag_depth(dag), dag@term_env$n_children, dag@term_env$n_parents, dag@terms)
	dag@tpl_sorted = seq_len(n_terms)[tpl_order]
	dag@tpl_pos = seq_len(n_terms)
	
	# pos has the same order as terms, but the value are the positions on `sorted`, i.e. the rank
	dag@tpl_pos[tpl_order] = seq_len(n_terms)

	dag@annotation = list(list = vector("list", 0), names = character(0))
	dag = add_annotation(dag, annotation)

	depth = dag_depth(dag)
	tb1 = table(depth)
	aspect_ratio1 = max(tb1)/quantile(depth, 0.99)
	aspect_ratio1 = round(aspect_ratio1, 2)

	dist = cpp_dag_dist_from_root(dag)
	tb2 = table(dist)
	aspect_ratio2 = max(tb2)/quantile(dist, 0.99)
	aspect_ratio2 = round(aspect_ratio2, 2)

	dag@aspect_ratio = c(aspect_ratio1, aspect_ratio2)

	iring = which(depth < 0)
	if(length(iring)) {
		if(remove_rings) {
			if(verbose) {
				message_wrap(qq("Removed @{length(iring)} terms in isolated rings."))
			}
			return(dag[[dag_root(dag)]])
		} else {
			for(i in iring) {
				cyclic_paths = cpp_check_cyclic_node(dag, i)
				print_cyclic_paths(cyclic_paths[1], terms)
				stop_wrap("Found isolated rings (one path is listed above). Set `remove_rings = TRUE` to remove them.")
			}	
		}
	}

	return(dag)
}

add_annotation = function(dag, annotation) {
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
	dag
}

print_cyclic_paths = function(cyclic_paths, terms) {
	for(i in seq_len(length(cyclic_paths))) {
		message("[", appendLF = FALSE)
		msg = paste(terms[cyclic_paths[[i]]], collapse = " ~ ")
		message(msg, appendLF = FALSE)
		message("]", appendLF = TRUE)
	}
}

#' Print the ontology_DAG object
#' 
#' @param object An `ontology_DAG` object.
#' @returns No value is returned.
#' @exportMethod show
setMethod("show",
	signature = "ontology_DAG",
	definition = function(object) {

		n_terms = object@n_terms
		n_relations = object@n_relations

		cat("An ontology_DAG object:\n")
		cat("  Source:", object@source, "\n")
		if(n_terms == n_relations + 1) {
			cat("  ", n_terms, " terms / ", n_relations, " relations / a tree\n", sep = "")
		} else {
			cat("  ", n_terms, " terms / ", n_relations, " relations\n", sep = "")
		}

		cat("  Root:", paste(object@terms[object@root], collapse = ", "), "\n")
		txt = strwrap(paste(c(object@terms[seq_len(min(4, object@n_terms))], "..."), collapse = ", "), width = 60)
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
			cat("  Avg number of parents: ", sprintf("%.2f", mean(n_parents(object)[-dag_root(object, in_labels = FALSE)])), "\n", sep = "")
			cat("  Avg number of children: ", sprintf("%.2f", mean(n_parents(object)[-dag_leaves(object, in_labels = FALSE)])), "\n", sep = "")
			cat("  Aspect ratio: ", object@aspect_ratio[1], ":1 (based on the longest distance from root)\n", sep = "")
			cat("                ", object@aspect_ratio[2], ":1 (based on the shortest distance from root)\n", sep = "")
		}

		if(length(object@lt_children_relations)) {
			rel_levels = attr(object@lt_children_relations, "levels")
			txt = strwrap(paste(rel_levels, collapse = ", "), width = 60)
			txt[1] = paste0( "  Relations: ", txt[1])
			txt[-1] = paste0("             ", txt[-1])
			cat(txt, sep = "\n")

			if(!is.null(object@relations_DAG)) {
				if(length(intersect(object@relations_DAG@terms, attr(object@lt_children_relations, "levels")))) {
					cat("  Relation types may have hierarchical relations.\n")
				}
			}
		}

		if(length(object@annotation$list)) {
			cat("  Annotations:", length(object@annotation$names), "items\n")
			txt = strwrap(paste(c(object@annotation$names[seq_len(min(4, length(object@annotation$names)))], "..."), collapse = ", "), width = 60)
			txt = paste0("               ", txt)
			cat(txt, sep = "\n")
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
#'            to the object. For other non-model organisms, consider to use the **AnnotationHub** package to find one.
#' @param evidence_code A vector of evidence codes for gene annotation to GO terms. See \url{https://geneontology.org/docs/guide-go-evidence-codes/}.
#' @param verbose Whether to print messages.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom utils getFromNamespace
#' @examples
#' dag = create_ontology_DAG_from_GO_db()
#' dag
create_ontology_DAG_from_GO_db = function(namespace = "BP", relations = "part of", org_db = NULL, 
	evidence_code = NULL, verbose = simona_opt$verbose) {

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
	relations = normalize_relation_type(relations)
	if("regulates" %in% relations) {
		relations = c(relations, "negatively_regulates", "positively_regulates")
	}
	df = df[df[, 3] %in% relations, , drop = FALSE]

	if(verbose) {
		message("relations: ", paste(relations, collapse = ", "))
	}

	if(!is.null(org_db)) {
		if(is.character(org_db)) {
			check_pkg(org_db, bioc = TRUE)
			org_db = getFromNamespace(org_db, ns = org_db)
		}
	
		suppressMessages(tb <- AnnotationDbi::select(org_db, keys = AnnotationDbi::keys(org_db), columns = c("GO", "EVIDENCE", "ONTOLOGY")))
		tb = tb[tb$ONTOLOGY == namespace, , drop = FALSE]
		if(!is.null(evidence_code)) {
			tb = tb[tb$EVIDENCE %in% evidence_code, , drop = FALSE]
			if(row(tb) == 0) {
				stop("No annotation left after filtering by the evidence codes.")
			}
		}
		tb = tb[, c(1, which(colnames(tb) == "GO")), drop = FALSE]
		annotation = split(tb[, 1], tb[, 2])
		annotation = lapply(annotation, unique)
		annotation = lapply(annotation, as.character)
	} else {
		annotation = NULL
	}

	relations_DAG = create_ontology_DAG(c("regulates", "regulates"), c("negatively regulates", "positively regulates"))

	go_db_version = read.dcf(system.file("DESCRIPTION", package = "GO.db"))[1, "Version"]
	dag = create_ontology_DAG(parents = df[, 2], children = df[, 1], relations = df[, 3], relations_DAG = relations_DAG,
		annotation = annotation, source = paste0("GO ", namespace, " / GO.db package ", go_db_version))

	go = GO.db::GOTERM[dag@terms]
	meta = data.frame(id = AnnotationDbi::GOID(go),
		        name = AnnotationDbi::Term(go),
		        definition = AnnotationDbi::Definition(go))
	rownames(meta) = dag@terms

	mcols(dag) = meta

	dag
}

#' Names of all terms
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @return `dag_all_terms()` returns a vector of term names. `dag_n_terms()` returns
#'      a single iteger.
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

#' @rdname dag_all_terms
#' @export
dag_n_relations = function(dag) {
	sum(sapply(dag@lt_children, length))
}

#' @rdname dag_all_terms
#' @export
dag_n_leaves = function(dag) {
	length(dag@leaves)
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
	parents = rep(seq_along(dag@terms), times = vapply(lt_children, length, FUN.VALUE = integer(1)))
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


#' Get or set meta columns on DAG
#' 
#' @param x An `ontology_DAG` object.
#' @param use.names Please ignore.
#' @param ... Other argument. For `mcols()`, it can be a vector of column names in the meta data frame.
#' @exportMethod mcols
#' @importMethodsFrom S4Vectors mcols
#' @rdname mcols_ontology_DAG
#' @returns A data frame.
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


#' @param x An `ontology_DAG` object.
#' @param ... Other argument. For `mcols()`, it can be a vector of column names in the meta data frame.
#' @param value A data frame or a matrix where rows should correspond to terms in `x@terms`.
#' @exportMethod 'mcols<-'
#' @importMethodsFrom S4Vectors 'mcols<-'
#' @rdname mcols_ontology_DAG
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
		if(dag_root(x) == SUPER_ROOT) {
			df2 = df[seq_len(nrow(df)+1), , drop = FALSE]
			rownames(df2)[nrow(df2)] = SUPER_ROOT
			
			if(! (nrow(df) == x@n_terms - 1 || nrow(df) == x@n_terms) ) {
				stop("value should be a table with the same number of rows as the number of terms in DAG.")
			}

			df = df2
		} else {
			if(nrow(df) != x@n_terms) {
				stop("value should be a table with the same number of rows as the number of terms in DAG.")
			}
		}
	}

	x@elementMetadata = df
	invisible(x)
})


dag_namespaces = function(dag) {
	df = mcols(dag)
	if(is.null(df)) {
		return(NULL)
	}

	if(!"namespace" %in% colnames(df)) {
		return(NULL)
	}

	unique(df$namespace)
}

#' Create the ontology_DAG object from the igraph object
#' 
#' @param g An [`igraph::igraph`] object.
#' @param relations A vector of relation types. The length of the vector should be the same as the number of edges in `g`.
#' @param verbose Whether to print messages.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @import igraph
create_ontology_DAG_from_igraph = function(g, relations = NULL, verbose = simona_opt$verbose) {
	edges = get.edgelist(g)

	create_ontology_DAG(as.character(edges[, 1]), as.character(edges[, 2]), relations = relations, 
		source = "igraph object", verbose = verbose)
}


#' Whether the terms exist in the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term IDs.
#' 
#' @return A logical vector.
#' @export
dag_has_terms = function(dag, terms) {
	terms %in% dag@terms
}
