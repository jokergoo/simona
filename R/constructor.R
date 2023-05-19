

OntologyDAG = setClass("OntologyDAG",
	slots = c("lt_parents" = "list",
		      "lt_children" = "list",
		      "lt_offspring" = "list",
		      "lt_ancestor" = "list",
		      "dist_offspring" = "list",
		      "dist_ancestor" = "list",
		      "n_terms" = "numeric",
		      "n_relations" = "numeric",
		      "name" = "character",
		      "terms" = "character",
		      "root" = "integer",
		      "leaves" = "integer",
		      "stats_env" = "environment",
		      "graph" = "ANY"
		      ),
	prototype = list(
			stats_env = new.env(parent = emptyenv())
		)
)

createOntologyDAG = function(parents, children, relations = NULL, name = "Ontology DAG", verbose = TRUE,
	add_super_root = NULL) {

	if(length(unique(parents)) > length(unique(children))) {
		stop("Found the number of parents is larger than the number of children. Maybe you mis-set the arguments `parents` and `children`?")
	}

	if(length(intersect(parents, children)) == 0) {
		stop("No overlap between parents and children.")
	}

	if(!is.character(parents)) {
		warning("`parents` is better in character mode.")
	}

	if(!is.character(children)) {
		warning("`children` is better in character mode.")
	}

	if(verbose) cat("* creating the DAG as an igraph object...\n")
	g = graph.edgelist(cbind(parents, children), directed = TRUE)
	if(!is.null(relations)) {
		relations = gsub("^is\\W*a$", "isa", relations, ignore.case = TRUE)
		relations = gsub("^part\\W*of$", "part of", relations, ignore.case = TRUE)
		E(g)$relation = relations
	}

	root = unname(which(igraph::degree(g, mode = "in") == 0))
	leaves = unname(which(igraph::degree(g, mode = "out") == 0))
	terms = V(g)$name
	n_terms = vcount(g)
	if(is.null(terms)) {
		terms = as.character(seq_len(vcount(g)))
	}
	g = delete_vertex_attr(g, "name")

	if(!is_dag(g)) {		
		stop("The ontology is not a DAG.")
	}

	if(length(root) > 1) {
		if(identical(add_super_root, TRUE)) {
			super_root = n_terms + 1L
			g = add_vertices(g, 1)
			for(r in root) {
				g = add_edges(g, c(super_root, r))
			}
			terms = c(terms, "_super_root_")
			root = super_root
		} else if(identical(add_super_root, NULL)) {
			opt = readline(prompt = "There are more than one root. Do you want to add a super root? [y|n] ")
			if(opt %in% c("y", "")) {
				super_root = n_terms + 1L
				g = add_vertices(g, 1)
				for(r in root) {
					g = add_edges(g, c(super_root, r))
				}
				terms = c(terms, "_super_root_")
				root = super_root
			} else {
				warning("There are more than one root.")
			}
		} else {
			warning("There are more than one root.")
		}
	}

	if(verbose) cat("* extracting parents/children/ancestors/offspring relations...\n")
	lt_parents = ego(g, order = 1, mode = "in", mindist = 1)
	lt_children = ego(g, order = 1, mode = "out", mindist = 1)
	lt_ancestor = ego(g, order = 1e8, mode = "in", mindist = 1)
	lt_offspring = ego(g, order = 1e8, mode = "out", mindist = 1)

	lt_parents = lapply(lt_parents, as.vector)
	lt_children = lapply(lt_children, as.vector)
	lt_ancestor = lapply(lt_ancestor, as.vector)
	lt_offspring = lapply(lt_offspring, as.vector)

	# the node itself is always the first element
	for(i in seq_along(lt_ancestor)) {
		lt_ancestor[[i]] = c(i, lt_ancestor[[i]])
	}
	for(i in seq_along(lt_offspring)) {
		lt_offspring[[i]] = c(i, lt_offspring[[i]])
	}

	if(verbose) cat("* calculating distance to ancestor/offspring terms....\n")
	d = distances(g, V(g), V(g), mode = "out", weights = NA)
	dist_ancestor = lapply(seq_len(n_terms), function(i) {
		ind = lt_ancestor[[i]]
		d[ind, i]
	})
	dist_offspring = lapply(seq_len(n_terms), function(i) {
		ind = lt_offspring[[i]]
		d[i, ind]
	})
	rm(d); gc(verbose = FALSE)
	
	if(verbose) cat("* creating the OntologyDAG object...\n")
	dag = OntologyDAG(
		lt_parents = lt_parents,
		lt_children = lt_children,
		lt_offspring = lt_offspring,
		lt_ancestor = lt_ancestor,
		dist_offspring = dist_offspring,
		dist_ancestor = dist_ancestor,
		name = name,
		terms = terms,
		root = root,
		leaves = leaves,
		stats_env = new.env(parent = emptyenv()),
		graph = g
	)

	dag@n_terms = vcount(g)
	dag@n_relations = ecount(g)

	dag@stats_env$n_parents = sapply(lt_parents, length)
	dag@stats_env$n_children = sapply(lt_children, length)
	dag@stats_env$n_offspring = sapply(lt_offspring, length)
	dag@stats_env$n_ancestor = sapply(lt_ancestor, length)

	return(dag)
}

setMethod("show",
	signature = "OntologyDAG",
	definition = function(object) {

		cat("An OntologyDAG object:\n")
		cat("  Name:", object@name, "\n")
		cat("  #Terms:", object@n_terms, "\n")
		cat("  #Relations:", object@n_relations, "\n")
		cat("  Root:", paste(object@terms[object@root], collapse = ", "), "\n")

		if(!is.null(dag@stats_env$dag_depth)) {
			depth = dag@stats_env$dag_depth
			cat("  Max depth:", max(depth), "\n")
			cat("  Width:", max(table(depth)), "\n")
		}
	}
)

createSubDAG = function(dag, term) {
	g = dag@graph
	V(g)$name = dag@terms

	if(is.character(term)) {
		term = which(dag@terms == term)
	}

	lt_offspring = dag@lt_offspring
	offspring = lt_offspring[[term]]

	if(length(offspring) == 0) {
		return(NULL)
	}

	subg = induced_subgraph(g, offspring, impl = "copy_and_delete")

	el = get.edgelist(subg)
	parents = el[, 1]
	children = el[, 2]
	relations = E(subg)$relation

	createOntologyDAG(parents, children, relations)
}

