

polar2Cartesian = function(d) {
    theta = d[, 1]/180*pi
    rou = d[, 2]
    x = rou * cos(theta)
    y = rou * sin(theta)
    return(cbind(x, y))
}


#' Visualize the DAG
#' 
#' @param dag An `ontology_Dag` object.
#' 
#' @details
#' `dag_circular_viz()` uses a circular layout for visualizing large DAGs. `dag_graphviz()`
#' uses a hierarchical layout for visualizing small DAGs.
#' 
#' @import grid
#' @importFrom Polychrome alphabet.colors
#' @rdname dag_viz
#' @export
dag_circular_viz = function(dag) {

	term_pos = cpp_term_pos_on_circle(dag) ## in polar coordinate

	level1 = setdiff(unique(term_pos$level1_group), 0)
	all_colors = alphabet.colors(length(level1))
	names(all_colors) = level1
	all_colors["0"] = "black"
	col = rep("black", nrow(term_pos))
	col = all_colors[as.character(term_pos$level1_group)]
	
	lt_children = dag@lt_children
	n = dag@n_terms
	df_edge = data.frame(from = rep(seq_len(n), times = sapply(lt_children, length)),
		                 to = unlist(lt_children))

	term_pos2 = polar2Cartesian(term_pos) ## xy coordinate

	x1 = term_pos2[df_edge$from, 1]
	y1 = term_pos2[df_edge$from, 2]
	x2 = term_pos2[df_edge$to, 1]
	y2 = term_pos2[df_edge$to, 2]

	max_depth = max(abs(term_pos$rho))

	n_children = dag@term_env$n_children

	grid.newpage()
	pushViewport(viewport(x = unit(0, "npc"), just = "left", width = unit(1, "snpc"), height = unit(1, "snpc"), 
		xscale = c(-max_depth, max_depth), yscale = c(-max_depth, max_depth)))
	grid.points(term_pos2[, 1], term_pos2[, 2], default.units = "native", pch = 16, 
		gp = gpar(col = add_transparency(col, 0.5)), size = unit(4, "pt"))

	grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(col = "#00000010"))

}

#' @importFrom circlize rand_color
dag_as_DOT = function(dag) {

	name = NULL
	if(!is.null(mcols(dag))) {
		meta = mcols(dag)
		if("name" %in% colnames(meta)) {
			name = meta[, "name"]
		}
	}

	if(is.null(name)) {
		nodes = paste(
			"node [shape = box, fontname = Helvetical]",
			paste( paste0("\"", dag@terms, "\""), collapse = ";\n"),
			sep = "\n"
		)
	} else {
		nodes = paste(
			"node [shape = box, fontname = Helvetical]",
			paste( paste0("\"", dag@terms, "\"", " [tooltip = \"", name, "\"]"), collapse = ";\n"),
			sep = "\n"
		)
	}

	lt_children = dag@lt_children
	children = unlist(lt_children)
	parents = rep(seq_along(dag@terms), times = sapply(lt_children, length))
	children = dag@terms[children]
	parents = dag@terms[parents]

	if(length(dag@lt_children_relations) > 0) {
		relation_levels = attr(dag@lt_children_relations, "levels")
		v_relations = unlist(dag@lt_children_relations)
		v_relations = relation_levels[v_relations]
		n_levels = length(relation_levels)

		if(n_levels <= 8) {
			relation_col = structure(seq_len(n_levels), names = relation_levels)
		} else if(n_levels <= 26) {
			relation_col = structure(alphabet.colors(n_levels), names = relation_levels)	
		} else {
			relation_col = structure(rand_color(n_levels), names = relation_levels)
		}
	} else {
		v_relations = NULL
	}

	if(is.null(v_relations)) {
		edges = paste(
			"edge [dir=\"back\"]",
			paste( paste0("\"", parents, "\""), "->", paste0("\"", children, "\""), sep = "", collapse = ";\n"),
			sep = "\n"
		)
	} else {
		edges = paste(
			"edge [dir=\"back\"]",
			paste( paste0("\"", parents, "\""), "->", paste0("\"", children, "\""), " [tooltip=\"", v_relations, "\" color=\"", relation_col[v_relations], "\"]", sep = "", collapse = ";\n"),
			sep = "\n"
		)
	}

	DOT = paste(
		"digraph {",
		"graph [overlap = true, fontsize = 10]",
		nodes,
		"\n",
		edges,
		"}",
		"\n",
		sep = "\n"
	)

	DOT
}

#' @rdname dag_viz
#' @export
dag_graphviz = function(dag) {
	if(dag@n_terms > 100) {
		warning("Gviz is only effcient for visualizing small graphs.")
	}
	if(dag@n_terms > 500) {
		stop("Too many terms.")
	}
	check_pkg("DiagrammeR", bioc = FALSE)
	DiagrammeR::grViz(dag_as_DOT(dag))
}
