

polar2Cartesian = function(d) {
    theta = d[, 1]/180*pi
    rou = d[, 2]
    x = rou * cos(theta)
    y = rou * sin(theta)
    return(cbind(x, y))
}

calc_n_neighbours_on_circle = function(theta, width = 1) {
	od = order(order(theta))
	theta2 = sort(theta)

	k = cpp_calc_n_neighbours_on_circle(theta2, width)
	k[od]
}

#' Visualize the DAG
#' 
#' @param dag An `ontology_Dag` object.
#' @param highlight A vector of terms to be highlighted on the DAG.
#' @param start Start of the circle, measured in degree.
#' @param end End of the circle, measured in degree.
#' @param reorder_level Whether to reorder child terms. See [`dag_reorder()`].
#' @param partition_level If `node_col` is not set, users can cut the DAG into clusters with different node colors.
#' @param node_col Colors of nodes. If the value is a vector, the order should correspond to terms in [`dag_all_terms()`].
#' @param node_transparency Transparency of nodes. The same format as `node_col`.
#' @param node_size Size of nodes. The same format as `node_col`.
#' @param edge_col A named vector where names correspond to relation types.
#' @param edge_transparency A named vector where names correspond to relation types.
#' @param legend_labels_from If partitioning is applied on the DAG, a legend is generated showing different top
#'         terms. By default, the legend labels are the term IDs. If there are additionally column stored
#'         in the meta data frame of the DAG object, the column name can be set here to replace the term IDs as
#'         legend labels.
#' @param legend_labels_max_width Maximal width of legend labels. Labels are wrapped into
#'        multiple lines if the widths exceed it.
#' 
#' @details
#' `dag_circular_viz()` uses a circular layout for visualizing large DAGs. `dag_graphviz()`
#' uses a hierarchical layout for visualizing small DAGs.
#' 
#' @import grid
#' @importFrom Polychrome alphabet.colors
#' @import ComplexHeatmap
#' @rdname dag_viz
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db()
#' dag_circular_viz(dag)
#' }
dag_circular_viz = function(dag, highlight = NULL, start = 0, end = 360,
	reorder_level = NA, partition_level = 1,
	node_col = NULL, node_transparency = 0.5, node_size = NULL, 
	edge_col = NULL, edge_transparency = 0.98,
	legend_labels_from = NULL, legend_labels_max_width = 50) {

	if(end < start) {
		stop_wrap("`end` should be larger than `start`.")
	}

	n_terms = dag@n_terms

	if(!is.na(reorder_level)) {
		message(qq("reordering children on depth <= @{reorder_level} in DAG..."))
		dag = dag_reorder(dag, max_level = reorder_level)
	}

	has_highlight = FALSE
	if(!is.null(highlight)) {
		has_highlight = TRUE
		l_highlight = dag@terms %in% highlight
	}

	if(dag_is_tree(dag)) {
		tree = dag
	} else {
		message("converting DAG to a tree...")
		tree = dag_treelize(dag)
	}
	message("calculating term positions on the DAG...")
	term_pos = cpp_term_pos_on_circle(tree, n_offspring(tree), start, end) ## in polar coordinate

	term_pos$n_neighbours = 0
	for(level in sort(unique(term_pos$rho))) {
		l = term_pos$rho == level
		message(qq("calculating numbers of neighbours within 1 degree neighbourhood on level @{level}, @{sum(l)} terms..."))
		term_pos[l, "n_neighbours"] = calc_n_neighbours_on_circle(term_pos$theta[l], width = 0.5)
	}

	node_col_map = NULL
	if(is.null(node_col)) {
		group = partition_by_level(dag, level = partition_level, term_pos = term_pos)
		level1 = unique(group); level1 = level1[!is.na(level1)]
		n_levels = length(level1)
		
		ind = which( dag@terms %in% level1)
		level1 = dag@terms[ind][order(term_pos[ind, "theta"])]

		if(n_levels <= 7) {
			default_col = c("#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "#9E9E9E")
			node_col_map = default_col[seq_len(n_levels)]
		} else if(n_levels <= 26) {
			node_col_map = alphabet.colors(n_levels)
		} else {
			node_col_map = rand_color(n_levels)
		}

		names(node_col_map) = level1
		node_col = node_col_map[as.character(group)]
		node_col[is.na(node_col)] = "black"
	}
	if(is.null(node_size)) {
		n_children = dag@term_env$n_children
		node_size = .scale(n_children, c(0, quantile(n_children, 0.99)+1), c(2, 10))
	}

	if(!length(node_col) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `node_col` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(node_transparency) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `node_transparency` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(node_size) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `node_size` should be one or @{n_terms} (number of terms in the DAG)."))
	}

	node_transparency = .scale(term_pos$n_neighbours, c(1, quantile(term_pos$n_neighbours, 0.9)+1), c(0.9, 0.5))

	node_col = add_transparency(node_col, node_transparency)

	if(length(dag@lt_children_relations) > 0) {
		relation_levels = attr(dag@lt_children_relations, "levels")
		v_relations = unlist(dag@lt_children_relations)
		v_relations = relation_levels[v_relations]
		n_levels = length(relation_levels)

		if(is.null(edge_col)) {
			if(n_levels <= 8) {
				default_col = c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "#9E9E9E")
				edge_col = structure(default_col[seq_len(n_levels)], names = relation_levels)
			} else if(n_levels <= 26) {
				edge_col = structure(alphabet.colors(n_levels), names = relation_levels)	
			} else {
				edge_col = structure(rand_color(n_levels), names = relation_levels)
			}
		} else if(is.atomic(edge_col)) {
			if(length(edge_col) == 1) {
				edge_col = structure(rep(edge_col, n_levels), names = relation_levels)
			}
		}

		if(length(setdiff(names(edge_col), relation_levels))) {
			stop("names in `edge_col` should cover all relation types.")
		}

		if(is.atomic(edge_transparency)) {
			if(length(edge_transparency) == 1) {
				edge_transparency = structure(rep(edge_transparency, n_levels), names = relation_levels)
			}
		}

		if(length(setdiff(names(edge_transparency), relation_levels))) {
			stop("names in `edge_transparency` should cover all relation types.")
		}

		edge_col_v = add_transparency(edge_col[v_relations], edge_transparency)
	} else {
		if(is.null(edge_col)) {
			edge_col_v = "black"
		}
		edge_col_v = add_transparency(edge_col_v, edge_transparency)
	}

	if(has_highlight) {
		node_col[!l_highlight] = add_transparency("black", 0.95)
		node_col[l_highlight] = add_transparency(node_col[l_highlight], 0)
		node_size[!l_highlight] = 2
	}
	
	lt_children = dag@lt_children
	n = dag@n_terms
	df_edge = data.frame(from = rep(seq_len(n), times = sapply(lt_children, length)),
		                 to = unlist(lt_children))
	# adjust rho
	tb = log(table(term_pos$rho))
	rho2 = cumsum(tb - tb["0"])/sum(tb - tb["0"])
	term_pos$rho = rho2[as.character(term_pos$rho)]

	term_pos2 = polar2Cartesian(term_pos) ## xy coordinate

	x1 = term_pos2[df_edge$from, 1]
	y1 = term_pos2[df_edge$from, 2]
	x2 = term_pos2[df_edge$to, 1]
	y2 = term_pos2[df_edge$to, 2]

	max_depth = max(abs(term_pos$rho))

	message("making plot...")
	grid.newpage()
	pushViewport(viewport(x = unit(0, "npc"), just = "left", width = unit(1, "snpc"), height = unit(1, "snpc"), 
		xscale = c(-max_depth, max_depth), yscale = c(-max_depth, max_depth)))
	if(!is.null(node_col_map)) {
		for(nm in names(node_col_map)) {
			i = which(dag@terms == nm)
			draw_sector(term_pos[i, "theta"] - term_pos[i, "width"]/2,
				        term_pos[i, "theta"] + term_pos[i, "width"]/2,
				        max_depth, gp = gpar(fill = add_transparency(node_col_map[[nm]], 0.9, FALSE), col = NA))
		}
	}

	if(has_highlight) {
		grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(col = edge_col_v))
		grid.points(term_pos2[!l_highlight, 1], term_pos2[!l_highlight, 2], default.units = "native", pch = 16, 
			gp = gpar(col = node_col[!l_highlight]), size = unit(node_size[!l_highlight], "pt"))
		grid.points(term_pos2[l_highlight, 1], term_pos2[l_highlight, 2], default.units = "native", pch = 16, 
			gp = gpar(col = node_col[l_highlight]), size = unit(node_size[l_highlight], "pt"))
	} else {
		grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(col = edge_col_v))
		grid.points(term_pos2[, 1], term_pos2[, 2], default.units = "native", pch = 16, 
			gp = gpar(col = node_col), size = unit(node_size, "pt"))
	}

	lgd_list = list()
	if(!is.null(node_col_map)) {

		if(has_highlight) {
			node_col_map = node_col_map[intersect(names(node_col_map), group[l_highlight])]
		} else {
			sector_width = term_pos[sapply(names(node_col_map), function(x) which(dag@terms == x)), "width"]
			node_col_map = node_col_map[sector_width > 1]

		}
		if(length(node_col_map) > 20) {
			n_offspring = n_offspring(dag)
			node_col_map = node_col_map[ order(-n_offspring[names(node_col_map)])[1:20] ]
		}
		if(is.null(legend_labels_from)) {
			legend_labels = names(node_col_map)
		} else {
			legend_labels = mcols(dag)[names(node_col_map), legend_labels_from]
		}
		legend_labels = sapply(legend_labels, function(x) {
			x = strwrap(x, width = legend_labels_max_width)
			if(length(x) > 1) {
				x[-1] = paste0("  ", x[-1])
			}
			paste(x, collapse = "\n")
		})
		lgd_list = c(lgd_list, list(Legend(title = "Top terms", labels = legend_labels,
			type = "points", pch = 16,
			background = add_transparency(node_col_map, 0.9, FALSE), 
			legend_gp = gpar(col = node_col_map))))
	}
	if(!is.null(edge_col)) {
		legend_labels = names(edge_col)
		legend_labels = sapply(legend_labels, function(x) {
			x = strwrap(x, width = legend_labels_max_width)
			if(length(x) > 1) {
				x[-1] = paste0("  ", x[-1])
			}
			paste(x, collapse = "\n")
		})
		lgd_list = c(lgd_list, list(Legend(title = "Relations", labels = legend_labels, type = "lines", legend_gp = gpar(col = edge_col))))
	}

	if(length(lgd_list)) {
		lgd = packLegend(list = lgd_list)
		draw(lgd, x = unit(1, "snpc") + unit(2, "mm"), just = "left")
	}

}

.scale = function(v, range, map) {
	v[v > range[2]] = range[2]
	v[v < range[1]] = range[1]

	(map[2] - map[1])/(range[2] - range[1])*(v - range[1]) + map[1]
}

#' @param color A vector of colors. If the value is a vector, the order should correspond to terms in [`dag_all_terms()`].
#' @param style Style of the nodes. See https://graphviz.org/docs/attr-types/style/ for possible values.
#' @param fontcolor Color of labels.
#' @param fontsize Font size of labels.
#' @param shape Shape of nodes. See https://graphviz.org/doc/info/shapes.html#polygon for possible values.
#' @param edge_color A named vector where names correspond to relation types.
#' @param edge_style A named vector where names correspond to relation types. See https://graphviz.org/docs/attr-types/style/ for possible values for edges.
#'
#' @seealso http://magjac.com/graphviz-visual-editor/ is nice place to try the DOT code.
#' @details `dag_as_DOT()` generates the DOT code of the DAG.
#' @importFrom circlize rand_color
#' @export
#' @rdname dag_viz
dag_as_DOT = function(dag, color = "black", style = "solid",
	fontcolor = "black", fontsize = 10, shape = "box",
	edge_color = NULL, edge_style = NULL) {

	n_terms = dag@n_terms

	if(!length(color) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `color` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(style) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `style` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(fontcolor) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `fontcolor` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(fontsize) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `fontsize` should be one or @{n_terms} (number of terms in the DAG)."))
	}
	if(!length(shape) %in% c(1, n_terms)) {
		stop_wrap(qq("Length of `shape` should be one or @{n_terms} (number of terms in the DAG)."))
	}

	name = NULL
	if(!is.null(mcols(dag))) {
		meta = mcols(dag)
		if("name" %in% colnames(meta)) {
			name = meta[, "name"]
			name = gsub("\"", "\\\\\"", name)
			name[is.na(name)] = ""
		}
	}

	if(is.null(name)) {
		nodes = paste0(
			qq("  node [fontname=Helvetical]\n"),
			qq("  \"@{dag@terms}\" [color=\"@{color}\", style=\"@{style}\", shape=\"@{shape}\", fontsize=@{fontsize}, fontcolor=\"@{fontcolor}\"];\n", collapse = TRUE)
		)
	} else {
		nodes = paste0(
			qq("  node [fontname=Helvetical]\n"),
			qq("  \"@{dag@terms}\" [color=\"@{color}\", style=\"@{style}\", shape=\"@{shape}\", tooltip=\"@{name}\", fontsize=@{fontsize}, fontcolor=\"@{fontcolor}\"];\n", collapse = TRUE)
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

		if(is.null(edge_color)) {
			if(n_levels <= 8) {
				default_col = c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "#9E9E9E")
				edge_color = structure(default_col[seq_len(n_levels)], names = relation_levels)
			} else if(n_levels <= 26) {
				edge_color = structure(alphabet.colors(n_levels), names = relation_levels)	
			} else {
				edge_color = structure(rand_color(n_levels), names = relation_levels)
			}
		}

		if(length(setdiff(names(edge_color), relation_levels))) {
			stop("names in `edge_color` should cover all relation types.")
		}

		if(is.null(edge_style)) {
			edge_style = structure(rep("solid", n_levels), names = relation_levels)
		}

		if(length(setdiff(names(edge_style), relation_levels))) {
			stop("names in `edge_type` should cover all relation types.")	
		}
	} else {
		v_relations = NULL
	}

	if(is.null(v_relations)) {
		edges = paste0(
			"  edge [dir=\"back\"]\n",
			qq("  \"@{parents}\" -> \"@{children}\";\n", collapse = TRUE)
		)
	} else {
		edges = paste0(
			"  edge [dir=\"back\"]\n",
			qq("  \"@{parents}\" -> \"@{children}\" [color=\"@{edge_color[v_relations]}\", style=\"@{edge_style[v_relations]}\", tooltip=\"@{v_relations}\"];\n", collapse = TRUE)
		)
	}

	DOT = paste0(
		"digraph {\n",
		"  graph [overlap = true]\n",
		"\n",
		nodes,
		"\n",
		edges,
		"}\n"
	)

	DOT
}

#' @param ... Pass to [`DiagrammeR::grViz()`].
#' @rdname dag_viz
#' @details
#' `dag_graphviz()` visualizes the DAG with the **DiagrammeR** package.
#' @export
dag_graphviz = function(dag, color = "black", 
	fontcolor = "black", fontsize = 10, shape = "box",
	edge_color = NULL, edge_style = NULL, ...) {

	if(dag@n_terms > 100) {
		warning("graphviz is only efficient for visualizing small graphs.")
	}
	if(dag@n_terms > 500) {
		stop("Too many terms.")
	}
	check_pkg("DiagrammeR", bioc = FALSE)

	dot = dag_as_DOT(dag, color = color, fontcolor = fontcolor, fontsize = fontsize, shape = shape, edge_color = edge_color, edge_style = edge_style)
	DiagrammeR::grViz(dot, ...)
}

# center at (0, 0)
draw_sector = function(start, end, radius, ...) {
	if(end < start) {
		end = end + 360
	}

	theta = c(0, seq(start, end, by = 0.5), 0)
	rho = c(0, rep(radius, length(theta)-2), 0)
	
	pos = polar2Cartesian(cbind(theta, rho))
	grid.polygon(pos[, 1], pos[, 2], default.units = "native", ...)	
}
