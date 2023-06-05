

polar2Cartesian = function(d) {
    theta = d[, 1]/180*pi
    rou = d[, 2]
    x = rou * cos(theta)
    y = rou * sin(theta)
    return(cbind(x, y))
}


#' Visualize the DAG
#' 
#' @param x
#' 
plot.ontology_DAG = function(x, ...) {

	dag = x
	term_pos = cpp_term_pos_on_circle(dag) ## in polar coordinate

	level1 = setdiff(unique(term_pos$level1_group), 0)
	all_colors = Polychrome::alphabet.colors(length(level1))
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
		gp = gpar(col = add_transparency(col, 0.5)), size = unit(sqrt(n_children)+2, "pt"))

	grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(col = "#00000002"))
}
