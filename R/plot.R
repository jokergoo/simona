

polar2Cartesian = function(d) {
    theta = d[, 1]/180*pi
    rou = d[, 2]
    x = rou * cos(theta)
    y = rou * sin(theta)
    return(cbind(x, y))
}


#' Visualize the DAG
#' 
#' @param x An `ontology_Dag` object.
#' @param ... Other parguments.
#' 
#' @import grid
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
		gp = gpar(col = add_transparency(col, 0.5)), size = unit(4, "pt"))

	grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(col = "#00000010"))

}


# df1 = term_pos[term_pos$rho == 1, ]
# bg_col = Polychrome::alphabet.colors(nrow(df1))
# bg_col = add_transparency(bg_col, 0.7)

# col = colorRamp2(seq(1, 18,length = 11), rev(brewer.pal(11, "Spectral")))(1:18)

# rho = 1:100
# hc = HilbertCurve(0, 360, level = 7, reference = TRUE, arrow = FALSE, reference_gp = gpar(lty = 1, col = "#EEEEEE"))
# hc_polygon(hc, x1 = df1$theta - df1$width/2, x2 = df1$theta + df1$width/2, gp = gpar(fill = bg_col, col = NA), end_type = "average")
# l = term_pos$rho %in% rho
# hc_points(hc, x1 = term_pos$theta[l], x2 = term_pos$theta[l], np = 1, pch = 16, size = unit(4, "pt"), gp = gpar(col = col[dag_depth(dag)[l]]))


# lt_children = dag@lt_children
# n = dag@n_terms
# df_edge = data.frame(from = rep(seq_len(n), times = sapply(lt_children, length)),
# 	                 to = unlist(lt_children))

# x1 = term_pos[df_edge$from, 1]
# y1 = term_pos[df_edge$from, 2]
# x2 = term_pos[df_edge$to, 1]
# y2 = term_pos[df_edge$to, 2]

# library(spiralize)
# spiral_initialize(xlim = c(0, 360), end = 360*7)
# spiral_track(ylim = c(0, 15))
# # spiral_points(term_pos$theta, term_pos$rho, pch = 16, size = unit(2, "pt"))
# spiral_segments(x1, y1, x2, y2, gp = gpar(col = col[df_edge$from]))
