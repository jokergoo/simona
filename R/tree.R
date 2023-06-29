
dag_treelize = function(dag, based_on = "depth") {

	if(dag_is_tree(dag)) {
		return(dag)
	}

	n = dag@n_terms
	lt_children = dag@lt_children
	lt_children2 = rep(list(integer(0)), n)

	if(based_on == "depth") {
		depth = dag_depth(dag)
	} else {
		depth = cpp_dag_dist_from_root(dag)
	}

	current_depth = 0
	current_term = dag@root
	l_visited = rep(FALSE, n)
	while(1) {
		if(length(current_term) == 0) {
			break
		}
		current_depth = current_depth + 1
		current_term2 = integer(0)
		for(t in current_term) {
			children = lt_children[[t]]
			for(ch in children) {
				if(depth[ch] == current_depth && !l_visited[ch]) {
					lt_children2[[t]] = c(lt_children2[[t]], ch)
					l_visited[ch] = TRUE
					current_term2 = c(current_term2, ch)
				}
			}
		}
		current_term = current_term2
	}

	n = length(lt_children2)
	parents = rep(seq_len(n), times = sapply(lt_children2, length))
	children = unlist(lt_children2)

	tree = create_ontology_DAG(dag@terms[parents], dag@terms[children])
	tree@annotation = dag@annotation

	tree
}


dag_as_dendrogram = function(dag) {

	if(!dag_is_tree(dag)) {
		stop("dag is not a tree.")
	}

	add_dend = function(node, lt_children, height) {
		if(length(lt_children[[node]]) == 0) {			
			dend = node
			attributes(dend) = list(members = 1, 
				                    height = height, 
				                    label = dag@terms[node], 
				                    term_id = node, 
				                    n_nodes = 1,
				                    max_dist_to_leaf = 0,
				                    leaf = TRUE, 
				                    class = "dendrogram")
		} else {
			dend = lapply(lt_children[[node]], add_dend, lt_children, height - 1)
			attributes(dend) = list(members = length(unlist(dend)), 
				                    height = height, 
				                    label = dag@terms[node], 
				                    term_id = node, 
				                    n_nodes = 1 + sum(sapply(dend, function(x) attr(x, "n_nodes"))),
				                    max_dist_to_leaf = 1 + max(sapply(dend, function(x) attr(x, "max_dist_to_leaf"))),
				                    leaf = FALSE, 
				                    class = "dendrogram")
		}
		dend
	}

	max_depth = max(dag_depth(dag))
	dend = add_dend(dag@root, dag@lt_children, max_depth)

	dend = dend_set_midpoint(dend)
	class(dend) = c("ontology_tree", "dendrogram")

	dend
}


print.ontology_tree = function(x, ...) {
	cat("This is a normal `dendrogram` object, but with several more attributes attached on each node:\n")
	cat("  - label: term name\n")
	cat("  - term_id: the integer index of the term\n")
	cat("  - n_nodes: number of offspring terms, including itself\n")
	cat("  - members: this is the standard attribute, which is the number of leaves the term can reach\n")
	cat("  - max_dist_to_leaf: max distance to the term's leaves\n")
	cat("\n")
	cat("Two special functions can be applied:\n")
	cat("  - dend_sort(): sort the dendrogram based on an attributes. The following attributes can be used:\n")
	cat("                   'members', 'n_nodes', 'max_dist_to_leaf'\n")
	cat("  - dend_tpl(): topologically sort the nodes on the dendrogram via depth first search.")
	cat(" \n")
	cat(" \n")
	getFromNamespace("print.dendrogram", "stats")(x)
}

dend_add_x = function(dend) {

	env = new.env(parent = emptyenv())
	i_leaf = 0

	add_x = function(dend) { # current sub-dend
		if(attr(dend, "leaf")) {
			i_leaf <<- i_leaf + 1
			x = i_leaf
			attr(dend, "x") = x
			attr(dend, "width") = 0
		} else {
			k = length(dend) # number of sub-dend
			for(i in seq_len(k)) {
				dend[[i]] = add_x(dend[[i]])
			}
			x = (attr(dend[[1]], "x") + attr(dend[[k]], "x"))/2
			attr(dend, "x") = x
			attr(dend, "width") = attr(dend[[k]], "x") - attr(dend[[1]], "x")
		}
		dend
	}

	add_x(dend)
}

dend_set_midpoint = function(dend) {

	dend = dend_add_x(dend)
	
	add_midpoint = function(dend) {
		if(attr(dend, "leaf")) {
			attr(dend, "midpoint") = 0
		} else {
			k = length(dend)
			for(i in seq_len(k)) {
				dend[[i]] = add_midpoint(dend[[i]])
			}
			midpoint = attr(dend[[1]], "midpoint") + attr(dend, "width")/2
			
			attr(dend, "midpoint") = midpoint
		}
		dend
	}

	add_midpoint(dend)
}


dend_sort = function(dend, by = "n_nodes", decreasing = TRUE) {

	if(!by %in% names(attributes(dend))) {
		stop("cannot find ...")
	}

	reorder = function(dend) {
		k = length(dend)
		if(k > 1) {
			v = sapply(dend, function(x) attr(x, by))
			attr = attributes(dend)

			dend = dend[order(v, decreasing = decreasing)]
			attributes(dend) = attr
			attr(dend, "width") = attr(dend[[k]], "x") - attr(dend[[1]], "x")

			for(i in seq_len(k)) {
				dend[[i]] = reorder(dend[[i]])
			}
		}
		dend
	}

	dend = reorder(dend)
	dend_set_midpoint(dend)
}

dend_tpl_sort = function(dend, in_labels = TRUE) {

	env = new.env(parent = emptyenv())
	n_nodes = attr(dend, "n_nodes")

	if(in_labels) {
		env$tpl_sorted = character(n_nodes)
	} else {
		env$tpl_sorted = integer(n_nodes)
	}
	env$tpl_index = 0

	add_tpl_index = function(dend) {
		if(!attr(dend, "leaf")) {
			for(i in seq_along(dend)) {
				add_tpl_index(dend[[i]])
			}
		}

		env$tpl_index = env$tpl_index + 1
		if(in_labels) {
			env$tpl_sorted[env$tpl_index] = attr(dend, "label")
		} else {
			env$tpl_sorted[env$tpl_index] = attr(dend, "term_id")
		}
	}

	add_tpl_index(dend)

	env$tpl_sorted
}

dend_tpl_rank = function(dend) {

	env = new.env(parent = emptyenv())
	n_nodes = attr(dend, "n_nodes")

	env$tpl_rank = integer(n_nodes)
	env$tpl_index = 0

	add_tpl_rank = function(dend) {
		if(!attr(dend, "leaf")) {
			for(i in seq_along(dend)) {
				add_tpl_rank(dend[[i]])
			}
		}

		env$tpl_index = env$tpl_index + 1
		env$tpl_rank[attr(dend, "term_id")] = env$tpl_index
	}

	add_tpl_rank(dend)

	env$tpl_rank
}


dend_as_data_frame = function(dend) {
	env = new.env(parent = emptyenv())
	n_nodes = attr(dend, "n_nodes")

	env$list = vector("list", n_nodes)
	message("extracting attributes on each node...")
	add_node_attr = function(dend) {
		if(!attr(dend, "leaf")) {
			for(i in seq_along(dend)) {
				add_node_attr(dend[[i]])
			}
		}

		at = attributes(dend)
		at$class = NULL

		id = at$term_id
		env$list[[id]] = at
	}

	add_node_attr(dend)

	df = do.call(rbind, lapply(env$list, data.frame))
	df = df[, setdiff(colnames(df), c("class"))]

	message("extracting parent id...")
	parent_id = integer(n_nodes)
	add_parent_id = function(dend, parent) {
		if(!attr(dend, "leaf")) {
			for(i in seq_along(dend)) {
				add_parent_id(dend[[i]], attr(dend, "term_id"))
			}
		}

		parent_id[attr(dend, "term_id")] <<- parent
	}

	parent_id[attr(dend, "term_id")] = NA_integer_
	for(i in seq_along(dend)) {
		add_parent_id(dend[[i]], attr(dend, "term_id"))
	}

	message("extracting level1 ancestor id...")
	level1_ancestor = integer(n_nodes)
	level1_ancestor[attr(dend, "term_id")] = NA_integer_
	for(i in seq_along(dend)) {
		id = attr(dend[[i]], "term_id")
		dendrapply(dend[[i]], function(d) {
			level1_ancestor[attr(d, "term_id")] <<- id
		})
	}

	df$parent_id = parent_id
	df$level1_ancestor = level1_ancestor

	message("extracting tpl rank")
	df$tpl_rank = dend_tpl_rank(dend)

	df
}

#' @import HilbertCurve
dend_hc = function(df) {
	n = nrow(df)
	df = df[order(df$tpl_rank), ]
	n_bg = table(df$level1_ancestor)
	
	bg_col = Polychrome::alphabet.colors(length(n_bg))
	bg_col = add_transparency(bg_col, 0.7)

	col = colorRamp2(seq(1, 19,length = 11), brewer.pal(11, "Spectral"))(1:19)

	hc = HilbertCurve(0, n, level = 7, reference = TRUE, arrow = FALSE, reference_gp = gpar(lty = 1, col = "#EEEEEE"))
	
	hc_polygon(hc, x1 = c(0, cumsum(n_bg[-length(n_bg)])), x2 = cumsum(n_bg), gp = gpar(fill = NA, col = "black"), end_type = "average")
	
	hc_points(hc, x1 = 1:n, x2 = 1:n, np = 1, pch = 16, size = unit(2, "pt"), gp = gpar(col = col[df$height+1]))
}

#' @import spiralize
dend_spiralize = function(dend) {
	spiral_initialize(xlim = c(0, nobs(dend)), start = 45, end = 360*7)
	spiral_track(height = 1, background_gp = gpar(fill = "#EEEEEE", col = NA))
	spiral_dendrogram(dend)
}


#' @importFrom stats dendrapply
fill_treemap = function(dend) {

	term_coor = matrix(-1, nrow = attr(dend, "n_nodes"), ncol = 4)
	grid.newpage()

	fill_rect = function(dend, bottom_x, bottom_y, top_x, top_y) {
		if(attr(dend, "leaf")) {
			leaf_coor[attr(dend, "term_id"), ] <<- 	c(bottom_x, bottom_y, top_x, top_y)
		
			return( c((bottom_x + top_x)/2, (bottom_y + top_y)/2))
		}

		children_size = sapply(dend, function(x) {
			attr(x, "n_nodes")
		})

		od = order(children_size, decreasing = TRUE)
		children_size = children_size[od]
		coor = split_rect(children_size, bottom_x, bottom_y, top_x, top_y, od)
		for(i in seq_along(children_size)) {
			fill_rect(dend[[ coor[i, 5] ]], coor[i, 1], coor[i, 2], coor[i, 3], coor[i, 4])
		}

		term_coor[ coor[, 5], ] = coor
	}

	fill_rect(dend, 0, 0, 1, 1)

	leave = unlist(dendrapply(dend, function(d) {
		if(attr(d, "leaf")) {
			attr(d, "term_id")
		} else {
			NA
		}
	}))
	leaf_coor = coor[leave, ]
	grid.rect( (leaf_coor[,1]+leaf_coor[,3])/2, (leaf_coor[,2]+leaf_coor[,4])/2, (leaf_coor[,3] - leaf_coor[,1]), leaf_coor[,4] - leaf_coor[,2], gp = gpar(col = "white", fill = rand_color(nrow(leaf_coor))))

	coor
}

split_rect = function(value, bottom_x, bottom_y, top_x, top_y, od) {

	n = length(value)
	w = top_x - bottom_x
	h = top_y - bottom_y

	if(n == 1) {

		return(matrix(c(bottom_x, bottom_y, top_x, top_y, od), nrow = 1))
	}
	
	if(value[1] > ceiling(sum(value)/2)) {
		l_group1 = rep(FALSE, n)
		l_group1[1] = TRUE
	} else {
		partition = adagio::bpp_approx(value, ceiling(sum(value)/2), method = "bestfit")
		l_group1 = partition$xbins == 1
	}

	p = sum(value[l_group1])/sum(value)

	if(w/h < 1) {
		rbind(
			split_rect(value[l_group1],  bottom_x, bottom_y, top_x, bottom_y + h*p, od[l_group1]),
			split_rect(value[!l_group1], bottom_x, bottom_y + h*p, top_x, top_y, od[!l_group1])
		)
	} else {
		rbind(
			split_rect(value[l_group1],  bottom_x, bottom_y, bottom_x + w*p, top_y, od[l_group1]),
			split_rect(value[!l_group1], bottom_x + w*p, bottom_y, top_x, top_y, od[!l_group1])
		)
	}
}

