
# theta: in degree
.calc_force = function(theta, rho) {
	theta*rho
}

get_force = function(dag, window_size_fun = n_connected_leaves, verbose = verbose) {
	tree = dag_treelize(dag, verbose = FALSE)
	term_pos = exec_under_message_condition({
		cpp_node_pos_in_tree(tree, window_size_fun(tree), 0, 360) ## in polar coordinate
	}, verbose = FALSE)
	theta = term_pos$theta
	depth = dag_depth(dag)
	force = numeric(dag@n_terms)
	for(i in seq_len(dag@n_terms)) {

		pa = dag_parents(dag, i, in_labels = FALSE)
		pa = setdiff(pa, dag@root)
		ch = dag_offspring(dag, i, in_labels = FALSE, include_self = TRUE)

		if(length(ch) == 0) { # leaf
			terms = pa
		} else {
			terms = c(pa, ch)
		}

		if(length(ch)) {
			## ch's other parents that is not i
			cp = setdiff(dag_parents(dag, ch, in_labels = FALSE), i)
			terms = unique(c(terms, cp))
		}

		t = theta[terms]

		t1 = t - theta[i]
		t2 = t + 360 - theta[i]
		t3 = t - 360 - theta[i]

		l = t1 < -180
		t1[l] = t2[l]
		l = t1 > 180
		t1[l] = t3[l]

		t1[abs(t1) < 5] = 0
		force[i] = sum(.calc_force(t1, depth[i]))
		
	}
	force[is.na(force)] = 0
	force
}

reorder_on_circle = function(dag, term_pos, times = 1, window_size_fun = n_connected_leaves, verbose = TRUE) {

	theta = term_pos$x
	width = term_pos$width
	n_off = n_offspring(dag, include_self = TRUE)

	depth = dag_depth(dag)
	total_force = 0
	for(k in seq_len(times)) {
		if(verbose) message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\badjusting branches on the DAG... ", k, ".", appendLF = FALSE)
		force = get_force(dag, window_size_fun = window_size_fun, verbose = verbose)

		if(k == 1) {
			total_force = sum(abs(force))
		} else {
			r = abs(sum(abs(force)) - total_force)/total_force
			if( r > 0.01) {
				total_force = sum(abs(force))
			} else {
				if(verbose) message(" Done.")
				return(theta)
			}
		}

		for(i in seq_len(dag@n_terms)) {

			ch = dag@lt_children[[i]]
			if(length(ch) > 1) {
				dag@lt_children[[i]] = .adjust_children(ch, force, width, depth, n_off)
			}
		}

		tree = dag_treelize(dag, verbose = FALSE)
		term_pos = exec_under_message_condition({
			cpp_node_pos_in_tree(tree, window_size_fun(tree), 0, 360)
		}, verbose = FALSE)
		theta = term_pos$x
		width = term_pos$width
	}
	if(verbose) message(" Done.")

	theta
}

.move_index = function(x, i, j) {  # i: from, j : to
	n = length(x)
	if(i == n && j == 1) {
		x[c(n, seq(1, n - 1))]
	} else if(i == n) {
		x[c(seq(1, j-1), n, seq(j, n-1))]
	} else if(j == 1) {
		x[c(i, seq(1, i-1), seq(i+1, n))]
	} else if(i > j) {
		x[c(seq(1, j-1), i, seq(j, i-1), seq(i+1, n))]
	} else if(i == 1 && j == n) {
		x[c(seq(2, n), 1)]
	} else if(i == 1) {
		x[c(seq(2, j-1), i, seq(j, n))]
	} else if(j == n) {
		x[c(seq(1, i-1), seq(i+1, n), i)]
	} else if(i < j) {
		x[c(seq(1, i-1), seq(i+1, j), i, seq(j+1, n))]
	}
}

.adjust_children = function(ch, force, width, depth, n_off, min_force = 0) {
	f = force[ch]
	n = length(ch)

	od = order(abs(f), decreasing = TRUE)

	ind = seq_len(n)  # the original index order
	cur_pos_left = 1
	cur_pos_right = n

	force1 = f
	ch2 = ch
	for(i in 1:n) {
		
		if(abs(f[od[i]]) < min_force) next
		ind2 = ind
		if(f[od[i]] > 0) {
			ind = .move_index(ind, od[i], cur_pos_right)
			cur_pos_right = cur_pos_right - 1
		} else {
			ind = .move_index(ind, od[i], cur_pos_left)
			cur_pos_left = cur_pos_left - 1
		}
		
		ch2 = ch[ind]
		force2 = f[ind]
	
		offset = (cumsum(width[ch2]) + cumsum(c(0, width[ch2][-n])))/2 - 
		         (cumsum(width[ch]) + cumsum(c(0, width[ch][-n])))/2
		f2 = sapply(seq_along(ch2), function(x) {
			force2[x] + n_off[ch2[x]]*.calc_force(offset[x], depth[ch2[x]])
		})

		if(abs(sum(force1)) > abs(sum(f2))) {
			return(ch[ind2])
		}	
		
		force1 = f2  ## prevous sum of force
	}
	return(ch2)
}
