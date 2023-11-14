

#' Generate a random DAG
#' 
#' @param n_children Range of the numbers of children on a term to sample.
#' @param max Maximal number of terms.
#' @param p Probability to connect to other terms.
#' @param power The power.
#' @param tree Whether the dag is a tree. 
#' @param depth Depth of the tree. Note `tree` and `depth` should be both set. `dag(n_children = nc, tree = TRUE, depth = d)`
#'         is a shortcut of `dag(n_children = c(nc, nc), max = nc^(d+1)/(nc-1), p = 0)`. The setting of `dag(n_children = 2, tree = TRUE, depth = d)`
#'         generates a binary tree with max depth `d` and `2^d` leaf nodes all have the same depth `d`.
#' 
#' @details
#' It first generates a tree growing from the root. On each node, it has
#' a random number of children sampled from `n_children` and each child node
#' can extend with a probability of `1/(height^power)`.
#' 
#' In the second step, each node on the tree has a probability of `p` to link to other
#' terms in the tree.
#' 
#' The growing stops when the number of terms reaches `max` or there is no term can further extend.
#' 
#' @returns An `ontology_DAG` object.
#' @importFrom stats runif
#' @export
#' @examples
#' dag_random()
#' dag_random(p = 0)  # a tree
dag_random = function(n_children = c(2, 5), max = 1000, p = 0.1, power = 0.2,
	tree = FALSE, depth = NULL) {

	lt_children = list()

	n_children = unique(n_children)
	if(tree) {
		if(!is.null(depth)) {
			n_children = n_children[1]
			max = n_children[1]^(depth+1) - 1
		} else {
			stop("If `tree = TRUE`, `depth` should also be set.")
		}
		p = 0
	}

	current_nodes = 1
	current_i = 1
	current_depth = 0
	depth = 1
	while(TRUE) {

		current_depth = current_depth + 1
		current_nodes2 = NULL
		for(ni in current_nodes) {
			if(runif(1) < 1/current_depth^power) {
				if(length(n_children) == 1) {
					nc = n_children
				} else {
					nc = sample(seq(n_children[1], n_children[2]), 1)
				}
				
				if(current_i + nc > max) {
					break
				}
				lt_children[[ni]] = current_i + seq_len(nc)
				depth[lt_children[[ni]]] = current_depth + 1
				for(k in current_i + seq_len(nc)) {
					lt_children[[k]] = integer(0)
				}

				current_i = current_i + nc

				current_nodes2 = c(current_nodes2, lt_children[[ni]])
			}
		}
		current_nodes = current_nodes2
		if(length(current_nodes) == 0) {
			break
		}
	}

	n = length(lt_children)
	if(p > 0 && !tree) {
		for(i in 1:n) {
			if(runif(1) < p) {
				# add a link to nodes with depth >= itself
				l = depth > depth[i] 
				if(length(lt_children[[i]])) {
					l[lt_children[[i]]] = FALSE
				}

				if(any(l)) {
					lt_children[[i]] = c(lt_children[[i]], sample(which(l), 1))
				}
			}
		}
	}

	n = length(lt_children)
	parents = rep(1:n, times = sapply(lt_children, length))
	children = unlist(lt_children)

	parents = as.character(parents)
	children = as.character(children)

	create_ontology_DAG(parents, children, remove_cyclic_paths = TRUE)
}
