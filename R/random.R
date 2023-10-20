

#' Generate a random DAG
#' 
#' @param n_children Range of the numbers of children on a term to sample.
#' @param max Maximal number of terms.
#' @param p Probability to connect to other terms.
#' @param power The power.
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
dag_random = function(n_children = c(2, 5), max = 1000, p = 0.1, power = 0.5) {
	lt_children = list()

	current_nodes = 1
	current_i = 1
	depth = 0
	while(TRUE) {
		depth = depth + 1
		current_nodes2 = NULL
		for(ni in current_nodes) {
			if(runif(1) < 1/depth^power) {
				nc = sample(seq(n_children[1], n_children[2]), 1)
				
				if(current_i + nc > max) {
					break
				}
				lt_children[[ni]] = current_i + seq_len(nc)
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
		
	if(p > 0) {
		for(i in 1:n) {
			if(runif(1) < p) {
				lt_children[[i]] = c(lt_children[[i]], sample(setdiff(1:n, c(lt_children[[i]], i)), 1))
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
