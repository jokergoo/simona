

#' Generate a random DAG
#' 
#' @param n_children Number of children of a term. The value can also be a vector of
#'     length two representing the range of the number of child terms.
#' @param p_stop The probability of a term to stop growing.
#' @param max Maximal number of terms.
#' @param verbose Whether to print messages.
#' 
#' @details
#' `dag_random_tree()` generates a random DAG tree from the root term. In a certain step of
#' the growing, let's denote the set of all leaf terms as L, then in the next round of growing,
#' `floor(length(L)*p_stop)` leaf terms stop growing, and for the remaining leaf terms that
#' continue to grow, each term will add child terms with number in uniformly sampled within `[ n_children[1], n_children[2] ]`.
#' The growing stops when the total number of terms in the DAG exceeds `max`.
#' 
#' @returns An `ontology_DAG` object.
#' @rdname dag_random
#' @export
#' @examples
#' tree = dag_random_tree()
#' dag = dag_random()
dag_random_tree = function(n_children = 2, p_stop = 0, max = 2^10-1,
	verbose = simona_opt$verbose) {

	lt_children = list()

	n_children = unique(range(n_children))

	current_nodes = 1
	current_i = 1
	while(TRUE) {

		current_nodes2 = NULL

		k = length(current_nodes)

		ind_stop = sample(k, floor(k*p_stop))
		ind_continue = setdiff(sample(k, k), ind_stop)  # randomly permutated

		# number of children for each "continued" node
		if(length(n_children) == 1) {
			nc = rep(n_children, length(ind_continue))
		} else {
			nc = sapply(ind_continue, function(x) {
				sample(seq(n_children[1], n_children[2]), 1)
			})
		}

		# compare to max
		ll = cumsum(nc) <= max - current_i

		ind_stop = c(ind_stop, ind_continue[!ll])
		ind_continue = ind_continue[ll]
		nc = nc[ll]
		
		current_nodes2 = current_nodes[ind_stop]

		if(length(ind_continue) == 0) {
			break
		}

		for(i in seq_along(ind_continue)) {
			ni = current_nodes[ ind_continue[i] ]

			
			lt_children[[ni]] = current_i + seq_len(nc[i])

			# initialize new children of lt_children[[ni]]
			for(k in lt_children[[ni]]) {
				lt_children[[k]] = integer(0)
			}

			current_i = current_i + nc[i]

			current_nodes2 = c(current_nodes2, lt_children[[ni]])
		}
		current_nodes = current_nodes2
		
	}

	n = length(lt_children)
	parents = rep(1:n, times = sapply(lt_children, length))
	children = unlist(lt_children)

	parents = as.character(parents)
	children = as.character(children)

	create_ontology_DAG(parents, children, source = "dag_random_tree", verbose = verbose)
}


#' @param dag An `ontology_DAG` object.
#' @param p_add The probability to add children on each term.
#' @param new_children The number or range of numbers of new children if a term is selected to add more children.
#' @param add_random_children_fun A function to randomly add children from the DAG.
#' 
#' @details
#' `dag_add_random_children()` adds more links in a DAG. Each term is associated with a probability `p_add`
#' to add new links where the term, if it is selected, is as a parent term, linking to other terms in the DAG.
#' The number of new child terms is controlled by `new_children` which can be a single number of a range. By default,
#' new child terms of a term `t` are randomly selected from other terms that are lower than the term `t`
#' (check the function `simona:::add_random_children`). The way how to randomly select new child terms for `t`
#' can be controlled by a self-defined function for the `add_random_children_fun` argument.
#' 
#' @rdname dag_random
#' @importFrom stats runif
#' @export
dag_add_random_children = function(dag, p_add = 0.1, new_children = c(1, 4),
	add_random_children_fun = NULL, verbose = simona_opt$verbose) {

	if(p_add < 0) {
		return(dag)
	}

	if(length(new_children) == 1) {
		new_children = c(new_children, new_children)
	}

	depth = dag_depth(dag)
	lt_children = dag@lt_children
	n_terms = dag_n_terms(dag)

	if(!is.null(add_random_children_fun)) {
		for(i in 1:n_terms) {
			if(runif(1) < p_add) {
				lt_children[[i]] = c(lt_children[[i]], add_random_children_fun(dag, i))
			}
		}
	} else {
		for(i in 1:n_terms) {
			if(runif(1) < p_add) {
				lt_children[[i]] = c(lt_children[[i]], add_random_children(dag, i, new_children))
			}
		}
	}

	lt_children = lapply(lt_children, unname)

	parents = rep(1:n_terms, times = sapply(lt_children, length))
	children = unlist(lt_children)

	parents = dag@terms[parents]
	children = dag@terms[children]

	dag2 = create_ontology_DAG(parents, children, verbose = verbose)
	if(has_annotation(dag)) {
		dag2@annotation = dag@annotation
	}
	mcols(dag2) = mcols(dag)

	dag2@source = "dag_add_random_children"

	dag2
}


add_random_children = function(dag, i, new_children = c(1, 4)) {
	# add a link to nodes with depth > itself
	depth = dag_depth(dag)
	l = depth > depth[i] 
	if(length(dag@lt_children[[i]])) {
		l[dag@lt_children[[i]]] = FALSE
	}

	candidates = which(l)
	n_candidates = length(candidates)
	if(n_candidates) {
		if(n_candidates < new_children[1]) {
			integer(0)
		} else {
			sample(candidates, min(n_candidates, sample(seq(new_children[1], new_children[2]), 1)))
		}
	} else {
		integer(0)
	}	
}


#' @details
#' `dag_random()`: it simply wraps `dag_random_tree()` and `dag_add_random_children()`.
#' 
#' @rdname dag_random
#' @export
dag_random = function(n_children = 2, p_stop = 0, max = 2^10-1,
	p_add = 0.1, new_children = c(1, 4), verbose = simona_opt$verbose) {

	tree = dag_random_tree(n_children = n_children, p_stop = p_stop, max = max, verbose = verbose)
	dag_add_random_children(tree, p_add = p_add, new_children = new_children, verbose = verbose)
}

