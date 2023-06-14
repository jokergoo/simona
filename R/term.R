


#' IC_annotation
#' 
#' @section method:
#' ## IC_annotation
#' Say `k` is the number of items annotated to a term and `N` is the number of items annotated to the root (there is only one
#' global root of the DAG), the IC for the term is calculated as `-log(k/N)`. Due to the DAG structure, if an item is annoated
#' to a term, it is also annotated to all this term's ancestor terms, so `k` is the number of unique terms after merging all
#' its offspring terms. Similarly, all items will be finally annotated to the root term, so `N` is the size of the full set
#' of items that will be annoated to the terms in DAG.
#' 
#' @rdname temp__IC_annotation
IC_annotation = function(dag, use_cache = TRUE) {

	if(is.null(dag@term_env$IC_annotation)) {
		use_cache = FALSE
	} 
	if(!use_cache) {

		n = n_annotations(dag, use_cache)
		n_all_anno = attr(n, "N")
		
		p = n/n_all_anno
		ic = ifelse(n == 0, 0, -log(p))
		
		dag@term_env$IC_annotation = ic
	}
	
	dag@term_env$IC_annotation
}
ADD_IC_METHOD("IC_annotation")



##########################
### IC universe

#' IC_universal
#' 
#' @section method:
#' ## IC_universal
#' 
#' It measures the probability of a term getting full transmission from the root. The probability is calculated
#' recursively. For example, an intermediate term `t` has two parents `parent1` and `parent2`, also assume `parent1` has `k1` children
#' and `parent2` has `k2` children. Assume a parent transmits information equally to its childre, then `parent1` only transmits `1/k1` and
#' `parent2` only transmits `1/k2` to term `t`. Let's say `p1` is the accmulated content from the root and `p2` is the accmulated content
#' from the root, then the content `parent1` transmits to `t` is `p1/k1` and the content `parent2` transmitting to `t` is `p2/k2`. If saying
#' `t` recieves content from both parents, the value is `p1/k1 * p2/k2`. This is the content `t` recieves from all its ancestors.
#' 
#' Let's say `p` is the accumulated content of `t`, the information content is `-log(p)`.
#' 
#' Paper link: <https://doi.org/10.1155/2012/975783>.
#' 
#' @rdname temp__IC_universal
IC_universal = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_universal) || !use_cache) {
		lt_parents = dag@lt_parents

		n_children = dag@term_env$n_children

		n_terms = dag@n_terms

		e = new.env()
		d = rep(NA, n_terms)

		root = dag@root
		d[root] = 0

		depth = dag_depth(dag)

		current_depth = 0
		while(1) {
			current_depth = current_depth + 1
			l = depth == current_depth
			if(sum(l) == 0) {
				break
			}

			for(i in which(l)) {
				p = lt_parents[[i]]
				d[i] = sum(d[p] + log(n_children[p]))
			}
		}
		dag@term_env$IC_universal = d
	}

	dag@term_env$IC_universal
}
ADD_IC_METHOD("IC_universal")


##############################################
### reachability is the number of ways for a node to reach the leaves
reachability = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$reachability) || !use_cache) {
		lt_children = dag@lt_children

		n_terms = dag@n_terms

		reach = rep(NA, n_terms)

		all_leaves = dag@leaves
		reach[all_leaves] = 1

		height = dag_height(dag)

		current_height = 0
		while(1) {
			current_height = current_height + 1
			l = height == current_height
			if(sum(l) == 0) {
				break
			}

			for(i in which(l)) {
				reach[i] = sum(reach[ lt_children[[i]] ])
			}
		}
		dag@term_env$reachability = reach
	}
	dag@term_env$reachability
}

########################################
### Zhang et al

#' IC_Zhang_2006
#' 
#' @section method:
#' ## IC_Zhang_2006
#' 
#' It measures the number of ways from the term to reach a leaf. The information content
#' is calculated as `-log(k/N)` where `k` is the number of ways for term `t` and `N` is the number
#' of ways for the root term, which is the largest number of ways.
#' 
#' Paper link: <https://doi.org/10.1186/1471-2105-7-135>.
#' 
#' @rdname temp__IC_Zhang_2006
IC_Zhang_2006 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Zhang_2006) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@term_env$IC_Zhang_2006 = -log(re/max(re))
	}
	dag@term_env$IC_Zhang_2006
}
ADD_IC_METHOD("IC_Zhang_2006")


########################################
### Seco et al

#' IC_Seco_2004
#' 
#' @section method:
#' ## IC_Seco_2004
#' 
#' Similar as *IC_Zhang_2006*, the information content is `1 - log(k+1)/log(N+1)`.
#' 
#' Paper link: <https://dl.acm.org/doi/10.5555/3000001.3000272>.
#' 
#' @rdname temp__IC_Seco_2004
IC_Seco_2004 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Seco_2004) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@term_env$IC_Seco_2004 = 1 - log(re + 1)/log(max(re) + 1)
	}
	dag@term_env$IC_Seco_2004
}
ADD_IC_METHOD("IC_Seco_2004")


########################################
### Zhou et al

#' IC_Zhou_2008
#' 
#' @section method:
#' ## IC_Zhou_2008
#' 
#' It is calculated as `0.5*IC_Seco + 0.5*log(depth)/log(max_depth)`, where `depth` is the depth of term `t`
#' in the DAG, defined as the maximal distance to root, `max_depth` is the largest depth in the DAG.
#' 
#' Paper link: <https://doi.org/10.1109/FGCNS.2008.16>.
#' 
#' @rdname temp__IC_Zhou_2008
IC_Zhou_2008 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Zhou_2008) || !use_cache) {
		depth = dag_depth(dag, use_cache)
		ic_seco = IC_Seco_2004(dag, use_cache)
		
		sigma = 0.5
		dag@term_env$IC_Zhou_2008 = sigma*ic_seco + (1-sigma)*log(ifelse(depth == 0, 1, depth))/log(max(depth))
	}
	dag@term_env$IC_Zhou_2008
}
ADD_IC_METHOD("IC_Zhou_2008")


########################################
### Seddiqui et al

#' IC_Seddiqui_2010
#' 
#' @section method:
#' ## IC_Seddiqui_2010
#' 
#' It is similar as *IC_Zhou_2008*, defined as `(1-sigma)*IC_Seco + sigma*log((n_parents + n_children + 1)/log((total_edges + 1))`.
#' Here `sigma` is defined as `log(total_edges+1)/(log(total_edges) + log(total_terms))`, where `total_edges` is the number of all relations
#' and `total_terms` is the number of all terms in the DAG. `n_parents` and `n_children` are the number of parents and children of term `t`.
#' 
#' Paper link: <https://dl.acm.org/doi/10.5555/1862330.1862343>.
#' 
#' @rdname temp__IC_Seddiqui_2010
IC_Seddiqui_2010 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Seddiqui_2010) || !use_cache) {
		n_relations = dag@term_env$n_parents + dag@term_env$n_children
		n_edges = sum(dag@term_env$n_parents)
		n_nodes = dag@n_terms

		ic_seco = IC_Seco_2004(dag, use_cache)

		sigma = log(n_edges + 1)/( log(n_edges) + log(n_nodes) )
		dag@term_env$IC_Seddiqui_2010 = (1 - sigma)*ic_seco + sigma*log(n_relations + 1)/log(n_edges + 1)
	}
	dag@term_env$IC_Seddiqui_2010
}
ADD_IC_METHOD("IC_Seddiqui_2010")



#########################################
### Sanchez et al: information pass to leaves

#' IC_Sanchez_2011
#' 
#' @section method:
#' ## IC_Sanchez_2011
#' 
#' It measures average contribution to leaf terms, calculated as `-log(zeta/n_ancestors/total_leaves)`.
#' In the formula, `zeta` is the number of leaves that can be connected from term `t`, since `t`'s ancestors can also
#' reach `t`'s leaves. `zeta` is scaled by `n_ancestors` which is the number of `t`'s ancestor terms. Then it is normalized
#' to the maximal number of leaves in the DAG.
#' 
#' Paper link: <https://doi.org/10.1016/j.knosys.2010.10.001>.
#' 
#' @rdname temp__IC_Sanchez_2011
IC_Sanchez_2011 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Sanchez_2011) || !use_cache) {
		nl = length(dag@leaves)
		n_connected_leaves = n_leaves(dag)

		na = n_ancestors(dag)
		na[na == 0] = 1

		dag@term_env$IC_Sanchez_2011 = -log( n_connected_leaves/na/nl )
	}
	dag@term_env$IC_Sanchez_2011
}
ADD_IC_METHOD("IC_Sanchez_2011")



#########################################
### Meng et al

#' IC_Meng_2012
#' 
#' @section method:
#' ## IC_Meng_2012
#' 
#' It has a complex form which takes depth, 
#' 
#' Paper link: <http://article.nadiapub.com/IJGDC/vol5_no3/6.pdf>.
#' 
#' @rdname temp__IC_Meng_2012
IC_Meng_2012 = function(dag, correct = FALSE, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Meng_2012) || !use_cache) {
		dag@term_env$IC_Meng_2012 = cpp_ic_meng(dag, correct)
	}
	dag@term_env$IC_Meng_2012
}
ADD_IC_METHOD("IC_Meng_2012")



###############################
### totipotency
totipotency = function(dag, use_cache = TRUE) {

	if(is.null(dag@term_env$totipotency) || !use_cache) {
		lt_parents = dag@lt_parents

		n_offspring = n_offspring(dag) + 1
		n_terms = dag@n_terms

		t = rep(NA, n_terms)

		root = dag@root
		t[root] = 1

		depth = dag_depth(dag)

		current_depth = 0
		while(1) {
			current_depth = current_depth + 1
			l = depth == current_depth
			if(sum(l) == 0) {
				break
			}

			for(i in which(l)) {
				p = lt_parents[[i]]
				t[i] = mean(n_offspring[i]/n_offspring[p] * t[p])
			}
		}

		dag@term_env$totipotency = t
	}
	dag@term_env$totipotency
}


#' IC_Wang_2007
#' 
#' @section method:
#' ## IC_Wang_2007
#' 
#' It first calculates a "S" value between from `t`'s ancestor term to `t`, defined as the maximal distance by taking semantic weight.
#' Then the information content of `t` is the sum of `S` values of all its ancestor terms.
#' 
#' Paper link: <https://doi.org/10.1093/bioinformatics/btm087>.
#' 
#' @rdname temp__IC_Wang_2007
IC_Wang_2007 = function(dag, contribution_factor = c("isa" = 0.8, "part of" = 0.6), use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Wang_2007) || !use_cache) {
		if(length(dag@lt_children_relations) == 0) {
			stop("`relations` is not set when creating the ontology_DAG object.")
		}
		relation_levels = attr(dag@lt_children_relations, "levels")
		if(is.null(names(contribution_factor))) {
			stop("`contribution_factor` should be a named numeric vector where names should correspond to all relations.")
		}
		if(length(setdiff(relation_levels, names(contribution_factor)))) {
			stop("Contribution factor should be provided for all relations.")
		}
		dag@term_env$IC_Wang_2007 = cpp_ic_wang(dag, unname(contribution_factor[relation_levels]))
	}
	dag@term_env$IC_Wang_2007
}
ADD_IC_METHOD("IC_Wang_2007")


