

#' IC_offspring
#' 
#' @section Methods:
#' ## IC_offspring
#' Denote `k` as the number of offspring terms plus the term itself and `N` is such value for root (or the total number of terms in the DAG), the information
#' content is calculated as:
#' 
#' ```
#' IC = -log(k/N)
#' ```
#' 
#' @rdname temp__IC_offspring
IC_offspring = function(dag, use_cache = simone_opt$use_cache) {

	if(is.null(dag@term_env$IC_offspring)) {
		use_cache = FALSE
	} 
	if(!use_cache) {

		n = n_offspring(dag, include_self = TRUE)
		
		p = n/max(n)
		ic = -log(p)
		dag@term_env$IC_offspring = ic
		
	}
	
	
	dag@term_env$IC_offspring
}
ADD_IC_METHOD("IC_offspring", "use_cache")


#' IC_depth
#' 
#' @section Methods:
#' ## IC_height
#' For a term `t` in the DAG, denote `d` as the maximal distance from root (i.e. the depth) and `h` as the maximal distance to leaves (i.e. the height),
#' the relative position `p` on the longest path from root to leaves via term `t` is calculated as:
#' 
#' ```
#' p = (h + 1)/(h + d + 1)
#' ```
#' 
#' In the formula where 1 is added gets rid of `p = 0`. Then the information content is:
#' 
#' ```
#' IC = -log(p) 
#'    = -log((h+1)/(h+d+1))
#' ```
#' 
#' @rdname temp__IC_height
IC_height = function(dag, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_height)) {
		use_cache = FALSE
	} 
	if(!use_cache) {

		ic = -log((dag_height(dag) + 1)/(dag_depth(dag) + dag_height(dag) + 1))
		dag@term_env$IC_height = ic
		
	}
	
	
	dag@term_env$IC_height
}
ADD_IC_METHOD("IC_height", "use_cache")


#' IC_annotation
#' 
#' @section Methods:
#' ## IC_annotation
#' Denote `k` as the number of items annotated to a term `t`, and `N` is the number of items annotated to the root (which is
#' the total number of items annotated to the DAG), IC for term `t` is calculated as:
#' 
#' ```
#' IC = -log(k/N)
#' ```
#' 
#' In current implementations in other tools, there is an inconsistency of defining `k` and `N`. 
#' Please see [`n_annotations()`] for explanation.
#' 
#' `NA` is assigned to the terms with no item annotated.
#' 
#' @rdname temp__IC_annotation
IC_annotation = function(dag, uniquify = simone_opt$anno_uniquify, use_cache = simone_opt$use_cache) {

	if(!uniquify && is.null(dag@term_env$IC_annotation)) {
		use_cache = FALSE
	} else if(uniquify && is.null(dag@term_env$IC_annotations_unique)) {
		use_cache = FALSE
	}
	if(!use_cache) {

		n = n_annotations(dag, uniquify = uniquify, use_cache = use_cache)
		
		p = n/max(n)
		ic = ifelse(n == 0, NA_real_, -log(p))
		
		if(uniquify) {
			dag@term_env$IC_annotation_unique = ic
		} else {
			dag@term_env$IC_annotation = ic
		}
	}
	
	if(uniquify) {
		dag@term_env$IC_annotation_unique
	} else {
		dag@term_env$IC_annotation
	}
}
ADD_IC_METHOD("IC_annotation", c("uniquify", "use_cache"))



##########################
### IC universal

#' IC_universal
#' 
#' @section Methods:
#' ## IC_universal
#' 
#' It measures the probability of a term getting full transmission from the root. Each term is associated with a p-value and the root has
#' the p-value of 1.
#' 
#' For example, an intermediate term `t` has two parent terms `parent1` and `parent2`, also assume `parent1` has `k1` children
#' and `parent2` has `k2` children, assume a parent transmits information equally to all its children, then `parent1` only transmits `1/k1` and
#' `parent2` only transmits `1/k2` of its content to term `t`, or the probability of a parent to reach `t` is `1/k1` or `1/k2`. 
#' Let's say `p1` and `p2` are the accmulated contents from the root for `parnet1` and `parent2` respectively (or the probability 
#' of the two parent terms getting full transmission from root), then the probability of reaching `t` via a full transmission graph from `parent1`
#' is the multiplication of `p1` and `1/k1`, which is `p1/k1`, and it is similar for `p2/k2`. Then, for term `t`, if getting transmitted from `parent1` and
#' `parent2` are independent, the probability of `t` (denoted as `p_t`) to get transmitted from both parents is:
#' 
#' ```
#' p_t = (p1/k1) * (p2/k2)
#' ```
#' 
#' Since the two parents are the full set of `t`'s parents, `p_t` is the probability of `t` getting full transmission from root. And the final
#' information content is:
#' 
#' ```
#' IC = -log(p_t)
#' ```
#' 
#' Paper link: <https://doi.org/10.1155/2012/975783>.
#' 
#' @rdname temp__IC_universal
IC_universal = function(dag, use_cache = simone_opt$use_cache) {
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
ADD_IC_METHOD("IC_universal", "use_cache")


##############################################
### reachability is the number of ways for a node to reach the leaves
reachability = function(dag, use_cache = simone_opt$use_cache) {
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
#' @section Methods:
#' ## IC_Zhang_2006
#' 
#' It measures the number of ways from a term to reach leaf terms. E.g. in the following DAG:
#' 
#' ```
#'      a      upstream
#'     /|\
#'    b | c
#'      |/
#'      d      downstream
#' ```
#' 
#' term `a` has three ways to reach leaf, which are `a->b`, `a->d` and `a->c->d`.
#' 
#' Let's denote `k` as the number of ways for term `t` to reach leaves and `N` as the maximal value of `k` which
#' should be associated with the root term, the information content is calculated as 
#' 
#' ```
#' IC = -log(k/N) 
#'    = log(N) - log(k)
#' ``` 
#' 
#' Paper link: <https://doi.org/10.1186/1471-2105-7-135>.
#' 
#' @rdname temp__IC_Zhang_2006
IC_Zhang_2006 = function(dag, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Zhang_2006) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@term_env$IC_Zhang_2006 = -log(re/max(re))
	}
	dag@term_env$IC_Zhang_2006
}
ADD_IC_METHOD("IC_Zhang_2006", "use_cache")


########################################
### Seco et al

#' IC_Seco_2004
#' 
#' @section Methods:
#' ## IC_Seco_2004
#' 
#' It is based on the number of offspring terms of term `t`.
#' The information content is calculated as:
#' 
#' ```
#' IC = 1 - log(k+1)/log(N)
#' ```
#' 
#' where `k` is the number of offspring terms of `t`, or you can think `k+1` is the number of `t`'s offspring terms plus itself.
#' `N` is the total number of terms in the DAG.
#' 
#' Paper link: <https://dl.acm.org/doi/10.5555/3000001.3000272>.
#' 
#' @rdname temp__IC_Seco_2004
IC_Seco_2004 = function(dag, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Seco_2004) || !use_cache) {
		n = n_offspring(dag, include_self = TRUE)
		dag@term_env$IC_Seco_2004 = 1 - log(n)/log(dag@n_terms) # dag@n_terms == max(n)
	}
	dag@term_env$IC_Seco_2004
}
ADD_IC_METHOD("IC_Seco_2004", "use_cache")


########################################
### Zhou et al

#' IC_Zhou_2008
#' 
#' @section Methods:
#' ## IC_Zhou_2008
#' 
#' It is a correction of *IC_Seco_2004* which considers the depth of a term in the DAG.
#' The information content is calculated as:
#' 
#' ```
#' IC = 0.5*IC_Seco + 0.5*log(depth)/log(max_depth)
#' ```
#' 
#' where `depth` is the depth of term `t` in the DAG, defined as the maximal distance from root. `max_depth` is the largest depth in the DAG.
#' So IC is composed with two parts: the numbers of offspring terms and positions in the DAG.
#' 
#' Paper link: <https://doi.org/10.1109/FGCNS.2008.16>.
#' 
#' @rdname temp__IC_Zhou_2008
IC_Zhou_2008 = function(dag, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Zhou_2008) || !use_cache) {
		depth = dag_depth(dag, use_cache = use_cache)
		ic_seco = IC_Seco_2004(dag, use_cache = use_cache)
		
		sigma = 0.5
		dag@term_env$IC_Zhou_2008 = sigma*ic_seco + (1-sigma)*log(ifelse(depth == 0, 1, depth))/log(max(depth))
	}
	dag@term_env$IC_Zhou_2008
}
ADD_IC_METHOD("IC_Zhou_2008", "use_cache")


# ########################################
# ### Seddiqui et al

#' IC_Seddiqui_2010
#' 
#' @section Methods:
#' ## IC_Seddiqui_2010
#' 
#' It is also a correction to *IC_Seco_2004*, but considers number of relations connecting a term (i.e. number of parent terms and child terms).
#' The information content is defined as:
#' 
#' ```
#' (1-sigma)*IC_Seco + sigma*log((n_parents + n_children + 1)/log((total_edges + 1))
#' ```
#' 
#' where `n_parents` and `n_children` are the numbers of parents and children of term `t`. The tuning factor `sigma` is defined as 
#' 
#' ```
#' sigma = log(total_edges+1)/(log(total_edges) + log(total_terms))
#' ```
#' 
#' where `total_edges` is the number of all relations (all parent-child relations)
#' and `total_terms` is the number of all terms in the DAG. 
#' 
#' Paper link: <https://dl.acm.org/doi/10.5555/1862330.1862343>.
#' 
#' @rdname temp__IC_Seddiqui_2010
IC_Seddiqui_2010 = function(dag, use_cache = simone_opt$use_cache) {
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
ADD_IC_METHOD("IC_Seddiqui_2010", "use_cache")



#########################################
### Sanchez et al: information pass to leaves

#' IC_Sanchez_2011
#' 
#' @section Methods:
#' ## IC_Sanchez_2011
#' 
#' It measures average contribution of term `t` on leaf terms. First denote `zeta` as the number of leaf terms that
#' can be reached from term `t` (or `t`'s offspring that are leaves.). Since all `t`'s ancestors can also
#' reach `t`'s leaves, the contribution of `t` on leaf terms is scaled by `n_ancestors` which is the number of `t`'s ancestor terms.
#' The final information content is normalized by the total number of leaves in the DAG, which is the possible maximal value of `zeta`.
#' The complete definition of information content is:
#' 
#' ```
#' IC = -log( (zeta/n_ancestor) / n_all_leaves)
#' ```
#' 
#' Paper link: <https://doi.org/10.1016/j.knosys.2010.10.001>.
#' 
#' @rdname temp__IC_Sanchez_2011
IC_Sanchez_2011 = function(dag, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Sanchez_2011) || !use_cache) {
		nl = length(dag@leaves)
		n_connected_leaves = n_connected_leaves(dag)

		na = n_ancestors(dag)
		na[na == 0] = 1

		dag@term_env$IC_Sanchez_2011 = -log( (n_connected_leaves+1)/na/nl )
	}
	dag@term_env$IC_Sanchez_2011
}
ADD_IC_METHOD("IC_Sanchez_2011", "use_cache")



#########################################
### Meng et al

#' IC_Meng_2012
#' 
#' @section Methods:
#' ## IC_Meng_2012
#' 
#' It has a complex form which takes account of the term depth and the downstream of the term.
#' The first factor is calculated as:
#' 
#' ```
#' f1 = log(depth)/long(max_depth)
#' ``` 
#' 
#' The second factor is calculated as:
#' 
#' ```
#' f1 = 1 - log(1 + sum_{x => t's offspring}(1/depth_x))/log(total_terms)
#' ```
#' 
#' In the equation, the summation goes over `t`'s offspring terms.
#' 
#' The final information content is the multiplication of `f1` and `f2`:
#' 
#' ```
#' IC = f1 * f2
#' ```
#' 
#' Paper link: <http://article.nadiapub.com/IJGDC/vol5_no3/6.pdf>.
#' 
#' There is one parameter `correct`. If it is set to `TRUE`, the first factor `f1` is calculated as:
#' 
#' ```
#' f1 = log(depth + 1)/long(max_depth + 1)
#' ```
#' 
#' `correct` can be set as:
#' 
#' ```
#' term_IC(dag, method = "IC_Meng_2012", control = list(correct = TRUE))
#' ```
#' 
#' @rdname temp__IC_Meng_2012
IC_Meng_2012 = function(dag, correct = FALSE, use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Meng_2012) || !use_cache) {
		dag@term_env$IC_Meng_2012 = cpp_ic_meng(dag, correct)
	}
	dag@term_env$IC_Meng_2012
}
ADD_IC_METHOD("IC_Meng_2012", c("correct", "use_cache"))



###############################
### totipotency
totipotency = function(dag, use_cache = simone_opt$use_cache) {

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
#' @section Methods:
#' ## IC_Wang_2007
#' 
#' Each relation is weighted by a value less than 1 based on the semantic relation, i.e. 0.8 for "is_a" and 0.6 for "part_of".
#' For a term `t` and one of its ancestor term `a`, it first calculates an "S-value" which corresponds to a path from `a` to `t` where
#' the accumulated multiplication of weights along the path reaches maximal:
#' 
#' ```
#' S(a->t) = max(prod(w))
#' ```
#' 
#' Here `max` goes over all possible paths from `a` to `t`, and `prod()` multiplies edge weights in a certain path.
#' 
#' The formula can be transformed as (we simply rewrite `S(a->t)` to `S`):
#' 
#' ```
#' 1/S = min(prod(1/w))
#' log(1/S) = log(min(prod(1/w)))
#'          = min(sum(log(1/w)))
#' ```
#' 
#' Since `w < 1`, `log(1/w)` is positive. According to the equation, the path (`a->...->t`) is actually the shortest path from `a` to `t` by taking
#' `log(1/w)` as the weight, and `log(1/S)` is the weighted shortest distance.
#' 
#' If `S(a->t)` can be thought as the maximal semantic contribution from `a` to `t`, the information content is calculated
#' as the sum from all `t`'s ancestors (including `t` itself):
#' 
#' ```
#' IC = sum_{a => t's ancestors + t}(S(a->t))
#' ```
#' 
#' Paper link: <https://doi.org/10.1093/bioinformatics/btm087>.
#' 
#' The contribution of different semantic relations can be set with the `contribution_factor` parameter. The value should be a named numeric
#' vector where names should cover the relations defined in `relations` set in [`create_ontology_DAG()`]. For example, if there are two relations
#' "relation_a" and "relation_b" set in the DAG, the value for `contribution_factor` can be set as:
#' 
#' ```
#' term_IC(dag, method = "IC_Wang", 
#'     control = list(contribution_factor = c("relation_a" = 0.8, "relation_b" = 0.6)))
#' ```
#' 
#' @rdname temp__IC_Wang_2007
IC_Wang_2007 = function(dag, contribution_factor = c("is_a" = 0.8, "part_of" = 0.6), use_cache = simone_opt$use_cache) {
	if(is.null(dag@term_env$IC_Wang_2007) || !use_cache) {
		if(length(dag@lt_children_relations) == 0) {
			stop("`relations` is not set when creating the ontology_DAG object.")
		}
		relation_levels = attr(dag@lt_children_relations, "levels")
		if(is.null(names(contribution_factor))) {
			stop("`contribution_factor` should be a named numeric vector where names should correspond to all relations.")
		}
		
		contribution_factor = extend_contribution_factor(dag@relations_DAG, contribution_factor)
		if(length(setdiff(relation_levels, names(contribution_factor)))) {
			stop("Contribution factor should be provided for all relations.")
		}
		dag@term_env$IC_Wang_2007 = cpp_ic_wang(dag, unname(contribution_factor[relation_levels]))
	}
	dag@term_env$IC_Wang_2007
}
ADD_IC_METHOD("IC_Wang_2007", c("contribution_factor", "use_cache"))


