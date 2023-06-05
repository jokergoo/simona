

#' Depth and height of the DAG
#' 
#' @param dag A `ontology_DAG` object.
#' @param use_cache Internally used.
#' @details
#' The depth of a term in the DAG is defined as the maximal distance to the root. The height
#' of a term in the DAG is the maximal finite distance to the leaf nodes.
#' 
#' @return An integer vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' dag_depth(dag)
#' dag_height(dag)
dag_depth = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_depth) || !use_cache) {
		d = cpp_dag_depth_bfs(dag)
		dag@term_env$dag_depth = d
	}
	dag@term_env$dag_depth
}

dag_depth_R = function(dag) {
	lt_parents = dag@lt_parents

	n_terms = dag@n_terms

	e = new.env()
	e$d = rep(NA_integer_, n_terms)
	calc_d = function(x) {
		if(x %in% dag@root) {
			e$d[x] = 0L
			return(0L)
		}
		if(!is.na(e$d[x])) {
			return(e$d[x])
		}
		v = max(sapply(lt_parents[[x]], calc_d)) + 1L
		e$d[x] = v
		return(v)
	}

	for(i in seq_len(n_terms)) {
		calc_d(i)
	}

	e$d 
}

dag_depth_R_bfs = function(dag) {
	lt_children = dag@lt_children

	n_terms = dag@n_terms

	d = rep(0L, n_terms)
	current_nodes = dag@root
	while(length(current_nodes)) {
		for(cr in current_nodes) {
			children = lt_children[[cr]]
			if(length(children)) {
				d[children] = base::pmax(d[children], d[cr] + 1L)
			}
		}
		current_nodes = unique(unlist(lt_children[current_nodes]))
	}
	d
}

#' @rdname dag_depth
#' @export
dag_height = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$dag_height) | !use_cache) {
		d = cpp_dag_height_bfs(dag)
		dag@term_env$dag_height = d
	}
	dag@term_env$dag_height
}

dag_height_R = function(dag) {
	lt_children = dag@lt_children

	n_terms = dag@n_terms

	e = new.env()
	e$d = rep(NA_integer_, n_terms)
	calc_d = function(x) {
		if(length(lt_children[[x]]) == 0) {
			e$d[x] = 0L
			return(0L)
		}
		if(!is.na(e$d[x])) {
			return(e$d[x])
		}
		v = max(sapply(lt_children[[x]], calc_d)) + 1L
		e$d[x] = v
		return(v)
	}

	for(i in seq_len(n_terms)) {
		calc_d(i)
	}

	e$d
}

dag_height_R_bfs = function(dag) {
	lt_parents = dag@lt_parents

	n_terms = dag@n_terms

	d = rep(0L, n_terms)
	current_nodes = dag@leaves
	while(length(current_nodes)) {
		for(cr in current_nodes) {
			parents = lt_parents[[cr]]
			if(length(parents)) {
				d[parents] = base::pmax(d[parents], d[cr] + 1L)
			}
		}
		current_nodes = unique(unlist(lt_parents[current_nodes]))
	}
	d
}

#' Number of annotated items
#' 
#' @param dag A `ontology_DAG` object.
#' @param use_cache Internally used.
#' 
#' @details
#' `annotation` argument should be set in [create_ontology_DAG()].
#' @returns An integer vector.
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db")
#' head(n_annotations(dag))
#' }
n_annotations = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$n_annotations)) {
		use_cache = FALSE
	} 
	if(!use_cache) {
		if(length(dag@annotation$list) == 0) {
			stop("You should specify `annotation` in `create_ontology_DAG()`.")
		}
		n_all_anno = length(dag@annotation$names)
	
		n = cpp_n_annotations(dag)
		
		attr(n, "N") = n_all_anno

		dag@term_env$n_annotations = n
	}
	
	dag@term_env$n_annotations
}

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
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_annotation")



##########################
### IC universe
IC_universal_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@term_env$IC_universal) || !use_cache) {

		lt_parents = dag@lt_parents
		n_children = dag@term_env$n_children

		n_terms = dag@n_terms

		e = new.env()
		e$d = rep(NA, n_terms)
		calc_ic = function(x) {
			if(x %in% dag@root) {	
				e$d[x] = 0
				return(0)
			}
			if(!is.na(e$d[x])) {
				return(e$d[x])
			}
			v = sum(sapply(lt_parents[[x]], function(t) {
				calc_ic(t) + log(n_children[t])
			}))
			e$d[x] = v
			return(v)
		}

		for(i in seq_len(n_terms)) {
			calc_ic(i)
		}

		dag@term_env$IC_universal = e$d
	}
	dag@term_env$IC_universal
}

IC_universal_bfs = function(dag, use_cache = TRUE) {
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

IC_universal = IC_universal_bfs
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_universal")


##############################################
### reachability is the number of ways for a node to reach the leaves
reachability_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@term_env$reachability) || !use_cache) {
		lt_parents = dag@lt_parents
		lt_children = dag@lt_children

		terms = dag@terms

		e = new.env()
		e$reach = rep(NA, length(terms))
		calc_re = function(x) {
			if(length(lt_children[[x]]) == 0) {
				e$reach[x] = 1
				return(1)
			}
			if(!is.na(e$reach[x])) {
				return(e$reach[x])
			}
			v = sum(sapply(lt_children[[x]], calc_re))
			e$reach[x] = v
			return(v)
		}

		for(i in seq_along(terms)) {
			calc_re(i)
		}

		dag@term_env$reachability = e$reach
	} 
	dag@term_env$reachability
}

reachability_bfs = function(dag, use_cache = TRUE) {
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

reachability = reachability_bfs  # total number of ways/paths to reach leaves

########################################
### Zhang et al
IC_Zhang_2006 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Zhang_2006) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@term_env$IC_Zhang_2006 = -log(re/max(re))
	}
	dag@term_env$IC_Zhang_2006
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Zhang_2006")


########################################
### Seco et al
IC_Seco_2004 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Seco_2004) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@term_env$IC_Seco_2004 = 1 - log(re + 1)/log(max(re) + 1)
	}
	dag@term_env$IC_Seco_2004
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Seco_2004")


########################################
### Zhou et al
IC_Zhou_2008 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Zhou_2008) || !use_cache) {
		depth = dag_depth(dag, use_cache)
		ic_seco = IC_Seco_2004(dag, use_cache)
		
		sigma = 0.5
		dag@term_env$IC_Zhou_2008 = sigma*ic_seco + (1-sigma)*log(ifelse(depth == 0, 1, depth))/log(max(depth))
	}
	dag@term_env$IC_Zhou_2008
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Zhou_2008")


########################################
### Seddiqui et al
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
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Seddiqui_2010")



#########################################
### Sanchez et al: information pass to leaves
IC_Sanchez_2011 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Sanchez_2011) || !use_cache) {
		n_leaves = length(dag@leaves)
		n_connected_leaves = n_leaves(dag)

		n_ancestor = dag@term_env$n_ancestor
		n_ancestor[n_ancestor == 0] = 1

		dag@term_env$IC_Sanchez_2011 = -log( (n_connected_leaves/n_ancestor + 1)/(n_leaves + 1) )
	}
	dag@term_env$IC_Sanchez_2011
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Sanchez_2011")



#########################################
### Meng et al
IC_Meng_2012 = function(dag, use_cache = TRUE, correct = FALSE) {
	if(is.null(dag@term_env$IC_Meng_2012) || !use_cache) {
		dag@term_env$IC_Meng_2012 = cpp_ic_meng(dag, correct + 0)
	}
	dag@term_env$IC_Meng_2012
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Meng_2012")



###############################
### totipotency
totipotency_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@term_env$totipotency) || !use_cache) {
		lt_parents = dag@lt_parents

		n_offspring = n_offspring(dag) + 1
		n_terms = dag@n_terms

		e = new.env()
		e$t = rep(NA, n_terms)
		calc_t = function(x) {
			if(x %in% dag@root) {	
				e$t[x] = 1
				return(1)
			}
			if(!is.na(e$t[x])) {
				return(e$t[x])
			}

			v = sapply(lt_parents[[x]], FUN = function(t, x) {
					n_offspring[x]/n_offspring[t] * calc_t(t)
				}, x)
			e$t[x] = mean(v)
			return(e$t[x])
		}

		for(i in seq_len(n_terms)) {
			calc_t(i)
		}

		dag@term_env$totipotency = e$t
	}
	dag@term_env$totipotency
}

totipotency_bfs = function(dag, use_cache = TRUE) {

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

totipotency = totipotency_bfs



# IC_Wang_2007 = function(dag, use_cache = TRUE, contribution_factor = c("isa" = 0.8, "part of" = 0.6)) {
# 	if(is.null(dag@term_env$IC_Wang_2007) || !use_cache) {
# 		if(length(dag@lt_children_relations) == 0) {
# 			stop("`relations` is not set when creating the ontology_DAG object.")
# 		}
# 		relation_levels = attr(dag@lt_children_relations, "levels")
# 		if(length(setdiff(relation_levels, names(contribution_factor)))) {
# 			stop("Contribution factor should be provided for all relations.")
# 		}
# 		dag@term_env$IC_Wang_2007 = cpp_ic_wang(dag, unname(contribution_factor[relation_levels]))
# 	}
# 	dag@term_env$IC_Wang_2007
# }
# ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Wang_2007")

