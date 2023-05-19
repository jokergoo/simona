


IC_annotation = function(dag, annotation, merge_to_ancestor = FALSE) {

	annotation2 = rep(list(character(0)), dag@n_terms)
	term2ind = structure(seq_along(dag@terms), names = dag@terms)

	cn = intersect(dag@terms, names(annotation))
	if(length(cn) == 0) {
		stop("No overlap between annotation names and terms.")
	}
	annotation2[ term2ind[cn] ] = annotation[cn]

	if(merge_to_ancestor) {
		annotation2 = lapply(dag@lt_offspring, function(ind) {
			unique(unlist(annotation2[ind]))
		})
	} else {
		lt_parents = dag@lt_parents
		ind = sample(length(lt_parents), min(length(lt_parents), 10))
		sapply(ind, function(i) {
			if(length(lt_parents[[i]]) > 0) {
				if(length(setdiff(annotation2[[ sample(lt_parents[[i]], 1) ]], annotation2[[i]])) > 0) {
					stop("Found a child term has more annotated items than its parent. Please set `merge_to_ancestor = TRUE`.")
				}
			}
		})
		annotation2 = lapply(annotation2, unique)
	}

	n = sapply(annotation2, length)
	p = n/max(n)
	ic = ifelse(n == 0, 0, -log(p))
	attr(ic, "N") = max(n)

	dag@stats_env$IC_annotation = ic
	
	dag@stats_env$IC_annotation
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_annotation")


###########################################
### max distance to root
dag_depth = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$dag_depth) || !use_cache) {
		lt_parents = dag@lt_parents
		n_terms = dag@n_terms

		d = cpp_dag_depth(lt_parents, n_terms)
		dag@stats_env$dag_depth = d
	}
	dag@stats_env$dag_depth
}

dag_depth_R = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$dag_depth) || !use_cache) {
		lt_parents = dag@lt_parents

		n_terms = dag@n_terms

		e = new.env()
		e$d = rep(NA, n_terms)
		calc_d = function(x) {
			if(x %in% dag@root) {
				e$d[x] = 0
				return(0)
			}
			if(!is.na(e$d[x])) {
				return(e$d[x])
			}
			v = max(sapply(lt_parents[[x]], calc_d)) + 1
			e$d[x] = v
			return(v)
		}

		for(i in seq_len(n_terms)) {
			calc_d(i)
		}

		dag@stats_env$dag_depth = e$d
	}
	dag@stats_env$dag_depth
}

###############################
### max distance to leaves
dag_height = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$dag_height) | !use_cache) {
		lt_children = dag@lt_children
		n_terms = dag@n_terms

		d = cpp_dag_height(lt_children, n_terms)
		dag@stats_env$dag_height = d
	}
	dag@stats_env$dag_height
}

dag_height_R = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$dag_height) || !use_cache) {
		lt_children = dag@lt_children

		n_terms = dag@n_terms

		e = new.env()
		e$d = rep(NA, n_terms)
		calc_d = function(x) {
			if(length(lt_children[[x]]) == 0) {
				e$d[x] = 0
				return(0)
			}
			if(!is.na(e$d[x])) {
				return(e$d[x])
			}
			v = max(sapply(lt_children[[x]], calc_d)) + 1
			e$d[x] = v
			return(v)
		}

		for(i in seq_len(n_terms)) {
			calc_d(i)
		}

		dag@stats_env$dag_height = e$d
	}
	dag@stats_env$dag_height
}

##########################
### IC universe
IC_universal_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@stats_env$IC_universal) || !use_cache) {

		lt_parents = dag@lt_parents
		n_children = dag@stats_env$n_children

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

		dag@stats_env$IC_universal = e$d
	}
	dag@stats_env$IC_universal
}

IC_universal_bfs = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_universal) || !use_cache) {
		lt_parents = dag@lt_parents

		n_children = dag@stats_env$n_children

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
		dag@stats_env$IC_universal = d
	}

	dag@stats_env$IC_universal
}

IC_universal = IC_universal_bfs
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_universal")


##############################################
### reachability is the number of ways for a node to reach the leaves
reachability_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@stats_env$reachability) || !use_cache) {
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

		dag@stats_env$reachability = e$reach
	} 
	dag@stats_env$reachability
}

reachability_bfs = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$reachability) || !use_cache) {
		lt_parents = dag@lt_parents
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
		dag@stats_env$reachability = reach
	}
	dag@stats_env$reachability
}

reachability = reachability_bfs  # total number of ways/paths to reach leaves

########################################
### Zhang et al
IC_Zhang_2006 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Zhang_2006) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@stats_env$IC_Zhang_2006 = -log(re/max(re))
	}
	dag@stats_env$IC_Zhang_2006
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Zhang_2006")


########################################
### Seco et al
IC_Seco_2004 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Seco_2004) || !use_cache) {
		re = reachability(dag, use_cache)
		dag@stats_env$IC_Seco_2004 = 1 - log(re + 1)/log(max(re) + 1)
	}
	dag@stats_env$IC_Seco_2004
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Seco_2004")


########################################
### Zhou et al
IC_Zhou_2008 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Zhou_2008) || !use_cache) {
		depth = dag_depth(dag, use_cache)
		ic_seco = IC_Seco_2004(dag, use_cache)
		
		sigma = 0.5
		dag@stats_env$IC_Zhou_2008 = sigma*ic_seco + (1-sigma)*log(ifelse(depth == 0, 1, depth))/log(max(depth))
	}
	dag@stats_env$IC_Zhou_2008
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Zhou_2008")


########################################
### Seddiqui et al
IC_Seddiqui_2010 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Seddiqui_2010) || !use_cache) {
		n_relations = dag@stats_env$n_parents + dag@stats_env$n_children
		n_edges = dag@n_relations
		n_nodes = dag@n_terms

		ic_seco = IC_Seco_2004(dag, use_cache)

		sigma = log(n_edges + 1)/( log(n_edges) + log(n_nodes) )
		dag@stats_env$IC_Seddiqui_2010 = (1 - sigma)*ic_seco + sigma*log(n_relations + 1)/log(n_edges + 1)
	}
	dag@stats_env$IC_Seddiqui_2010
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Seddiqui_2010")



#########################################
### Sanchez et al: information pass to leaves
IC_Sanchez_2011 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Sanchez_2011) || !use_cache) {
		n_leaves = length(dag@leaves)
		n_connected_leaves = sapply(dag@lt_offspring, function(x) length(intersect(x, dag@leaves)))

		n_ancestor = dag@stats_env$n_ancestor
		n_ancestor[n_ancestor == 0] = 1

		dag@stats_env$IC_Sanchez_2011 = -log( (n_connected_leaves/n_ancestor + 1)/(n_leaves + 1) )
	}
	dag@stats_env$IC_Sanchez_2011
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Sanchez_2011")


#########################################
### Meng et al
IC_Meng_2012 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Meng_2012) || !use_cache) {
		depth = dag_depth(dag, use_cache)
		max_depth = max(depth)
		lt_offspring = dag@lt_offspring
		n_terms = dag@n_terms
		dag@stats_env$IC_Meng_2012 = log(depth)/log(max_depth)*
		       (1 - log(sapply(lt_offspring, function(x) ifelse(length(x) == 1, 1, sum(1/depth[x[-1]]+1))))/log(n_terms))
	}
	dag@stats_env$IC_Meng_2012
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Meng_2012")



###############################
### totipotency
totipotency_recursive = function(dag, use_cache = TRUE) {

	if(is.null(dag@stats_env$totipotency) || !use_cache) {
		lt_parents = dag@lt_parents

		n_offspring = dag@stats_env$n_offspring

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

		dag@stats_env$totipotency = e$t
	}
	dag@stats_env$totipotency
}

totipotency_bfs = function(dag, use_cache = TRUE) {

	if(is.null(dag@stats_env$totipotency) || !use_cache) {
		lt_parents = dag@lt_parents

		n_offspring = dag@stats_env$n_offspring

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

		dag@stats_env$totipotency = t
	}
	dag@stats_env$totipotency
}

totipotency = totipotency_bfs


IC_Wang_2007 = function(dag, use_cache = TRUE) {
	if(is.null(dag@stats_env$IC_Wang_2007) || !use_cache) {
		g = dag@graph

		if(is.null(E(g)$relation)) {
			stop("'relations' was not set in `createOntologyDAG()`.")
		}
		w = -log(ifelse(E(g)$relation == "isa", 0.8, ifelse(E(g)$relation == "part of", 0.6, 0.7)))

		d = distances(g, weights = w, mode = "out")
		sv = 1/exp(d)
		ic = colSums(sv)

		dag@stats_env$IC_Wang_2007 = ic
	}
	dag@stats_env$IC_Wang_2007
}
ALL_IC_METHODS = c(ALL_IC_METHODS, "IC_Wang_2007")


OntologyIC = function(dag, method, use_cache = TRUE, annotation, merge_to_ancestor = FALSE) {
	IC_fun = get_IC_method(method)
	if(method == "IC_annotation") {
		ic = IC_fun(dag, annotation = annotation, merge_to_ancestor = merge_to_ancestor)
	} else {
		ic = IC_fun(dag, use_cache = use_cache)
	}
	ic
}

get_IC_method = function(method) {
	if(!method %in% ALL_IC_METHODS) {
		stop("Supported IC method are in `ALL_IC_METHODS`")
	}

	get(method, envir = topenv(), inherits = FALSE)
}
