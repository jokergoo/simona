


#' IC_annotation
#' 
#' @section method:
#' ## what is IC_annotation
#' blablabla
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
#' ## what is IC_universal
#' blablabla
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
#' what is IC_Zhang_2006
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
#' what is IC_Seco_2004
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
#' what is IC_Zhou_2008
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
#' what is IC_Seddiqui_2010
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
#' what is IC_Sanchez_2011
#' @rdname temp__IC_Sanchez_2011
IC_Sanchez_2011 = function(dag, use_cache = TRUE) {
	if(is.null(dag@term_env$IC_Sanchez_2011) || !use_cache) {
		nl = length(dag@leaves)
		n_connected_leaves = n_leaves(dag)

		na = n_ancestors(dag)
		na[na == 0] = 1

		dag@term_env$IC_Sanchez_2011 = -log( (n_connected_leaves/na + 1)/(nl + 1) )
	}
	dag@term_env$IC_Sanchez_2011
}
ADD_IC_METHOD("IC_Sanchez_2011")



#########################################
### Meng et al

#' IC_Meng_2012
#' 
#' @section method:
#' what is IC_Meng_2012
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
#' what is IC_Wang_2007
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


