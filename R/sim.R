

#' Sim_Ling_1998
#' 
#' @section method:
#' ## Sim_Ling_1998
#' 
#' 
#' Paper link: <https://dl.acm.org/doi/10.5555/645527.657297>.
#' 
#' @rdname temp__Sim_Ling_1998
Sim_Ling_1998 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = 2*ic_mica/outer(ic, ic, "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_Ling_1998")


#' Sim_Resnik_1999
#' 
#' @section method:
#' ## Sim_Resnik_1999
#' 
#' 
#' Paper link: <https://doi.org/10.1613/jair.514>, <https://doi.org/10.1186/1471-2105-9-S5-S4>, <https://doi.org/10.1186/1471-2105-11-562>, <https://doi.org/10.1155/2013/292063>.
#' 
#' @rdname temp__Sim_Resnik_1999
Sim_Resnik_1999 = function(dag, terms, norm_method = "Nmax") {
	IC_method = "IC_annotation"

	if(!norm_method %in% c("Nunif", "Nmax", "Nunivers", "none")) {
		stop("`norm_method` can only be in 'Nunif', 'Nmax', 'Nunivers, 'none'.")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]
	
	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	if(norm_method == "Nunif") {
		max_n = attr(n_annotations(dag), "N")
		sim = ic_mica/log(max_n)
	} else if(norm_method == "Nmax") {
		sim = ic_mica/max(ic)
	} else if(norm_method == "Nunivers") {
		sim = ic_mica/(outer(ic, ic, pmax))
	} else  if(norm_method == "none") {
		sim = ic_mica
	}
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ADD_TERM_SIM_METHOD("Sim_Resnik_1999")

#' Sim_FaITH_2010
#' 
#' @section method:
#' ## Sim_FaITH_2010
#' 
#' <https://doi.org/10.1007/978-3-642-17746-0_39>.
#' @rdname temp__Sim_FaITH_2010
Sim_FaITH_2010 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = ic_mica/(outer(ic, ic, "+") - ic_mica)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ADD_TERM_SIM_METHOD("Sim_FaITH_2010")

#' Sim_PS_2008
#' 
#' @section method:
#' ## Sim_PS_2008
#' 
#' <https://doi.org/10.1016/j.datak.2009.06.008>.
#' 
#' @rdname temp__Sim_PS_2008
Sim_PS_2008 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]
	
	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = 3*ic_mica - outer(ic, ic, "+")
	# scale into 0~1
	max_ic = max(get(IC_method)(dag))
	sim = (sim + 2*max_ic)/3/max_ic
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_PS_2008")

#' Sim_Relevance_2006
#' 
#' @section method:
#' ## Sim_Relevance_2006
#' 
#' <https://doi.org/10.1186/1471-2105-7-302>
#' 
#' @rdname temp__Sim_Relevance_2006
Sim_Relevance_2006 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = get_term_sim_method("Sim_Ling_1998")(dag, id)

	eps = 1 - exp(-ic_mica)
	eps*sim
}
ADD_TERM_SIM_METHOD("Sim_Relevance_2006")

#' Sim_SimIC_2010
#' 
#' @section method:
#' ## Sim_SimIC_2010
#' 
#' <https://doi.org/10.48550/arXiv.1001.0958>.
#' 
#' @rdname temp__Sim_SimIC_2010
Sim_SimIC_2010 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = get_term_sim_method("Sim_Ling_1998")(dag, id)

	eps = 1 - 1/(1 + ic_mica)
	eps*sim
}
ADD_TERM_SIM_METHOD("Sim_SimIC_2010")

#' Sim_EISI_2015
#' 
#' @section method:
#' ## Sim_EISI_2015
#' 
#' <https://doi.org/10.1016/j.gene.2014.12.062>
#' 
#' @rdname temp__Sim_EISI_2015
Sim_EISI_2015 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic = term_IC(dag, IC_method)

	mean_ic = cpp_common_ancestor_mean_IC_EISI(dag, id, ic)
	sim = mean_ic/outer(ic[id], ic[id], "+")*2
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_EISI_2015")

#' Sim_XGraSM_2013
#' 
#' @section method:
#' ## Sim_XGraSM_2013
#' 
#' <https://doi.org/10.1186/1471-2105-14-284>
#' 
#' @rdname temp__Sim_XGraSM_2013
Sim_XGraSM_2013 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic = term_IC(dag, IC_method)

	mean_ic = cpp_common_ancestor_mean_IC_XGraSM(dag, id, ic)
	sim = mean_ic/outer(ic[id], ic[id], "+")*2
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_XGraSM_2013")

#' Sim_GraSM_2005
#' 
#' @section method:
#' ## Sim_GraSM_2005
#' 
#' <https://doi.org/10.1145/1099554.1099658>.
#' 
#' @rdname temp__Sim_GraSM_2005
Sim_GraSM_2005 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic = term_IC(dag, IC_method)

	mean_ic = cpp_common_ancestor_mean_IC_GraSM(dag, id, ic)
	sim = mean_ic/outer(ic[id], ic[id], "+")*2
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_GraSM_2005")

#' Sim_AIC_2014
#' 
#' @section method:
#' ## Sim_AIC_2014
#' 
#' <https://doi.org/10.1109/tcbb.2013.176>.
#' 
#' @rdname temp__Sim_AIC_2014
Sim_AIC_2014 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic = term_IC(dag, IC_method)
	sim = cpp_sim_aic(dag, id, ic)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_AIC_2014")

#' Sim_Zhang_2006
#' 
#' @section method:
#' ## Sim_Zhang_2006
#' 
#' 
#' @rdname temp__Sim_Zhang_2006
Sim_Zhang_2006 = function(dag, terms) {
	IC_method = "IC_Zhang_2006"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)[id]

	ic_mica = MICA_IC(dag, id, IC_method)
	
	sim = 2*ic_mica/outer(ic, ic, "+")
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_Zhang_2006")

#' Sim_universal
#' 
#' @section method:
#' what is Sim_universal
#' @rdname temp__Sim_universal
Sim_universal = function(dag, terms) {
	IC_method = "IC_universal"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)[id]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = 2*ic_mica/outer(ic, ic, pmax)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_universal")

#' Sim_Wang_2007
#' 
#' @section method:
#' ## Sim_Wang_2007
#' 
#' Paper link: <https://doi.org/10.1093/bioinformatics/btm087>.
#' 
#' @rdname temp__Sim_Wang_2007
Sim_Wang_2007 = function(dag, terms, contribution_factor = c("isa" = 0.8, "part of" = 0.6)) {
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

	id = term_to_node_id(dag, terms, strict = FALSE)

	sim = cpp_sim_wang(dag, id, unname(contribution_factor[relation_levels]))
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Wang_2007")

###########################################
#### edge-based
###########################################

#' Sim_Rada_1989
#' 
#' @section method:
#' ## Sim_Rada_1989
#' 
#' 
#' <https://doi.org/10.1109/21.24528>
#' 
#' @rdname temp__Sim_Rada_1989
Sim_Rada_1989 = function(dag, terms, distance = "longest_distances_via_LCA") {

	id = term_to_node_id(dag, terms, strict = FALSE)

	if(distance == "shortest_distances_via_CA") {
		dist = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dist = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	sim = 1/(1 + dist)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Rada_1989")

#' Sim_Resnik_edge_2005
#' 
#' @section method:
#' ## Sim_Resnik_edge_2005
#' 
#' <https://doi.org/10.1145/1097047.1097051>
#' 
#' @rdname temp__Sim_Resnik_edge_2005
Sim_Resnik_edge_2005 = function(dag, terms, distance = "longest_distances_via_LCA") {

	id = term_to_node_id(dag, terms)

	max_depth = max(dag_depth(dag))
	if(distance == "shortest_distances_via_CA") {
		dist = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dist = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	sim = 1 - dist/2/max_depth
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Resnik_edge_2005")

#' Sim_Leocock_1998
#' 
#' @section method:
#' ## Sim_Leocock_1998
#' 
#' <https://ieeexplore.ieee.org/document/6287675>
#' 
#' @rdname temp__Sim_Leocock_1998
Sim_Leocock_1998 = function(dag, terms, distance = "longest_distances_via_LCA") {
	
	id = term_to_node_id(dag, terms, strict = FALSE)

	max_depth = max(dag_depth(dag))
	if(distance == "shortest_distances_via_CA") {
		dist = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dist = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	sim = 1 - log(dist + 1)/log(2*max_depth + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Leocock_1998")

#' Sim_WP_1994
#' 
#' @section method:
#' what is Sim_WP_1994
#' @rdname temp__Sim_WP_1994
Sim_WP_1994 = function(dag, terms) {

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	lca_depth = LCA_depth(dag, id)
	sim = 2*lca_depth/(cpp_longest_distances_via_LCA(dag, id) + 2*lca_depth)
	
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	l = which(id %in% dag@root)
	if(any(l)) {
		sim[l, l] = 1
	}

	sim
}
ADD_TERM_SIM_METHOD("Sim_WP_1994")

#' Sim_Slimani_2006
#' 
#' @section method:
#' ## Sim_Slimani_2006
#' 
#' <https://zenodo.org/record/1075130>
#' 
#' @rdname temp__Sim_Slimani_2006
Sim_Slimani_2006 = function(dag, terms, distance = "longest_distances_via_LCA") {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)

	sim_wp = 2*lca_depth/(cpp_longest_distances_via_LCA(dag, id) + 2*lca_depth)

	lambda = cpp_is_reachable(dag, id)
	ltd = cpp_longest_distances_from_LCA(dag, id)
	dist_left = ltd$left
	dist_right = ltd$right

	cf = (1 - lambda)*(pmin(dist_left, dist_right) + 1) + 
	     lambda/(dist_left + dist_right)  ## the two terms are ancestor/offspring, and one of `dist_left` and dist_right is zero
	diag(cf) = 1

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Slimani_2006")

#' Sim_Shenoy_2012
#' 
#' @section method:
#' ## Sim_Shenoy_2012
#' 
#' <https://doi.org/10.48550/arXiv.1211.4709>.
#' 
#' @rdname temp__Sim_Shenoy_2012
Sim_Shenoy_2012 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)

	sim_wp = 2*lca_depth/(cpp_longest_distances_via_LCA(dag, id) + 2*lca_depth)

	lambda = cpp_is_reachable(dag, id)

	max_depth = max(dag_depth(dag))
	dist = cpp_longest_distances_via_LCA(dag, id)

	cf = exp(- lambda*dist/max_depth)
	diag(cf) = 1

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Shenoy_2012")

#' Sim_Pekar_2002
#' 
#' @section method:
#' ## Sim_Pekar_2002
#' 
#' <https://aclanthology.org/C02-1090/>.
#' 
#' @rdname temp__Sim_Pekar_2002
Sim_Pekar_2002 = function(dag, terms) {

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)

	sim = lca_depth/(cpp_longest_distances_via_LCA(dag, id) + lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Pekar_2002")

#' Sim_Stojanovic_2001
#' 
#' @section method:
#' ## Sim_Stojanovic_2001
#' 
#' <https://doi.org/10.1145/500737.500762>
#' 
#' @rdname temp__Sim_Stojanovic_2001
Sim_Stojanovic_2001 = function(dag, terms) {

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = which(id %in% dag@root)
	if(any(l)) {
		message("remove root term.")
		id = id[!l]
	}

	depth = dag_depth(dag)
	lca_term = max_ancestor_id(dag, id, depth, in_labels = FALSE)
	lca_depth = structure(depth[lca_term], dim = dim(lca_term))
	
	sim = lca_depth/(outer(depth[id], depth[id], "+") - lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Stojanovic_2001")

#' Sim_Wang_edge_2012
#' 
#' @section method:
#' ## Sim_Wang_edge_2012
#' 
#' <https://doi.org/10.1186/1477-5956-10-s1-s18>.
#' 
#' @rdname temp__Sim_Wang_edge_2012
Sim_Wang_edge_2012 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = id %in% dag@root
	if(any(l)) {
		message("remove root term.")
		id = id[!l]
	}

	sim = cpp_sim_wang_edge(dag, id)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Wang_edge_2012")

#' Sim_Zhong_2002
#' 
#' @section method:
#' ## Sim_Zhong_2002
#' 
#' <https://doi.org/10.1007/3-540-45483-7_8>
#' 
#' @rdname temp__Sim_Zhong_2002
Sim_Zhong_2002 = function(dag, terms, depth_via_LCA = TRUE) {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	sim = cpp_sim_zhong(dag, id, depth_via_LCA)

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Zhong_2002")

#' Sim_AlMubaid_2006
#' 
#' @section method:
#' ## Sim_AlMubaid_2006
#' 
#' <https://doi.org/10.1109/IEMBS.2006.259235>
#' 
#' @rdname temp__Sim_AlMubaid_2006
Sim_AlMubaid_2006 = function(dag, terms, distance = "longest_distances_via_LCA") {
		
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)
	max_depth = max(dag_depth(dag))

	if(distance == "shortest_distances_via_CA") {
		dsp = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	dist = log(1 + dsp*(max_depth - lca_depth))
	max_dist = log(1 + 2*max_depth*max_depth)
	sim = (max_dist - dist)/max_dist
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_AlMubaid_2006")

#' Sim_Li_2003
#' 
#' @section method:
#' ## Sim_Li_2003
#' 
#' <https://doi.org/10.1109/TKDE.2003.1209005>
#' 
#' @rdname temp__Sim_Li_2003
Sim_Li_2003 = function(dag, terms, distance = "longest_distances_via_LCA") {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)

	if(distance == "shortest_distances_via_CA") {
		dsp = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	alpha = 0.2
	beta = 0.8

	sim = exp(-alpha*dsp)*atan(beta*lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Li_2003")


###########################################
#### hybrid
###########################################

#' Sim_RSS_2013
#' 
#' @section method:
#' ## Sim_RSS_2013
#' 
#' <https://doi.org/10.1371/journal.pone.0066745>
#' 
#' @rdname temp__Sim_RSS_2013
Sim_RSS_2013 = function(dag, terms, distance = "longest_distances_via_LCA") {
	id = term_to_node_id(dag, terms, strict = FALSE)

	lca_depth = LCA_depth(dag, id)
	max_depth = max(dag_depth(dag))
	height = dag_height(dag)[id]

	if(distance == "shortest_distances_via_CA") {
		dsp = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}

	sim = max_depth/(max_depth + dsp) * lca_depth/(lca_depth + 0.5*outer(height, height, "+") + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_RSS_2013")

#' Sim_HRSS_2013
#' 
#' @section method:
#' ## Sim_HRSS_2013
#' 
#' <https://doi.org/10.1371/journal.pone.0066745>
#' 
#' @rdname temp__Sim_HRSS_2013
Sim_HRSS_2013 = function(dag, terms, IC_method = "IC_annotation") {
	id = term_to_node_id(dag, terms, strict = FALSE)
	if(IC_method == "IC_annotation") {
		l = validate_annotated_terms(dag, id)
		id = id[l]
	}
	ic = term_IC(dag, IC_method)
	ic_mica = MICA_IC(dag, id, IC_method)

	MIL_term = cpp_max_leaves_id(dag, id, ic)

	alpha = abs(ic[dag@root] - ic_mica)
	ic_dist_leaf = abs(ic[id] - ic[MIL_term])
	beta = outer(ic_dist_leaf, ic_dist_leaf, "+")/2

	sim = 1/(1 + abs(outer(ic[id], ic[id], "-")))*alpha/(alpha + beta + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_HRSS_2013")

#' Sim_Shen_2010
#' 
#' @section method:
#' ## Sim_Shen_2010
#' 
#' 
#' <https://doi.org/10.1109/BIBM.2010.5706623>
#' 
#' @rdname temp__Sim_Shen_2010
Sim_Shen_2010 = function(dag, terms, IC_method = "IC_annotation") {
	id = term_to_node_id(dag, terms, strict = FALSE)
	if(IC_method == "IC_annotation") {
		l = validate_annotated_terms(dag, id)
		id = id[l]
	}
	ic = term_IC(dag, IC_method)

	sim = cpp_sim_shen(dag, id, ic);
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Shen_2010")

#' Sim_SSDD_2013
#' 
#' @section method:
#' ## Sim_SSDD_2013
#' 
#' <https://doi.org/10.1016/j.ygeno.2013.04.010>
#' 
#' @rdname temp__Sim_SSDD_2013
Sim_SSDD_2013 = function(dag, terms) {
	id = term_to_node_id(dag, terms, strict = FALSE)
	t = totipotency(dag)

	sim = cpp_sim_SSDD(dag, id, t)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	sim
}
ADD_TERM_SIM_METHOD("Sim_SSDD_2013")

#' Sim_Jiang_1997
#' 
#' @section method:
#' ## Sim_Jiang_1997
#' 
#' <https://aclanthology.org/O97-1002/>
#' 
#' @rdname temp__Sim_Jiang_1997
Sim_Jiang_1997 = function(dag, terms, normalization = "max") {

	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)
	ic = term_IC(dag, IC_method)[id]
	max_ic = max(ic)

	dist = outer(ic, ic, "+") - 2*ic_mica
	
	if(normalization == "max") {
		sim = 1 - dist/2/max_ic
	} else if(normalization == "Couto") {
		sim = dist/max_ic
		sim[sim > 1] = 1
		sim = 1 - sim
	} else if(normalization == "Lin") {
		sim = 2*ic_mica/outer(ic, ic, "+")
	} else if(normalization == "Garla") {
		sim = 1 - log(dist + 1)/log(2*max_ic + 1)
	} else if(normalization == "log-Lin") {
		sim = 1 - log(dist + 1)/log(outer(ic, ic, "+") + 1)
	} else if(normalization == "Rada") {
		sim = 1/(1 + dist)
	} else {
		stop("wrong normalization method.")
	}

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Jiang_1997")

######################
## count-based
######################

#' Sim_Kappa
#' 
#' @section method:
#' ## Sim_Kappa
#' @rdname temp__Sim_Kappa
Sim_Kappa = function(dag, terms, anno_universe = NULL) {
	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "kappa")
}
ADD_TERM_SIM_METHOD("Sim_Kappa")

#' Sim_Jaccard
#' 
#' @section method:
#' ## Sim_Jaccard
#' @rdname temp__Sim_Jaccard
Sim_Jaccard = function(dag, terms, anno_universe = NULL) {
	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "jaccard")
}
ADD_TERM_SIM_METHOD("Sim_Jaccard")

#' Sim_Dice
#' 
#' @section method:
#' ## Sim_Dice
#' @rdname temp__Sim_Dice
Sim_Dice = function(dag, terms, anno_universe = NULL) {
	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "dice")
}
ADD_TERM_SIM_METHOD("Sim_Dice")

#' Sim_Overlap
#' 
#' @section method:
#' ## Sim_Overlap
#' @rdname temp__Sim_Overlap
Sim_Overlap = function(dag, terms, anno_universe = NULL) {
	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "overlap")
}
ADD_TERM_SIM_METHOD("Sim_Overlap")

#' @importFrom methods as
.sim_overlap = function(dag, id, anno_universe = NULL, method = c("kappa", "jaccard", "dice", "overlap")) {

	check_pkg("proxyC", bioc = FALSE)

	if(!is.null(anno_universe)) {
		anno_universe = which(dag@annotation$names %in% anno_universe)
	} else {
		anno_universe = seq_along(dag@annotation$names)
	}

	if(length(dag@annotation$list) == 0) {
		stop("`annotation` should be set in `create_ontology_DAG()`.")
	}

	mg = cpp_get_term_annotations(dag, id)
	mg = mg[, anno_universe, drop = FALSE]

	mg = as(mg, "sparseMatrix")

	method = match.arg(method)[1]
	if(method == "kappa") {
		mat = kappa_dist(mg)
	} else if(method == "overlap") {
		mat = overlap_dist(mg)
	} else {
		mat = proxyC::simil(mg, method = method)
	}

	mat = as.matrix(mat)
	diag(mat) = 1
	rownames(mat) = colnames(mat) = dag@terms[id]
	return(mat)
}

kappa_dist = function(m) {
	tab = ncol(m)
	oab = proxyC::simil(m, method = "simple matching")
	m1 = rowSums(m)
	m2 = abs(rowSums(m - 1))
	aab = (outer(m1, m1) + outer(m2, m2))/tab/tab
	k = (oab - aab)/(1 - aab)
	return(k)
}

overlap_dist = function(m) {
	n = rowSums(m)
	proxyC::simil(m, method = "dice")*outer(n, n, "+")/2/outer(n, n, pmin)
}


