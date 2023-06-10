
# remove terms with zero annotation and the root term
validate_annotated_terms = function(dag, id) {
	n_anno = n_annotations(dag)[id]

	l1 = n_anno == 0
	if(any(l1)) {
		message("removed", sum(l1), "terms with no annotation.")
		
		if(length(sum(!l1)) == 0) {
			stop("No term is lef.")
		}
	}
	l2 = id %in% dag@root
	if(any(l2)) {
		message("removed root term.")
		
		if(length(sum(!l1 & !l2)) == 0) {
			stop("No term is lef.")
		}
	}

	!l1 & !l2
}

Sim_Ling_1998 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)
	browser()
	sim = 2*ic_mica/outer(ic[id], ic[id], "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Ling_1998")


Sim_Resnik_1999 = function(dag, terms, norm_method = "Nmax") {
	IC_method = "IC_annotation"

	if(norm_method %in% c("Nunif", "Nmax", "Nunivers", "none")) {
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
		sim = ic_mica/(outer(ic, ic, max))
	} else  if(norm_method == "none") {
		sim = ic_mica
	}
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Resnik_1999")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_FaITH_2010")


Sim_PS_2008 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = get(IC_method)(dag)[id]
	
	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = 3*ic_mica - outer(ic, ic, "+")
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	sim[sim < 0] = 0
	diag(sim) = 1

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_PS_2008")


Sim_Relevance_2006 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = get_Sim_method("Sim_Ling_1998")(dag, id)

	eps = 1 - exp(-ic_mica)
	eps*sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Relevance_2006")


Sim_SimIC_2010 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = get_Sim_method("Sim_Ling_1998")(dag, id)

	eps = 1 - 1/(1 + ic_mica)
	eps*sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_SimIC_2010")


Sim_EISI_2015 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)
	ic = term_IC(dag, IC_method)

	sim = get_Sim_method("Sim_Ling_1998")(dag, id)

	eps = cpp_eps_EISI(dag, id, ic)

	eps*sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_EISI_2015")


Sim_XGraSM_2013 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method)
	ic = term_IC(dag, IC_method)

	sim = get_Sim_method("Sim_Ling_1998")(dag, id)
	eps = cpp_common_ancestor_mean_IC_XGraSM(dag, id, ic)/ic_mica

	eps*sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_XGraSM_2013")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_AIC_2014")


Sim_Zhang_2006 = function(dag, terms) {
	IC_method = "IC_Zhang_2006"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)[id]

	ic_mica = MICA_IC(dag, id, IC_method)
	
	sim = 2*ic_mica/outer(ic, ic, "+")
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Zhang_2006")


Sim_universal = function(dag, terms) {
	IC_method = "IC_universal"

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)[id]

	ic_mica = MICA_IC(dag, id, IC_method)

	sim = 2*ic_mica/outer(ic, ic, max)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_universal")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Wang_2007")

###########################################
#### edge-based
###########################################

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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Rada_1989")


Sim_Resnik_edge_2012 = function(dag, terms, distance = "longest_distances_via_LCA") {

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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Resnik_edge_2012")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Leocock_1998")


Sim_Wu_1994 = function(dag, terms) {

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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Wu_1994")


Sim_Slimani_2006 = function(dag, terms, distance = "longest_distances_via_LCA") {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)

	sim_wp = 2*lca_depth/(cpp_longest_distances_via_LCA(dag, id) + 2*lca_depth)

	lamda = cpp_is_reachable(dag, id)

	cf = (1 - lambda)*(outer(depth[id], depth[id], min) - lca_depth) +
	     lambda/(abs(outer(depth[id], depth[id], "-")) + 1)
	diag(cf) = 1

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Slimani_2006")


Sim_Shenoy_2012 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)

	sim_wp = 2*lca_depth/(cpp_longest_distances_via_LCA(dag, id) + 2*lca_depth)

	lamda = cpp_is_reachable(dag, id)

	max_depth = max(dag_depth(dag))
	dist = cpp_longest_distances_via_LCA(dag, id)

	cf = exp(- lambda*dist/max_depth)
	diag(cf) = 1

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Shenoy_2012")


Sim_Pekar_2002 = function(dag, terms) {

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id)

	sim = lca_depth/(cpp_longest_distances_via_LCA(dag, id) + lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Pekar_2002")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Stojanovic_2001")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Wang_edge_2012")


Sim_Zhong_2002 = function(dag, terms, depth_via_LCA = TRUE) {
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	sim = cpp_sim_zhong(dag, id, depth_via_LCA)

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Zhong_2002")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_AlMubaid_2006")


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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Li_2003")


###########################################
#### hybrid
###########################################

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

	sim = max_depth/(max_depth + dist) * lca_depth/(lca_depth + 0.5*outer(height, height, "+"))
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_RSS_2013")


Sim_HRSS_2013 = function(dag, terms, IC_method) {
	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)[id]
	ic_mica = MICA_IC(dag, id, IC_method)

	MIL_term = cpp_max_leaves_id(dag, id, ic)

	alpha = abs(ic[dag@root] - ic_mica)
	ic_dist_leaf = ic[id] - ic[MIL_term]
	beta = outer(ic_dist_leaf, ic_dist_leaf, "+")/2

	sim = 1/(1 + abs(outer(ic[id], ic[id], "-")))*alpha/(alpha + beta)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_HRSS_2013")


Sim_Shen_2010 = function(dag, terms, IC_method) {
	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method)

	sim = cpp_sim_shen(dag, id, ic);
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Shen_2010")


Sim_SSDD_2013 = function(dag, terms) {
	id = term_to_node_id(dag, terms, strict = FALSE)
	t = totipotency(dag)

	sim = cpp_sim_SSDD(dag, id, t)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	sim
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_SSDD_2013")


Sim_Jiang_1997 = function(dag, terms, IC_method, normalization = "max") {

	id = term_to_node_id(dag, terms)
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
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Jiang_1997")

######################
## count-based
######################

Sim_Kappa = function(dag, terms, universe = NULL) {
	.sim_overlap(dag, terms, universe, method = "kappa")
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Kappa")


Sim_Jaccard = function(dag, terms, universe = NULL) {
	.sim_overlap(dag, terms, universe, method = "jaccard")
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Jaccard")

Sim_Dice = function(dag, terms, universe = NULL) {
	.sim_overlap(dag, terms, universe, method = "dice")
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Dice")

Sim_Overlap = function(dag, terms, universe = NULL) {
	.sim_overlap(dag, terms, universe, method = "overlap")
}
ALL_TERM_SIM_METHODS = c(ALL_TERM_SIM_METHODS, "Sim_Overlap")

.sim_overlap = function(dag, terms, universe = NULL, method = c("kappa", "jaccard", "dice", "overlap")) {

	check_pkg("proxyC", bioc = FALSE)

	id = term_to_node_id(dag, terms, strict = FALSE)
	if(!is.null(universe)) {
		universe = term_to_node_id(dag, universe, strict = FALSE)
		id = intersect(id, universe)
	} else {
		universe = seq_along(dag@annotation$names)
	}

	if(length(dag@annotation$list) == 0) {
		stop("`annotation` should be set in `create_ontology_DAG()`.")
	}

	mg = cpp_get_term_annotations(dag, id)
	mg = mg[, universe, drop = FALSE]

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


