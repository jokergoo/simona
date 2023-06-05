


###########################################
#### node-based
###########################################

MICA = function(dag, term_id, IC_method, return_term_id = FALSE) {

	if(!IC_method %in% ALL_IC_METHODS) {
		stop("IC method: '", IC_method, "' is not supported.")
	}

	ic = calc_IC(dag, IC_method)

	if(any(ic < 0)) {
		stop("IC cannot be negative.")
	}

	if(return_term_id) {
		cm = cpp_max_ancestor_id(dag, term_id, ic)
	} else {
		cm = cpp_max_ancestor_v(dag, term_id, ic)
	}
	dimnames(cm) = list(dag@terms[term_id], dag@terms[term_id])

	cm
}

#' Most informative common ancestor (MICA) and lowest common ancestor (LCA)
#' 
#' @param dag A `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param IC_method An IC method. Valid values are in [ALL_IC_METHODS].
#' 
#' @return 
#' `MICA_term()` returns a character matrix of the MICA terms. 
#' `MICA_IC()` returns a numeric matrix of the IC of the MICA terms.
#' `LCA_term()` returns a character matrix of the LCA term.
#' `LCA_depth()` reutrns an integer matrix of the depth of the LCA terms.
MICA_term = function(dag, terms, IC_method) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	cm = MICA(dag, id, IC_method, return_term_id = TRUE)
	structure(dag@terms[cm], dim = dim(cm), dimnames = dimnames(cm))
}

#' @rdname MICA_term
MICA_IC = function(dag, terms, IC_method) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	MICA(dag, id, IC_method, return_term_id = FALSE)
}

LCA = function(dag, term_id, return_term_id = FALSE) {

	depth = dag_depth(dag)

	if(return_term_id) {
		cm = cpp_max_ancestor_id(dag, term_id, depth)
	} else {
		cm = cpp_max_ancestor_v(dag, term_id, depth)
	}
	dimnames(cm) = list(dag@terms[term_id], dag@terms[term_id])

	cm
}

#' @rdname MICA_term
LCA_term = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	cm = LCA(dag, id, return_term_id = TRUE)
	structure(dag@terms[cm], dim = dim(cm), dimnames = dimnames(cm))
}

#' @rdname MICA_term
LCA_depth = function(dag, terms) {
	if(is.character(terms)) {
		id = term_to_node_id(dag, terms, strict = FALSE)
	} else {
		id = terms
	}
	LCA(dag, id, return_term_id = FALSE)
}


Sim_Ling_1998 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = get(IC_method)(dag)[id]

	sim = 2*ic_mica/outer(ic[id], ic[id], "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}

Sim_Resnik_1999 = function(dag, terms, norm_method = "Nmax") {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]
	n_anno = n_annotations(dag)

	if(norm_method == "Nunif") {
		sim = ic_mica/log(attr(n_anno, "N"))
	} else if(norm_method == "Nmax") {
		sim = ic_mica/max(ic)
	} else if(norm_method == "Nunivers") {
		sim = ic_mica/(outer(ic, ic, max))
	} else {
		sim = ic_mica
	}
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Resnik_1999")


Sim_FaITH_2010 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = ic_mica/(outer(ic, ic, "+") - ic_mica)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_FaITH_2010")


Sim_PS_2008 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = 3*ic_mica - outer(ic, ic, "+")
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	diag(sim) = 1
	sim[sim < 0] = 0

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_PS_2008")


Sim_Relevance_2006 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = get_Sim_method("Sim_Ling_1998")(dag, terms)

	eps = 1 - exp(-ic_mica)
	eps*sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Relevance_2006")


Sim_SimIC_2010 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = get_Sim_method("Sim_Ling_1998")(dag, terms)

	eps = 1 - 1/(1 + ic_mica)
	eps*sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_SimIC_2010")



Sim_EICA = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = get_Sim_method("Sim_Ling_1998")(dag, terms)

	eps = cpp_eps_EICA(dag@lt_offspring, dag@lt_ancestor, dag@lt_children, ic)/ic_mica

	eps*sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_EICA")

Sim_XGraSM = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic = calc_IC(dag, IC_method)[id]

	sim = get_Sim_method("Sim_Ling_1998")(dag, terms)

	eps = cpp_eps_XGraSM(dag, id, ic)

	eps*sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_XGraSM")


Sim_AIC_2014 = function(dag, terms) {
	IC_method = "IC_annotation"

	id = term_to_node_id(dag, terms)
	ic = calc_IC(dag, IC_method)[id]

	sim = cpp_common_ancestor_aggregate_aic(dag, id, ic)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	diag(sim) = 1
	sim[sim < 0] = 0

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_AIC_2014")


Sim_Zhang_2006 = function(dag, terms) {
	IC_method = "IC_Zhang_2006"

	id = term_to_node_id(dag, terms)

	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = 2*ic_mica/outer(ic, ic, "+")
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Zhang_2006")


Sim_universal = function(dag, terms) {
	IC_method = "IC_universal"

	id = term_to_node_id(dag, terms)

	ic_mica = MICA_IC(dag, id, IC_method)
	ic = calc_IC(dag, IC_method)[id]

	sim = 2*ic_mica/outer(ic, ic, max)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_universal")


# Sim_Wang_2007 = function(dag) {

# 	g = dag@graph

# 	w = -log(ifelse(E(g)$relation == "isa", 0.8, ifelse(E(g)$relation == "part of", 0.6, 0.7)))

# 	d = distances(g, weights = w, mode = "out")
# 	sv = 1/exp(d)

# 	rm(d); gc(verbose = FALSE)

# 	cpp_ancestor_aggregate_wang(sv, dag@lt_offspring)
# }
# ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Wang_2007")

###########################################
#### edge-based
###########################################

Sim_Rada_1989 = function(dag, terms) {

	id = term_to_node_id(dag, terms)
	dist = dag_distance_undirected(dag, id)
	1/(1 + dist);
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Rada_1989")


Sim_Resnik_edge_2012 = function(dag, terms) {

	id = term_to_node_id(dag, terms)
	max_depth = max(dag_depth(dag))
	dist = dag_distance_undirected(dag, id)

	1 - dist/2/max_depth
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Resnik_edge_2012")


Sim_Leocock_1998 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	max_depth = max(dag_depth(dag))
	dist = dag_distance_undirected(dag, id)

	1 - log(dist)/log(2*max_depth)
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Leocock_1998")


Sim_Wu_1994 = function(dag, terms) {

	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)

	sim =2*lca_depth/dag_longest_distance_undirected_via_LCA(dag, id)

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Wu_1994")

Sim_Slimani_2006 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)

	sim_wp =2*lca_depth/(dag_longest_distance_undirected_via_LCA(dag, id))

	lamda = cpp_is_reachable(dag, id);

	cf = (1 - lambda)*(outer(depth[id], depth[id], min) - dc) +
	     lambda/(abs(outer(depth[id], depth[id], "-")) + 1)

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Slimani_2006")


Sim_Shenoy_2012 = function(dag, terms) {
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)
	max_depth = max(depth)

	dist = dag_longest_distance_undirected_via_LCA(dag, id)
	sim_wp =2*lca_depth/dist

	lamda = cpp_is_reachable(dag, id);

	cf = exp(- lambda*dist/max_depth)

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Shenoy_2012")


Sim_Pekar_2002 = function(dag, terms) {

	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)

	sim = lca_depth/(dag_longest_distance_undirected_via_LCA(dag, id) - lca_depth)

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Pekar_2002")

Sim_Stojanovic_2001 = function(dag, terms) {

	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]

	sim =lca_depth/(outer(depth, depth, "+") + lca_depth)

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Stojanovic_2001")


Sim_Wang_edge_2012 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]

	d1 = dag_longest_distance_undirected_via_LCA(dag, id, 2)
	d2 = dag_longest_distance_undirected_via_LCA(dag, id, 3)

	sim = lca_depth^2/(lca_depth + d1)/(lca_depth + d2)
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Wang_edge_2012")


Sim_Zhong_2002 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]

	k = 2
	x = 1/k^lca_depth
	y = 0.5/k^depth

	sim = 1 - (x - outer(y, y, "+"))

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Zhong_2002")


Sim_AlMubaid_2006 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]
	max_depth = max(depth)

	dist = dag_distance_undirected(dag, id)

	sim = log(1 + (dist - 1)*(max_depth - lca_depth))

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_AlMubaid_2006")


Sim_Li_2003 = function(dag, terms) {
	
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]
	max_depth = max(depth)

	dist = dag_distance_undirected(dag, id)

	alpha = 0.2
	beta = 0.8

	sim = exp(-alpha*dist) * atan(beta*lca_depth)

	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_Li_2003")


###########################################
#### hybrid
###########################################

Sim_RSS = function(dag, terms) {
	id = term_to_node_id(dag, terms)
	lca_depth = LCA_depth(dag, id)
	depth = dag_depth(dag)[id]
	max_depth = max(depth)
	height = dag_height(dag)[id]

	dist = dag_distance_undirected(dag, id)

	sim = max_depth/(max_depth + dist) * lca_depth/(lca_depth + 0.5*outer(height, height, "+"))
	sim
}
ALL_SIM_METHODS = c(ALL_SIM_METHODS, "Sim_RSS_2013")


Sim_HRSS = function(dag, terms, IC_method) {
	id = term_to_node_id(dag, terms)
	ic = calc_IC(dag, IC_method)[id]
	ic_mica = MICA_IC(dag, id, IC_method)

	MIL_term = cpp_most_informative_leaf(dag, id, ic)

	alpha = abs(ic[dag@root] - ic_mica)
	ic_dist_leaf = ic[id] - ic[MIL_term]
	beta = outer(ic_dist_leaf, ic_dist_leaf, "+")/2

	sim = 1/(1 + abs(outer(ic[id], ic[id], "-")))*alpha/(alpha + beta)
	sim
}


Sim_Shen_2010 = function(dag, terms, IC_method) {
	id = term_to_node_id(dag, terms)
	ic = calc_IC(dag, IC_method)[id]

	sim = cpp_sim_shen(dag, id, ic);
	sim
}


Sim_SSDD = function(dag, terms, IC_method) {
	id = term_to_node_id(dag, terms)
	ic = calc_IC(dag, IC_method)[id]
	depth = dag_depth(dag)

	sim = cpp_sim_SSDD(dag, id, depth, ic);
	sim
}



