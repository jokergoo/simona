
validate_termID = function(dag, termID) {

	if(is.null(termID)) {
		return(NULL)
	}

	if(is.numeric(termID)) {  # internally used
		return(termID)
	}

	termID = which(dag@terms %in% termID)

	if(length(termID) == 0) {
		stop("Cannot find terms in `termID`.")
	}

	termID
}

###########################################
#### node-based
###########################################

LCA = function(dag, termID = NULL, return_termID = FALSE) {
	lt_offspring = dag@lt_offspring
	depth = dag_depth(dag)

	if(!is.null(termID)) {
		if(!is.numeric(termID)) termID = which(dag@terms %in% termID)
		lt_offspring = lapply(lt_offspring, intersect, termID)
	}
	if(return_termID) {
		cm = cpp_common_ancestor_ID(lt_offspring, depth, dag@lt_parents, dag@dist_offspring)
	} else {
		cm = cpp_common_ancestor_int(lt_offspring, depth, dag@lt_parents)
	}
	dimnames(cm) = list(dag@terms, dag@terms)

	if(!is.null(termID)) {
		cm[termID, termID, drop = FALSE]
	} else {
		cm
	}
}

LCA_term = function(dag, termID = NULL) {
	LCA(dag, termID, return_termID = TRUE)
}

LCA_depth = function(dag, termID = NULL) {
	LCA(dag, termID, return_termID = FALSE)
}

MICA = function(dag, IC_method, termID = NULL, return_termID = FALSE) {
	lt_offspring = dag@lt_offspring

	if(!is.null(termID)) {
		if(!is.numeric(termID)) termID = which(dag@terms %in% termID)
		lt_offspring = lapply(lt_offspring, intersect, termID)
	}

	if(!IC_method %in% IC_methods) {
		stop("IC method: '", IC_method, "' is not supported.")
	}

	ic = get(IC_method)(dag)

	if(any(ic < 0)) {
		stop("IC cannot be negative.")
	}

	if(return_termID) {
		cm = cpp_common_ancestor_ID(lt_offspring, ic, dag@lt_parents, dag@dist_offspring)
	} else {
		cm = cpp_common_ancestor_double(lt_offspring, ic, dag@lt_parents)
	}
	dimnames(cm) = list(dag@terms, dag@terms)

	if(!is.null(termID)) {
		cm[termID, termID, drop = FALSE]
	} else {
		cm
	}
}

MICA_term = function(dag, IC_method, termID = NULL) {
	MICA(dag, IC_method, termID, return_termID = TRUE)
}

MICA_IC = function(dag, IC_method, termID = NULL) {
	MICA(dag, IC_method, termID, return_termID = FALSE)
}


Sim_Ling_1998 = function(dag, termID = NULL) {
	IC_method = "IC_annotation"

	termID = validate_termID(dag, termID)

	ic_mica = MICA_IC(dag, IC_method, termID = termID)
	ic = get(IC_method)(dag)

	if(is.null(termID)) {
		sim = 2*ic_mica/outer(ic, ic, "+")
		dimnames(sim) = list(dag@terms, dag@terms)
	} else {
		sim = 2*ic_mica/outer(ic[termID], ic[termID], "+")
		dimnames(sim) = list(dag@terms[termID], dag@terms[termID])
	}
	sim[dag@root, dag@root] = 1
	sim
}

Sim_Resnik_1999 = function(dag, termID = NULL, norm_method = "Nmax") {
	IC_method = "IC_annotation"

	termID = validate_termID(dag, termID)
	ic_mica = MICA_IC(dag, IC_method, termID = termID)
	ic = get(IC_method)(dag)

	norm_method = match.arg(norm_method)
	if(is.null(termID)) {
		if(norm_method == "Nunif") {
			sim = ic_mica/log(attr(ic, "N"))
		} else if(norm_method == "Nmax") {
			sim = ic_mica/max(ic)
		} else if(norm_method == "Nunivers") {
			sim = ic_mica/(outer(ic, ic, max))
		} else {
			sim = ic_mica
		}
		dimnames(sim) = list(dag@terms, dag@terms)
	} else {
		if(norm_method == "Nunif") {
			sim = ic_mica/log(attr(ic, "N"))
		} else if(norm_method == "Nmax") {
			sim = ic_mica/max(ic)
		} else if(norm_method == "Nunivers") {
			sim = ic_mica/(outer(ic[termID], ic[termID], max))
		} else {
			sim = ic_mica
		}
		dimnames(sim) = list(dag@terms[termID], dag@terms[termID])
	}
	sim[dag@root, dag@root] = 1
	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Resnik_1999")


Sim_FaITH_2010 = function(dag, termID = NULL) {
	IC_method = "IC_annotation"

	termID = validate_termID(dag, termID)
	ic_mica = MICA_IC(dag, IC_method, termID = termID)
	ic = get(IC_method)(dag)

	if(is.null(termID)) {
		sim = ic_mica/(outer(ic, ic, "+") - ic_mica)
		dimnames(sim) = list(dag@terms, dag@terms)
	} else {
		sim = ic_mica/(outer(ic[termID], ic[termID], "+") - ic_mica)
		dimnames(sim) = list(dag@terms[termID], dag@terms[termID])
	}
	sim[dag@root, dag@root] = 1
	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_FaITH_2010")


Sim_PS_2008 = function(dag, IC_method) {
	IC_method = "IC_annotation"

	termID = validate_termID(dag, termID)
	ic_mica = MICA_IC(dag, IC_method, termID = termID)
	ic = get(IC_method)(dag)

	if(is.null(termID)) {
		sim = 3*ic_mica - outer(ic, ic, "+")
		dimnames(sim) = list(dag@terms, dag@terms)
	} else {
		sim = 3*ic_mica - outer(ic[termID], ic[termID], "+")
		dimnames(sim) = list(dag@terms[termID], dag@terms[termID])
	}
	diag(sim) = 1
	sim[sim < 0] = 0

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_PS_2008")


Sim_Relevance_2006 = function(dag) {
	IC_method = "IC_annotation"

	sim = get(Sim_method)(dag, IC_method)
	ic_mica = MICA_IC(dag, IC_method, FALSE)

	eps = 1 - exp(-ic_mica)
	eps*sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Relevance_2006")


Sim_SimIC_2010 = function(dag) {
	IC_method = "IC_annotation"

	sim = get(Sim_method)(dag, IC_method)
	ic_mica = MICA_IC(dag, IC_method, FALSE)

	eps = 1 - 1/(1 + ic_mica)
	eps*sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_SimIC_2010")



Sim_EICA = function(dag) {
	IC_method = "IC_annotation"

	ic = get(IC_method)(dag)
	sim = get(Sim_method)(dag, IC_method)
	ic_mica = MICA_IC(dag, IC_method, FALSE)

	eps = cpp_eps_EICA(dag@lt_offspring, dag@lt_ancestor, dag@lt_children, ic)/ic_mica

	eps*sim
}

Sim_XGraSM = function(dag) {
	IC_method = "IC_annotation"

	ic = get(IC_method)(dag)
	sim = get(Sim_method)(dag, IC_method)
	ic_mica = MICA_IC(dag, IC_method, FALSE)

	eps = cpp_eps_XGraSM(dag@lt_offspring, ic)/ic_mica

	eps*sim
}

Sim_AIC_2014 = function(dag, IC_method) {
	lt_offspring = dag@lt_offspring
	lt_ancestor = dag@lt_ancestor
	ic = get(IC_method)(dag)

	sw = 1/(1 + exp(-1/ic))
	sv = sapply(lt_ancestor, function(x) sum(sw[x]))
	cpp_ancestor_aggregate_aic(lt_offspring, sv, sw)
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_AIC_2014")


Sim_Zhang = function(dag) {
	IC_method = "IC_Zhang_2006"

	ic_mica = MICA_IC(dag, IC_method, FALSE)
	ic = get(IC_method)(dag)

	2*ic_mica/outer(ic, ic, "+")
}

Sim_universal = function(dag) {
	IC_method = "IC_universal"

	ic_mica = MICA_IC(dag, IC_method, FALSE)
	ic = get(IC_method)(dag)

	2*ic_mica/outer(ic, ic, "max")
}

Sim_Wang_2007 = function(dag) {

	g = dag@graph

	w = -log(ifelse(E(g)$relation == "isa", 0.8, ifelse(E(g)$relation == "part of", 0.6, 0.7)))

	d = distances(g, weights = w, mode = "out")
	sv = 1/exp(d)

	rm(d); gc(verbose = FALSE)

	cpp_ancestor_aggregate_wang(sv, dag@lt_offspring)
}

###########################################
#### edge-based
###########################################

Sim_Rada_1989 = function(dag, termID = NULL) {

	termID = validate_termID(termID, dag)
	g = dag@graph

	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}
	d = distances(g, v = termID, to = termID, mode = "all")
	sim = 1/(1+d)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Rada_1989")


Sim_Resnik_edge_2012 = function(dag, termID = NULL) {

	termID = validate_termID(termID, dag)
	g = dag@graph

	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}
	d = distances(g, v = termID, to = termID, mode = "all")

	max_depth = max(dag_depth(dag))

	sim = 1 - d/2/max_depth
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Resnik_edge_2012")


Sim_Leocock_1998 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	g = dag@graph

	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}
	d = distances(g, v = termID, to = termID, mode = "all")

	max_depth = max(dag_depth(dag))

	sim = 1 - log(d)/log(2*max_depth)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Leocock_1998")


Sim_Wu_1994 = function(dag, termID = NULL) {

	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	lca_term = LCA_term(dag, termID)
	d = distances(dag, termID, termID, mode = "out")

	sim = cpp_sim_pekar(lca_term, d, depth, 2)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Wu_1994")

Sim_Slimani_2006 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)

	dc = LCA_depth(dag, termID)
	sim_wp = 2*dc/(outer(depth[termID], depth[termID], "+"))

	lamda = is.finite(distances(dag, termID, termID, mode = "out")) + 0

	cf = (1 - lambda)*(outer(depth[termID], depth[termID], min) - dc) +
	     lambda/(abs(outer(depth[termID], depth[termID], "-")) + 1)

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Slimani_2006")


Sim_Shenoy_2012 = function(dag, termID = NULL) {
	termID = validate_termID(termID, dag)

	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	dc = LCA_depth(dag, termID)
	sim_wp = 2*dc/(outer(depth[termID], depth[termID], "+"))

	lamda = is.finite(distances(dag, termID, termID, mode = "out")) + 0

	d = distances(dag, termID, termID, mode = "all")
	cf = exp(-(lambda*d)/max(depth))

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Shenoy_2012")


Sim_Pekar_2002 = function(dag, termID = NULL) {

	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	lca_term = LCA_term(dag, termID)
	d = distances(dag, termID, termID, mode = "out")

	sim = cpp_sim_pekar(lca_term, d, depth, 1)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Pekar_2002")

Sim_Stojanovic_2001 = function(dag, termID = NULL) {

	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	dc = LCA_depth(dag, termID)

	sim = dc/(outer(depth[termID], depth[termID], "+") - dc)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Stojanovic_2001")


Sim_Wang_edge_2012 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	lca_term = LCA_term(dag, termID)
	d = distances(dag, termID, termID, mode = "out")

	sim = cpp_sim_wang_edge(lca_term, d, depth)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Wang_edge_2012")


Sim_Zhong_2002 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	dc = LCA_depth(dag, termID)

	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)

	k = 2
	sim = 1 - (1/k^dc - outer(0.5/depth[termID], 0.5/depth[termID], "+"))

	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Zhong_2002")


Sim_AlMubaid_2006 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	dc = LCA_depth(dag, termID)

	g = dag@graph
	d = distances(g, termID, termID, mode = "all")

	max_depth = max(depth)
	sim = log(1 + (d - 1)*(max_depth - dc))/log(1+(2*max_depth - 1)*max_depth)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_AlMubaid_2006")


Sim_Li_2003 = function(dag, termID = NULL) {
	
	termID = validate_termID(termID, dag)
	
	if(is.null(termID)) {
		termID = seq_along(dag@terms)
	}

	depth = dag_depth(dag)
	dc = LCA_depth(dag, termID)

	g = dag@graph
	d = distances(g, termID, termID, mode = "all")

	alpha = 0.2
	beta = 0.8

	sim = exp(-alpha*d) * atan(beta*dc)
	dimnames(sim) = list(dag@terms[termID], dag@terms[termID])

	sim
}
ALL_Sim_methods = c(ALL_Sim_methods, "Sim_Li_2003")


###########################################
#### hybrid
###########################################

Sim_RSS = function(dag) {
	lca_term = LCA_term(dag)
	
	height = dag_height(dag)
	dist_to_leaf = sapply(dag@lt_offspring, function(x) {
		max(height[x])
	}) - height

	depth = dag_depth(dag)
	max_depth = max(depth)

	d = distances(dag@graph, mode = "all")
	alpha = LCA_depth(dag)
	beta = outer(dist_to_leaf, dist_to_leaf, "+")/2
	max_depth/(max_depth + d) * alpha/(alpha + beta)
}

Sim_HRSS = function(dag) {
	lt_offspring = dag@lt_offspring
	most_informative_leaf = sapply(lt_offspring, function(x) {
		max(ic[x])
	})
	ic_diff = abs(ic - most_informative_leaf)

	1/(1 + abs(outer(ic, ic, "-")))*ic_mica/(ic_mica + outer(ic_diff, ic_diff, "+")/2)
}


Sim_Shen = function(dag) {
	mica = common_ancestor(lt_offspring, ic, TRUE)
	d = distances(dag@graph, weight = 1/ic, FALSE)

	cpp_sim_shen(mica, d);
}




Sim_SSDD = function(dag) {
	mica = common_ancestor(lt_offspring, depth(dag), FALSE)
	t = totipotency_bfs(dag)

	g = dag@graph

	sp = list()
	lt_offspring = dag@lt_offspring
	for(i in seq_along(lt_offspring)) {
		o = lt_offspring[[i]][-1]
		if(length(o)) {
			sp[[i]] = as.vector(shortest_paths(g, i, o))
		}
	}

	cpp_sim_ssdd(mica, t, sp)
}




