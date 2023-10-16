

#' Sim_Lin_1998
#' 
#' @section Methods:
#' ## Sim_Lin_1998
#' 
#' The similarity between two terms `a` and `b` is calculated as the IC of their MICA term `c` normalized by the average of the IC of the two terms:
#' 
#' ```
#' sim = IC(c)/((IC(a) + IC(b))/2) 
#'     = 2*IC(c)/(IC(a) + IC(b))
#' ```
#' 
#' Although any IC method can be used here, in more applications, it is normally used together with the *IC_annotation* method.
#' 
#' Paper link: \doi{10.5555/645527.657297}.
#' 
#' @rdname temp__Sim_Lin_1998
Sim_Lin_1998 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Lin_1998")
	}
	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	if(IC_method == "IC_annotation") {
		l = validate_annotated_terms(dag, id)
		id = id[l]
		ic = ic[l]

		if(verbose) {
			if(any(!l)) {
				message(sum(!l), " terms are removed because of no annotation.")
			}
		}
	}

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	sim = 2*ic_mica/outer(ic, ic, "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_Lin_1998")


#' Sim_Resnik_1999
#' 
#' @section Methods:
#' ## Sim_Resnik_1999
#' 
#' The IC method is fixed to `IC_annotation`.
#' 
#' The original Resnik similarity is the IC of the MICA term. There are three ways to normalize the Resnik similarity into the scale of `[0, 1]`:
#' 
#' 1. *Nunif*
#' 
#' ```
#' sim = IC(c)/log(N)
#' ```
#' 
#' where `N` is the total number of items annotated to the whole DAG, i.e. number of items annotated to the root. Then the IC
#' of a term with only one item annotated is `-log(1/N)` = log(N)` which is the maximal IC value in the DAG. 
#' 
#' 2. *Nmax*
#' 
#' `IC_max` is the maximal IC of all terms. If there is a term with only one item annotated, `Nmax` is identical to the `Nunif* method.
#' 
#' ```
#' sim = IC(c)/IC_max
#' ``` 
#' 
#' 3. *Nunivers*
#' 
#' The IC is normalized by the maximal IC of term `a` and `b`.
#' 
#' ```
#' sim = IC(c)/max(IC(a), IC(b))
#' ```
#' 
#' Paper link: \doi{10.1613/jair.514}, \doi{10.1186/1471-2105-9-S5-S4}, \doi{10.1186/1471-2105-11-562}, \doi{10.1155/2013/292063}.
#' 
#' The normalization method can be set with the `norm_method` parameter:
#' 
#' ```
#' term_sim(dag, terms, control = list(norm_method = "Nmax"))
#' ```
#' 
#' Possible values for the `norm_method` parameter are "Nunif", "Nmax", "Nunivers" and "none".
#' 
#' @rdname temp__Sim_Resnik_1999
Sim_Resnik_1999 = function(dag, terms, norm_method = "Nmax", verbose = simona_opt$verbose) {
	IC_method = "IC_annotation"

	if(!norm_method %in% c("Nunif", "Nmax", "Nunivers", "none")) {
		stop("`norm_method` can only be in 'Nunif', 'Nmax', 'Nunivers, 'none'.")
	}

	if(verbose) {
		message("term_sim_method: ", "Sim_Resnik_1999 +", norm_method)
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]
	
	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

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
ADD_TERM_SIM_METHOD("Sim_Resnik_1999", require_anno = TRUE)


#' Sim_FaITH_2010
#' 
#' @section Methods:
#' ## Sim_FaITH_2010
#' 
#' It is calculated as:
#' 
#' ```
#' sim = IC(c)/(IC(a) + IC(b) - IC(c))
#' ```
#' 
#' The relation between *FaITH_2010* similarity and *Lin_1998* similarity is:
#' 
#' ```
#' sim_FaITH = sim_Lin/(2 - sim_Lin)
#' ```
#' 
#' Paper link: \doi{10.1007/978-3-642-17746-0_39}.
#' 
#' @rdname temp__Sim_FaITH_2010
Sim_FaITH_2010 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_FaITH_2010")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	sim = ic_mica/(outer(ic, ic, "+") - ic_mica)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim[is.na(sim)] = 1
	sim
}
ADD_TERM_SIM_METHOD("Sim_FaITH_2010")


#' Sim_Relevance_2006
#' 
#' @section Methods:
#' ## Sim_Relevance_2006
#' 
#' The IC method is fixed to `IC_annotation`.
#' 
#' If thinking *Lin_1998* is a measure of how close term `a` and `b` to their MICA term `c`, the relevance method corrects it by multiplying
#' a factor which considers the specificity of how `c` brings the information. The factor is calculated as `1-p(c)` where `p(c)` is the annotation-based
#' probability `p(c) = k/N` where `k` is the number of items annotated to `c` and `N` is the total number of items annotated to the DAG. Then
#' the Relevance semantic similarity is calculated as:
#' 
#' ```
#' sim = (1 - p(c)) * IC_Lin 
#'     = (1 - p(c)) * 2*IC(c)/(IC(a) + IC(b))
#' ```
#' 
#' Paper link: \doi{10.1186/1471-2105-7-302}.
#' 
#' @rdname temp__Sim_Relevance_2006
Sim_Relevance_2006 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_Relevance_2006")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	sim = 2*ic_mica/outer(ic, ic, "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	eps = 1 - exp(-ic_mica)
	eps*sim
}
ADD_TERM_SIM_METHOD("Sim_Relevance_2006")


#' Sim_SimIC_2010
#' 
#' @section Methods:
#' ## Sim_SimIC_2010
#' 
#' The IC method is fixed to `IC_annotation`.
#' 
#' The SimIC method is an improved correction method of the Relevance method because the latter works bad when `p(c)` is very small. The SimIC
#' correction factor for MICA term `c` is:
#' 
#' ```
#' 1 - 1/(1 + IC(c))
#' ```
#' 
#' Then the similarity is:
#' 
#' ```
#' sim = (1 - 1/(1 + IC(c))) * IC_Lin 
#'     = (1 - 1/(1 + IC(c))) * 2*IC(c)/(IC(a) + IC(b))
#' ```
#' 
#' Paper link: \doi{10.48550/arXiv.1001.0958}.
#' 
#' @rdname temp__Sim_SimIC_2010
Sim_SimIC_2010 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_SimIC_2010")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	sim = 2*ic_mica/outer(ic, ic, "+")
	sim[is.na(sim)] = 1
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	eps = 1 - 1/(1 + ic_mica)
	eps*sim
}
ADD_TERM_SIM_METHOD("Sim_SimIC_2010")


#' Sim_XGraSM_2013
#' 
#' @section Methods:
#' ## Sim_XGraSM_2013
#' 
#' The IC method is fixed to `IC_annotation`.
#' 
#' Being different from the "Relevance" and "SimIC_2010" methods that only use the IC of the MICA term, the *XGraSM_2013* uses IC of all common ancestor terms of `a` and `b`.
#' First it calculates the mean IC of all common ancestor terms with positive IC values:
#' 
#' ```
#' IC_mean = mean_t(IC(t)) where t is an ancestor of both a and b, and IC(t) > 0
#' ```
#' 
#' then similar to the *Lin_1998* method, normalize to the average IC of `a` and `b`:
#' 
#' ```
#' sim = IC_mean*2/(IC(a) + IC(b))
#' ```
#' 
#' Paper link: \doi{10.1186/1471-2105-14-284}.
#' 
#' @rdname temp__Sim_XGraSM_2013
Sim_XGraSM_2013 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_XGraSM_2013")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = verbose)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]
	
	mean_ic = exec_under_message_condition({
		cpp_common_ancestor_mean_IC_XGraSM(dag, id, ic)
	}, verbose = verbose)

	sim = mean_ic/outer(ic[id], ic[id], "+")*2
	sim[is.na(sim)] = 0
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_XGraSM_2013")



#' Sim_EISI_2015
#' 
#' @section Methods:
#' ## Sim_EISI_2015
#' 
#' The IC method is fixed to `IC_annotation`.
#' 
#' It also selects a subset of common ancestors of terms `a` and `b`. It only selects common ancestors which can reach `a` or `b` via one of its child terms
#' that does not belong to the common ancestors. In other words, from the common ancestor, there exist a path where
#' the information is uniquely transmitted to `a` or `b`, not passing the other.
#' 
#' Then the mean IC of the subset common ancestors is calculated and normalized by the *Lin_1998* method.
#' 
#' Paper link: \doi{10.1016/j.gene.2014.12.062}.
#' 
#' @rdname temp__Sim_EISI_2015
Sim_EISI_2015 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_EISI_2015")
	}
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = verbose)[id]

	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]

	mean_ic = exec_under_message_condition({
		cpp_common_ancestor_mean_IC_EISI(dag, id, ic)
	}, verbose = verbose)

	sim = mean_ic/outer(ic[id], ic[id], "+")*2
	sim[is.na(sim)] = 0
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_EISI_2015")


#' Sim_AIC_2014
#' 
#' @section Methods:
#' ## Sim_AIC_2014
#' 
#' It uses the aggregate information content from ancestors. First define the semantic weight (`Sw`) of a term `t` in the DAG:
#' 
#' ```
#' Sw = 1/(1 + exp(-1/IC(t)))
#' ```
#' 
#' Then calculate the aggregation only in the common ancestors and the aggregationn
#' in the ancestors of the two terms `a` and `b` separatedly:
#' 
#' ```
#' SV_{common ancestors} = sum_{t in common ancestors}(Sw(t))
#' SV_a = sum{a' in a's ancestors}(Sw(a'))
#' SV_b = sum{b' in b's ancestors}(Sw(b'))
#' ```
#' 
#' The similarity is calculated as the ratio between the aggregation on the common ancestors and the average on `a`'s ancestors and `b`'s ancestors separatedly.
#' 
#' ```
#' sim = 2*SV_{common_ancestors}/(SV_a + SV_b)
#' ```
#' 
#' Paper link: \doi{10.1109/tcbb.2013.176}.
#' 
#' @rdname temp__Sim_AIC_2014
Sim_AIC_2014 = function(dag, terms, IC_method = "IC_annotation", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_AIC_2014")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = verbose)[id]
	
	l = validate_annotated_terms(dag, id)
	id = id[l]
	ic = ic[l]
	ic[ic == 0] = 0  # get rid of -0

	sim = exec_under_message_condition({
		cpp_sim_aic(dag, id, ic)
	}, verbose = verbose)

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_AIC_2014")


#' Sim_Zhang_2006
#' 
#' @section Methods:
#' ## Sim_Zhang_2006
#' 
#' It uses the *IC_Zhang_2006* IC method and the *Lin_1998* method to calculate similarities:
#' 
#' ```
#' sim = 2*IC_zhang(c)/(IC_zhang(a) + IC_zhang(b))
#' ```
#' 
#' @rdname temp__Sim_Zhang_2006
Sim_Zhang_2006 = function(dag, terms, verbose = simona_opt$verbose) {
	IC_method = "IC_Zhang_2006"

	if(verbose) {
		message("term_sim_method: ", "Sim_Zhang_2006")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)
	
	sim = 2*ic_mica/outer(ic, ic, "+")
	
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_Zhang_2006", require_anno = TRUE)


#' Sim_universal
#' 
#' @section Methods:
#' ## Sim_universal
#' 
#' It uses the *IC_universal* IC method and the *Nunivers* method to calculate similarities:
#' 
#' ```
#' sim = IC_universal(c)/max(IC_universal(a), IC_universal(b))
#' ```
#' 
#' @rdname temp__Sim_universal
Sim_universal = function(dag, terms, verbose = simona_opt$verbose) {
	IC_method = "IC_universal"

	if(verbose) {
		message("term_sim_method: ", "Sim_universal")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	ic = term_IC(dag, IC_method, verbose = FALSE)[id]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	sim = 2*ic_mica/outer(ic, ic, pmax)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_universal")


#' Sim_Wang_2007
#' 
#' @section Methods:
#' ## Sim_Wang_2007
#' 
#' First, S-value of an ancestor term `c` on a term `a` (`S(c->a)`) is calculated (the definition of the S-value can be found in the help page of [`term_IC()`]).
#' Similar to the *Sim_AIC_2014*, aggregation only to common ancestors, to `a`'s ancestors and to `b`'s ancestors are calculated.
#' 
#' ```
#' SV_{common ancestors} = sum_{c in common ancestors}(S(c->a) + S(c->b))
#' SV_a = sum{a' in a's ancestors}(S(a'->a))
#' SV_b = sum{b' in b's ancestors}(S(b'->b))
#' ```
#' 
#' Then the similarity is calculated as:
#' 
#' ```
#' sim = SV_{common_ancestors}*2/(SV_a + SV_b)
#' ```
#' 
#' Paper link: \doi{10.1093/bioinformatics/btm087}.
#' 
#' The contribution of different semantic relations can be set with the `contribution_factor` parameter. The value should be a named numeric
#' vector where names should cover the relations defined in `relations` set in [`create_ontology_DAG()`]. For example, if there are two relations
#' "relation_a" and "relation_b" set in the DAG, the value for `contribution_factor` can be set as:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Wang_2007", 
#'     control = list(contribution_factor = c("relation_a" = 0.8, "relation_b" = 0.6)))
#' ```
#' 
#' @rdname temp__Sim_Wang_2007
#' @import igraph
Sim_Wang_2007 = function(dag, terms, contribution_factor = c("is_a" = 0.8, "part_of" = 0.6), 
	calc_by = "igraph", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Wang_2007")
	}

	if(length(dag@lt_children_relations) == 0) {
		stop("`relations` is not set when creating the ontology_DAG object.")
	}
	relation_levels = attr(dag@lt_children_relations, "levels")
	if(is.null(names(contribution_factor))) {
		stop("`contribution_factor` should be a named numeric vector where names should correspond to all relations.")
	}

	names(contribution_factor) = normalize_relation_type(names(contribution_factor))
	
	contribution_factor = extend_contribution_factor(dag@relations_DAG, contribution_factor)
	if(length(setdiff(relation_levels, names(contribution_factor)))) {
		stop("Contribution factor should be provided for all relations.")
	}

	if(any(contribution_factor >= 1)) {
		stop("All values in `contribution_factor` should be smaller than 1.")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	sim = matrix(0, nrow = length(id), ncol = length(id))
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	diag(sim) = 1

	if(length(id) <= 1) {
		return(sim)
	}

	if(calc_by == "igraph") {
		g = dag_as_igraph(dag)
		E(g)$weight = contribution_factor[E(g)$relation]
		all_ancestors = cpp_ancestors_of_a_group(dag, id, include_self = TRUE)
		d = distances(g, v = all_ancestors, to = id, mode = "out", weights = -log(E(g)$weight))
		s = exp(-d)  # rows are ancestors
		sim = cpp_wang_sv_to_sim(s)
	} else {
		sim = cpp_sim_wang(dag, id, unname(contribution_factor[relation_levels]))
		
	}
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_Wang_2007")


#' Sim_GOGO_2018
#' 
#' @section Methods:
#' ## Sim_GOGO_2018
#' 
#' It is very similar as _Sim_Wang_2007_, but with a corrected contribution factor when calculating the S-value.
#' From a parent term to a child term, _Sim_Wang_2007_ directly uses a weight for the relation between the parent
#' and the child, e.g. 0.8 for "is_a" relation type and 0.6 for "part_of" relation type. In _Sim_GOGO_2018_, the weight
#' is also scaled by the total number of children of that parent:
#' 
#' ```
#' w = 1/(c + nc) + w_0
#' ```
#' 
#' where w_0 is the original contribution factor, `nc` is the number of child terms of the parent, `c` is calculated to ensure that 
#' maximal value of `w` is no larger than 1, i.e. `c = max(w_0)/(1 - max(w_0))`, assuming minimal value of `nc` is 1. By default _Sim_GOGO_2018_
#' sets contribution factor of 0.4 for "is_a" and 0.3 for "part_of", then `w = 1/(2/3 + nc) + w_0`.
#' 
#' Paper link: \doi{10.1038/s41598-018-33219-y}.
#' 
#' The contribution of different semantic relations can be set with the `contribution_factor` parameter. The value should be a named numeric
#' vector where names should cover the relations defined in `relations` set in [`create_ontology_DAG()`]. For example, if there are two relations
#' "relation_a" and "relation_b" set in the DAG, the value for `contribution_factor` can be set as:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_GOGO_2018", 
#'     control = list(contribution_factor = c("relation_a" = 0.4, "relation_b" = 0.3)))
#' ```
#' 
#' @rdname temp__Sim_GOGO_2018
Sim_GOGO_2018 = function(dag, terms, contribution_factor = c("is_a" = 0.4, "part_of" = 0.3), 
	calc_by = "igraph", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_GOGO_2018")
	}

	if(length(dag@lt_children_relations) == 0) {
		stop("`relations` is not set when creating the ontology_DAG object.")
	}
	relation_levels = attr(dag@lt_children_relations, "levels")
	if(is.null(names(contribution_factor))) {
		stop("`contribution_factor` should be a named numeric vector where names should correspond to all relations.")
	}

	names(contribution_factor) = normalize_relation_type(names(contribution_factor))
	
	contribution_factor = extend_contribution_factor(dag@relations_DAG, contribution_factor)
	if(length(setdiff(relation_levels, names(contribution_factor)))) {
		stop("Contribution factor should be provided for all relations.")
	}

	if(any(contribution_factor >= 1)) {
		stop("All values in `contribution_factor` should be smaller than 1.")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	sim = matrix(0, nrow = length(id), ncol = length(id))
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	diag(sim) = 1

	if(length(id) <= 1) {
		return(sim)
	}

	if(calc_by == "igraph") {
		g = dag_as_igraph(dag)
		em = get.edgelist(g)
		c = max(contribution_factor)/(1 - max(contribution_factor))
		nc = n_children(dag)
		E(g)$weight = contribution_factor[E(g)$relation] + 1/(c + nc[em[, 1]])
		all_ancestors = cpp_ancestors_of_a_group(dag, id, include_self = TRUE)
		d = distances(g, v = all_ancestors, to = id, mode = "out", weights = -log(E(g)$weight))
		s = exp(-d)  # rows are ancestors
		sim = cpp_wang_sv_to_sim(s)
	} else {
		sim = cpp_sim_wang(dag, id, unname(contribution_factor[relation_levels]), TRUE)
		
	}
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	
	sim
}
ADD_TERM_SIM_METHOD("Sim_GOGO_2018")



###########################################
#### edge-based
###########################################

#' Sim_Rada_1989
#' 
#' @section Methods:
#' ## Sim_Rada_1989
#' 
#' It is based on the distance between term `a` and `b`. It is defined as:
#' 
#' ```
#' sim = 1/(1 + d(a, b))
#' ```
#' 
#' The distance can be the shortest distance between `a` and `b` or the longest distance via the LCA term.
#' 
#' Paper link: \doi{10.1109/21.24528}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Rada_1989",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_Rada_1989
Sim_Rada_1989 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Rada_1989")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	if(distance == "shortest_distances_via_NCA") {
		dist = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dist = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	sim = 1/(1 + dist)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Rada_1989")


#' Sim_Resnik_edge_2005
#' 
#' @section Methods:
#' ## Sim_Resnik_edge_2005
#' 
#' It is also based on the distance between term `a` and `b`:
#' 
#' ```
#' sim = 1 - d(a, b)/2/max_depth
#' ```
#' 
#' where `max_depth` is the maximal depth (maximal distance from root) in the DAG. Similarly, `d(a, b)` can be the shortest
#' distance or the longest distance via LCA.
#' 
#' Paper link: \doi{10.1145/1097047.1097051}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Resnik_edge_2005",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_Resnik_edge_2005
Sim_Resnik_edge_2005 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Resnik_edge_2005")
	}

	id = term_to_node_id(dag, terms)

	max_depth = max(dag_depth(dag))
	if(distance == "shortest_distances_via_NCA") {
		dist = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dist = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	sim = 1 - dist/2/max_depth
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Resnik_edge_2005")


#' Sim_Leocock_1998
#' 
#' @section Methods:
#' ## Sim_Leocock_1998
#' 
#' It is similar as the *Sim_Resnik_edge_2005* method, but it applies log-transformation on the distance and the depth:
#' 
#' ```
#' sim = 1 - log(d(a, b) + 1)/log(2*max_depth + 1)
#' ```
#' 
#' Paper link: \doi{10.1186/1471-2105-13-261}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Leocock_1998",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_Leocock_1998
Sim_Leocock_1998 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_Leocock_1998")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	max_depth = max(dag_depth(dag))
	if(distance == "shortest_distances_via_NCA") {
		dist = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dist = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	sim = 1 - log(dist + 1)/log(2*max_depth + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Leocock_1998")


#' Sim_WP_1994
#' 
#' @section Methods:
#' ## Sim_WP_1994
#' 
#' It is based on the depth of the LCA term `c` and the longest distance between term `a` and `b`:
#' 
#' ```
#' sim = 2*depth(c)/(len_c(a, b) + 2*depth(c))
#' ```
#' 
#' where `len_c(a, b)` is the longest distance between `a` and `b` via LCA `c`. The denominator in the equation can also be written as:
#' 
#' ```
#' len_c(a, b) + 2*depth(c) = depth(c) + len(c, a) + depth(c) + len(c, b)
#'                          = depth_c(a) + depth_c(b)
#' ```
#' 
#' where `depth_c(a)` is the longest distance from root to `a` passing through `c`.
#' 
#' Paper link: \doi{10.3115/981732.981751}.
#' 
#' @rdname temp__Sim_WP_1994
Sim_WP_1994 = function(dag, terms, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_WP_1994")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	
	lca_depth = LCA_depth(dag, id, verbose = verbose)
	sim = 2*lca_depth/(longest_distances_via_LCA(dag, id, verbose = verbose) + 2*lca_depth)
	
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
#' @section Methods:
#' ## Sim_Slimani_2006
#' 
#' It is a correction of the *Sim_WP_1994* method. The correction factor for term `a` and `b` regarding to their LCA `t` is:
#' 
#' ```
#' CF(a, b) = (1-lambda)*(min(depth(a), depth(b)) - depth(c)) + 
#'            lambda/(1 + abs(depth(a) - depth(b)))
#' ```
#' 
#' `lambda` takes value of 1 if `a` and `b` are in ancestor-offspring relation, or else it takes 0. 
#' 
#' Paper link: <https://zenodo.org/record/1075130>.
#' 
#' @rdname temp__Sim_Slimani_2006
Sim_Slimani_2006 = function(dag, terms, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Slimani_2006")
	}
	
	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id, verbose = verbose)
	depth = dag_depth(dag)

	sim_wp = 2*lca_depth/(longest_distances_via_LCA(dag, id, verbose = verbose) + 2*lca_depth)

	lambda = cpp_is_reachable(dag, id)
	ltd = longest_distances_from_LCA(dag, id, verbose = verbose)
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
#' @section Methods:
#' ## Sim_Shenoy_2012
#' 
#' It is a correction of the *Sim_WP_1994* method. The correction factor for term `a` and `b` is:
#' 
#' ```
#' CF(a, b) = exp(-lambda*d(a, b)/max_depth)
#' ```
#' 
#' `lambda` takes value of 1 if `a` and `b` are in ancestor-offspring relation, or else it takes 0. `d(a, b)
#' 
#' Paper link: \doi{10.48550/arXiv.1211.4709}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Leocock_1998",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_Shenoy_2012
Sim_Shenoy_2012 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {
		
	if(verbose) {
		message("term_sim_method: ", "Sim_Shenoy_2012")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id, verbose = verbose)
	depth = dag_depth(dag)

	if(distance == "longest_distances_via_LCA") {
		sim_wp = 2*lca_depth/(longest_distances_via_LCA(dag, id, verbose = verbose) + 2*lca_depth)
	} else if(distance == "shortest_distances_via_NCA") {
		sim_wp = 2*lca_depth/(shortest_distances_via_NCA(dag, id, verbose = verbose) + 2*lca_depth)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	lambda = cpp_is_reachable(dag, id)

	max_depth = max(dag_depth(dag))
	dist = longest_distances_via_LCA(dag, id, verbose = verbose)

	cf = exp(- lambda*dist/max_depth)
	diag(cf) = 1

	sim = cf * sim_wp
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Shenoy_2012")


#' Sim_Pekar_2002
#' 
#' @section Methods:
#' ## Sim_Pekar_2002
#' 
#' It is very similar to the *Sim_WP_1994* method:
#' 
#' ```
#' sim = depth(c)/(len_c(a, b) + depth(c))
#'     = d(root, c)/(d(c, a) + d(c, b) + d(root, c))
#' ```
#' 
#' where `d(a, b)` is the longest distance between `a` and `b`.
#' 
#' Paper link: <https://aclanthology.org/C02-1090/>.
#' 
#' @rdname temp__Sim_Pekar_2002
Sim_Pekar_2002 = function(dag, terms, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Pekar_2002")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id, verbose = verbose)

	sim = lca_depth/(longest_distances_via_LCA(dag, id, verbose = verbose) + lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Pekar_2002")


#' Sim_Stojanovic_2001
#' 
#' @section Methods:
#' ## Sim_Stojanovic_2001
#' 
#' It is purely based on the depth of term `a`, `b` and their LCA `c`.
#' 
#' ```
#' sim = depth(c)/(depth(a) + depth(b) - depth(c))
#' ```
#' 
#' The similarity value might be negative because there is no restrction that the path from root to `a` or `b` must pass `c`.
#' 
#' Paper link: \doi{10.1145/500737.500762}.
#' 
#' @rdname temp__Sim_Stojanovic_2001
Sim_Stojanovic_2001 = function(dag, terms, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Stojanovic_2001")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = which(id %in% dag@root)
	if(any(l)) {
		if(verbose) message("remove root term.")
		id = id[!l]
	}

	depth = dag_depth(dag)
	lca_term = max_ancestor_id(dag, id, depth, in_labels = FALSE, verbose = verbose)
	lca_depth = structure(depth[lca_term], dim = dim(lca_term))
	
	sim = lca_depth/(outer(depth[id], depth[id], "+") - lca_depth)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Stojanovic_2001")


#' Sim_Wang_edge_2012
#' 
#' @section Methods:
#' ## Sim_Wang_edge_2012
#' 
#' It is calculated as:
#' 
#' ```
#' sim = depth(c)^2/depth_c(a)/depth_c(b)
#' ```
#' 
#' where `depth_c(a)` is the longest distance between root to `a` passing through `c`.
#' 
#' Paper link: \doi{10.1186/1477-5956-10-s1-s18}.
#' 
#' @rdname temp__Sim_Wang_edge_2012
Sim_Wang_edge_2012 = function(dag, terms, verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_Wang_edge_2012")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	sim = cpp_sim_wang_edge(dag, id)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Wang_edge_2012")


#' Sim_Zhong_2002
#' 
#' @section Methods:
#' ## Sim_Zhong_2002
#' 
#' For a term `x`, it first calculates a "mile-stone" value as 
#' 
#' ```
#' m(x) = 0.5/2^depth(x)
#' ```
#' 
#' The the distance bewteen term `a` and `b` via LCA term `c` is:
#' 
#' ```
#' D(c, a) + D(c, b) = m(c) - m(a) + m(c) - m(b)
#'                   = 2*m(c) - m(a) - m(b)
#'                   = 1/2^depth(c) - 0.5/2^depth(a) - 0.5/2^depth(b)
#' ```
#' 
#' We change the original `depth(a)` to let it go through LCA term `c` when calculating the depth:
#' 
#' ```
#' 1/2^depth(c) - 0.5/2^depth(a) - 0.5/2^depth(b) 
#'     = 1/2^depth(c)- 0.5/2^(depth(c) + len(c, a)) - 0.5/2^(depth(c) + len(c, b))
#'     = 1/2^depth(c) * (1 - 1/2^(len(c, a) + 1) - 1/2^(len(c, b) + 1))
#'     = 2^-depth(c) * (1 - 2^-(len(c, a) + 1) - 2^-(len(c, b) + 1))
#' ```
#' 
#' And the final similarity is `1 - distance`:
#' 
#' ```
#' sim = 1 - 2^-depth(c) * (1 - 2^-(len(c, a) + 1) - 2^-(len(c, b) + 1))
#' ```
#' 
#' Paper link: \doi{10.1007/3-540-45483-7_8}.
#' 
#' There is a parameter `depth_via_LCA` that can be set to `TRUE` or `FALSE`. IF it is set to `TRUE`, `depth(a)` is re-defined
#' as should pass the LCA term `c`. If it is `FALSE`, it goes to the original similarity definition in the paper and note the 
#' similarity might be negative.
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Zhong_2002",
#'     control = list(depth_via_LCA = FALSE))
#' ```
#' 
#' @rdname temp__Sim_Zhong_2002
Sim_Zhong_2002 = function(dag, terms, depth_via_LCA = TRUE, verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_Zhong_2002")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	if(!depth_via_LCA) {
		d = dag_depth(dag)
		sim = 1 - ( 2^(-LCA_depth(dag, id, verbose = verbose)) - 0.5*outer(2^(-d[id]), 2^(-d[id]), "+") )
	} else {
		sim = cpp_sim_zhong(dag, id, depth_via_LCA)
	}

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Zhong_2002")


#' Sim_AlMubaid_2006
#' 
#' @section Methods:
#' ## Sim_AlMubaid_2006
#' 
#' It also takes accout of the distance between term `a` and `b`, and the depth of the LCA term `c` in the DAG.
#' The distance is calculated as:
#' 
#' ```
#' D(a, b) = log(1 + d(a, b)*(max_depth - depth(c)))
#' ```
#' 
#' Here `d(a, b)` can be the shortest distance between `a` and `b` or the longst distance via LCA `c`.
#' 
#' Then the distance is transformed into the similarity value scaled by the possible maximal and minimal values of `D(a, b)` from the DAG:
#' 
#' ```
#' D_max = log(1 + 2*max_depth * max_depth)
#' ```
#' 
#' And the minimal value of `D(a, b)` is zero when `a` is identical to `b`. Then the similarity value is scaled as:
#' 
#' ```
#' sim = 1 - D(a, b)/D_max
#' ```
#' 
#' Paper link: \doi{10.1109/IEMBS.2006.259235}.
#'
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_AlMubaid_2006",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_AlMubaid_2006
Sim_AlMubaid_2006 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_AlMubaid_2006")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id, verbose = verbose)
	max_depth = max(dag_depth(dag))

	if(distance == "shortest_distances_via_NCA") {
		dsp = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	dist = log(1 + dsp*(max_depth - lca_depth))
	max_dist = log(1 + 2*max_depth*max_depth)
	sim = 1 - dist/max_dist
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_AlMubaid_2006")


#' Sim_Li_2003
#' 
#' @section Methods:
#' ## Sim_Li_2003
#' 
#' It is similar to the *Sim_AlMubaid_2006* method, but uses a non-linear form:
#' 
#' ```
#' sim = exp(0.2*d(a, b)) * atan(0.6*depth(c))
#' ```
#' 
#' where `d(a, b)` can be the shortest distance or the longest distance via LCA.
#' 
#' Paper link: \doi{10.1109/TKDE.2003.1209005}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Li_2003",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_Li_2003
Sim_Li_2003 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {
	
	if(verbose) {
		message("term_sim_method: ", "Sim_Li_2003")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	lca_depth = LCA_depth(dag, id, verbose = verbose)

	if(distance == "shortest_distances_via_NCA") {
		dsp = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
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
#' @section Methods:
#' ## Sim_RSS_2013
#' 
#' The similarity is adjusted by the positions of term `a`, `b` and the LCA term `c` in the DAG. The similarity is defined as:
#' 
#' ```
#' sim = max_depth/(max_depth + d(a, b)) * alpha/(alpha + beta)
#' ```
#' 
#' where `d(a, b)` is the distance between `a` and `b` which can be the shortest distance or the longest distance via LCA.
#' 
#' In the tuning factor, `alpha` is the distance of LCA to root, which is `depth(c)`. `beta` is the distance to leaves, which
#' is the minimal distance (or the minimal height) of term `a` and `b`:
#' 
#' ```
#' alpha/(alpha + beta) = depth(c)/(depth(c) + min(height(a), height(b)))
#' ```
#' 
#' Paper link: \doi{10.1371/journal.pone.0066745}.
#' 
#' There is a parameter `distance` which takes value of "longest_distances_via_LCA" (the default) or "shortest_distances_via_NCA":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_RSS_2013",
#'     control = list(distance = "shortest_distances_via_NCA"))
#' ```
#' 
#' @rdname temp__Sim_RSS_2013
Sim_RSS_2013 = function(dag, terms, distance = "longest_distances_via_LCA", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_RSS_2013")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	lca_depth = LCA_depth(dag, id, verbose = verbose)
	max_depth = max(dag_depth(dag))
	height = dag_height(dag)[id]

	if(distance == "shortest_distances_via_NCA") {
		dsp = shortest_distances_via_NCA(dag, id, verbose = verbose)
	} else if(distance == "longest_distances_via_LCA") {
		dsp = longest_distances_via_LCA(dag, id, verbose = verbose)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_NCA' or 'longest_distances_via_LCA'.")
	}

	sim = max_depth/(max_depth + dsp) * lca_depth/(lca_depth + outer(height, height, pmin) + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_RSS_2013")


#' Sim_HRSS_2013
#' 
#' @section Methods:
#' ## Sim_HRSS_2013
#' 
#' It is similar as the *Sim_RSS_2013* method, but it uses information content instead of the distance to adjust the similarity.
#' 
#' It first defines the semantic distance between term `a` and `b` as the sum of the distance to their MICA term `c`:
#' 
#' ```
#' D(a, b) = D(c, a) + D(c, b)
#' ```
#' 
#' And the distance between an ancestor to a term is:
#' 
#' ```
#' D(c, a) = IC(a) - IC(c)  # if c is an ancestor of a
#' D(a, b) = D(c, a) + D(c, b) = IC(a) + IC(b) - 2*IC(c) # if c is the MICA of a and b
#' ```
#' 
#' Similarly, the similarity is also corrected by the position of MICA term and `a` and `b` in the DAG:
#' 
#' ```
#' 1/(1 + D(a, b)) * alpha/(alph + beta)
#' ```
#' 
#' Now `alpha` is the IC of the MICA term:
#' 
#' ```
#' alpha = IC(c)
#' ```
#' 
#' And `beta` is the average of the maximal semantic distance of `a` and `b` to leaves.
#' 
#' ```
#' beta = 0.5*(IC(l_a) - IC(a) + IC(l_b) - IC(b))
#' ```
#' 
#' where `l_a` is the leaf that `a` can reach with the highest IC (i.e. most informative leaf), and so is `l_b`.
#' 
#' Paper link: \doi{10.1371/journal.pone.0066745}.
#' 
#' @rdname temp__Sim_HRSS_2013
Sim_HRSS_2013 = function(dag, terms, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_HRSS_2013")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	IC_method = "IC_annotation"
	if(IC_method == "IC_annotation") {
		l = validate_annotated_terms(dag, id)
		id = id[l]
	}
	ic = term_IC(dag, IC_method, verbose = FALSE)
	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)

	MIL_term = cpp_max_leaves_id(dag, id, ic)

	alpha = abs(ic[dag@root] - ic_mica)
	ic_dist_leaf = abs(ic[id] - ic[MIL_term])
	beta = outer(ic_dist_leaf, ic_dist_leaf, "+")/2

	sim = 1/(1 + abs(outer(ic[id], ic[id], "-")))*alpha/(alpha + beta + 1)
	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_HRSS_2013", require_anno = TRUE)


#' Sim_Shen_2010
#' 
#' @section Methods:
#' ## Sim_Shen_2010
#' 
#' It is based on the information content of terms on the path connecting term `a` and `b` via their MICA term `c`.
#' 
#' Denote a list of terms `a, ..., c, ..., b` which are composed by the shortest path from `a` to `c` and from `b` to `c`, the difference
#' between `a` and `b` is the sum of `1/IC` of the terms on the path:
#' 
#' ```
#' sum_{x in the path}(1/IC(x))
#' ```
#' 
#' Then the distance is scaled into `[0, 1]` by an arctangent tarnsformation:
#' 
#' ```
#' atan(sum_{x in the path}(1/IC(x)))/(pi/2)
#' ```
#' 
#' And finally the similarity is:
#' 
#' ```
#' sim = 1 - atan(sum_{x in the path}(1/IC(x)))/(pi/2)
#' ``` 
#' 
#' Paper link: \doi{10.1109/BIBM.2010.5706623}.
#' 
#' @rdname temp__Sim_Shen_2010
Sim_Shen_2010 = function(dag, terms, IC_method = "IC_annotation", distance = "shortest_distances_via_NCA", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Shen_2010")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	if(IC_method == "IC_annotation") {
		l = validate_annotated_terms(dag, id)
		id = id[l]
	}
	ic = term_IC(dag, IC_method, verbose = FALSE)

	sim = max_ancestor_path_sum(dag, id, ic, 1/ic, distance = ifelse(distance == "longest_distances_via_LCA", "longest", "shortest"), verbose = verbose)
	sim = 1 - atan(sim)/pi*2

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Shen_2010")


#' Sim_SSDD_2013
#' 
#' @section Methods:
#' ## Sim_SSDD_2013
#' 
#' It is similar as the *Sim_Shen_2010* which also sums content along the path passing through LCA term.
#' Instead of summing the information content, the *Sim_SSDD_2013* sums up a so-called "T-value":
#' 
#' ```
#' sim = 1 - atan(sum_{x in the path}(T(x)))/(pi/2)
#' ``` 
#' 
#' Each term has a T-value and it measures the semantic content a term averagely inherited from its parents
#' and distributed to its offsprings. The T-value of root is 1. Assume a term `t` has two parents `p1` and `p1`,
#' The T-value for term `t` is averaged from its
#' 
#' ```
#' (w1*T(p1) + w2*T(p2))/2
#' ```
#' 
#' Since the parent may have other child terms, a factor `w1` or `w2` is multiplied to `T(p1)` and `T(p2)`. Taking
#' `p1` as an example, it has `n_p` offsprings (including itself) and `t` has `n_t` offsprings (including itself),
#' this means `n_t/n_p` of information is transmitted from `p1` to downstream via `t`, thus `w1` is defined as `n_t/n_p`.
#' 
#' Paper link: \doi{10.1016/j.ygeno.2013.04.010}.
#' 
#' @rdname temp__Sim_SSDD_2013
Sim_SSDD_2013 = function(dag, terms, distance = "shortest_distances_via_NCA", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_SSDD_2013")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)
	t = totipotency(dag)

	sim = 1 - atan(max_ancestor_path_sum(dag, id, dag_depth(dag), t, distance = ifelse(distance == "longest_distances_via_LCA", "longest", "shortest"), verbose = verbose))/pi*2
	dimnames(sim) = list(dag@terms[id], dag@terms[id])
	sim
}
ADD_TERM_SIM_METHOD("Sim_SSDD_2013")


#' Sim_Jiang_1997
#' 
#' @section Methods:
#' ## Sim_Jiang_1997
#' 
#' First semantic distance between term `a` and `b` via MICA term `c` is defined as:
#' 
#' ```
#' D(a, b) = IC(a) + IC(b) - 2*IC(c)
#' ```
#' 
#' Then there are several normalization method to change the distance to similarity and to scale it into the range of `[0, 1]`.
#' 
#' - max: `1 - D(a, b)/2/IC_max`
#' - Couto: `min(1, D(a, b)/IC_max)`
#' - Lin: `1 - D(a, b)/(IC(a) + IC(b))` which is the same as the *Sim_Lin_1998* method
#' - Garla: `1 - log(D(a, b) + 1)/log(2*IC_max + 1)`
#' - log-Lin: `1 - log(D(a, b) + 1)/log(IC(a) + IC(b) + 1)`
#' - Rada: `1/(1 + D(a, b))`
#' 
#' Paper link: <https://aclanthology.org/O97-1002/>.
#' 
#' There is a parameter `norm_method` which takes value in "max", "Couto", "Lin", "Carla", "log-Lin", "Rada":
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Jiang_1997",
#'     control = list(norm_method = "Lin"))
#' ```
#' 
#' @rdname temp__Sim_Jiang_1997
Sim_Jiang_1997 = function(dag, terms, IC_method = "IC_annotation", norm_method = "max", verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Jiang_1997 + ", norm_method)
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]

	ic_mica = MICA_IC(dag, id, IC_method, verbose = verbose)
	ic = term_IC(dag, IC_method, verbose = FALSE)
	max_ic = max(ic, na.rm = TRUE)  # IC_annotation generates NA
	ic = ic[id]

	dist = outer(ic, ic, "+") - 2*ic_mica
	if(norm_method == "max") {
		sim = 1 - dist/2/max_ic
	} else if(norm_method == "Couto") {
		sim = dist/max_ic
		sim[sim > 1] = 1
		sim = 1 - sim
	} else if(norm_method == "Lin") {
		sim = 2*ic_mica/outer(ic, ic, "+")
	} else if(norm_method == "Garla") {
		sim = 1 - log(dist + 1)/log(2*max_ic + 1)
	} else if(norm_method == "log-Lin") {
		sim = 1 - log(dist + 1)/log(outer(ic, ic, "+") + 1)
	} else if(norm_method == "Rada") {
		sim = 1/(1 + dist)
	} else {
		stop("wrong norm_method")
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
#' @section Methods:
#' ## Sim_Kappa
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is [the kappa coeffcient](https://en.wikipedia.org/wiki/Cohen%27s_kappa)
#' of the two sets. 
#' 
#' The universe or the background can be set via parameter `anno_universe`:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_kappa",
#'     control = list(anno_universe = ...))
#' ```
#' 
#' @rdname temp__Sim_Kappa
Sim_Kappa = function(dag, terms, anno_universe = NULL, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Kappa")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "kappa")
}
ADD_TERM_SIM_METHOD("Sim_Kappa", require_anno = TRUE)


#' Sim_Jaccard
#' 
#' @section Methods:
#' ## Sim_Jaccard
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is the Jaccard coeffcient
#' of the two sets, defined as `length(intersect(A, B))/length(union(A, B))`.
#' 
#' The universe or the background can be set via parameter `anno_universe`:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Jaccard",
#'     control = list(anno_universe = ...))
#' ```
#' 
#' @rdname temp__Sim_Jaccard
Sim_Jaccard = function(dag, terms, anno_universe = NULL, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Jaccard")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "jaccard")
}
ADD_TERM_SIM_METHOD("Sim_Jaccard", require_anno = TRUE)


#' Sim_Dice
#' 
#' @section Methods:
#' ## Sim_Dice
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is the Dice coeffcient
#' of the two sets, defined as `2*length(intersect(A, B))/(length(A) + length(B))`.
#' 
#' The universe or the background can be set via parameter `anno_universe`:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Dice",
#'     control = list(anno_universe = ...))
#' ```
#' 
#' @rdname temp__Sim_Dice
Sim_Dice = function(dag, terms, anno_universe = NULL, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Dice")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "dice")
}
ADD_TERM_SIM_METHOD("Sim_Dice", require_anno = TRUE)


#' Sim_Overlap
#' 
#' @section Methods:
#' ## Sim_Overlap
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is the overlap coeffcient
#' of the two sets, defined as `length(intersect(A, B))/min(length(A), length(B))`.
#' 
#' The universe or the background can be set via parameter `anno_universe`:
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Overlap",
#'     control = list(anno_universe = ...))
#' ```
#' 
#' @rdname temp__Sim_Overlap
Sim_Overlap = function(dag, terms, anno_universe = NULL, verbose = simona_opt$verbose) {

	if(verbose) {
		message("term_sim_method: ", "Sim_Overlap")
	}

	id = term_to_node_id(dag, terms, strict = FALSE)

	l = validate_annotated_terms(dag, id)
	id = id[l]
	.sim_overlap(dag, id, anno_universe, method = "overlap")
}
ADD_TERM_SIM_METHOD("Sim_Overlap", require_anno = TRUE)


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
	m1 = Matrix::rowSums(m)
	m2 = abs(Matrix::rowSums(m - 1))
	aab = (outer(m1, m1) + outer(m2, m2))/tab/tab
	k = (oab - aab)/(1 - aab)
	return(k)
}

overlap_dist = function(m) {
	n = Matrix::rowSums(m)
	proxyC::simil(m, method = "dice")*outer(n, n, "+")/2/outer(n, n, pmin)
}

##############

#' Sim_Ancestor
#' 
#' @section Methods:
#' ## Sim_Ancestor
#' 
#' Denote `S_a` and `S_b` are two sets of ancestor terms of term `a` and `b` (including `a` and `b`), the 
#' semantic similarity is defined as:
#' 
#' ```
#' length(intersect(S_a, S_b))/length(union(S_a, S_b))
#' ```
#' 
#' ```
#' term_sim(dag, terms, method = "Sim_Ancestor")
#' ```
#' 
#' @rdname temp__Sim_Ancestor
Sim_Ancestor = function(dag, terms, verbose = simona_opt$verbose) {
	id = term_to_node_id(dag, terms, strict = FALSE)
	sim = cpp_sim_ancestor(dag, id)

	dimnames(sim) = list(dag@terms[id], dag@terms[id])

	sim
}
ADD_TERM_SIM_METHOD("Sim_Ancestor")
