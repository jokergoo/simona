

# cat(paste("#' @inheritSection",ALL_IC_METHODS, "method") , sep = "\n")

#' Information content
#' 
#' @param dag An `ontology_DAG` object.
#' @param method An IC method. All available methods are in [`ALL_IC_METHODS`].
#' @param control A list of parameters passing to individual methods. See the subsections.
#' 
#' @inheritSection IC_annotation method
#' @inheritSection IC_universal method
#' @inheritSection IC_Zhang_2006 method
#' @inheritSection IC_Seco_2004 method
#' @inheritSection IC_Zhou_2008 method
#' @inheritSection IC_Seddiqui_2010 method
#' @inheritSection IC_Sanchez_2011 method
#' @inheritSection IC_Meng_2012 method
#' @inheritSection IC_Wang_2007 method
#' 
#' @return A numeric vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = 1:3,
#'     "b" = 3:4,
#'     "c" = 5,
#'     "d" = 7,
#'     "e" = 4:7,
#'     "f" = 8
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' term_IC(dag, "IC_annotation")
term_IC = function(dag, method, control = list()) {
	IC_fun = get_IC_method(method)
	if(method == "IC_Wang_2007") {
		if(!is.null(control$contribution_factor)) {
			ic = IC_fun(dag, contribution_factor = control$contribution_factor)
		}
	} else {
		ic = IC_fun(dag)
	}
	structure(ic, names = dag@terms)
}

get_IC_method = function(method) {
	if(!method %in% ALL_IC_METHODS) {
		stop("Supported IC method should be in `ALL_IC_METHODS`")
	}

	get(method, envir = topenv(), inherits = FALSE)
}


# cat(paste("#' @inheritSection",ALL_TERM_SIM_METHODS, "method") , sep = "\n")

#' Semantic similarity
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param method A term similarity method. All available methods are in [`ALL_TERM_SIM_METHODS`].
#' @param control A list of parameters passing to individual methods. See the subsections.
#' 
#' @inheritSection Sim_Ling_1998 method
#' @inheritSection Sim_Resnik_1999 method
#' @inheritSection Sim_FaITH_2010 method
#' @inheritSection Sim_PS_2008 method
#' @inheritSection Sim_Relevance_2006 method
#' @inheritSection Sim_SimIC_2010 method
#' @inheritSection Sim_EISI_2015 method
#' @inheritSection Sim_XGraSM_2013 method
#' @inheritSection Sim_GraSM_2005 method
#' @inheritSection Sim_AIC_2014 method
#' @inheritSection Sim_Zhang_2006 method
#' @inheritSection Sim_universal method
#' @inheritSection Sim_Wang_2007 method
#' @inheritSection Sim_Rada_1989 method
#' @inheritSection Sim_Resnik_edge_2005 method
#' @inheritSection Sim_Leocock_1998 method
#' @inheritSection Sim_WP_1994 method
#' @inheritSection Sim_Slimani_2006 method
#' @inheritSection Sim_Shenoy_2012 method
#' @inheritSection Sim_Pekar_2002 method
#' @inheritSection Sim_Stojanovic_2001 method
#' @inheritSection Sim_Wang_edge_2012 method
#' @inheritSection Sim_Zhong_2002 method
#' @inheritSection Sim_AlMubaid_2006 method
#' @inheritSection Sim_Li_2003 method
#' @inheritSection Sim_RSS_2013 method
#' @inheritSection Sim_HRSS_2013 method
#' @inheritSection Sim_Shen_2010 method
#' @inheritSection Sim_SSDD_2013 method
#' @inheritSection Sim_Jiang_1997 method
#' @inheritSection Sim_Kappa method
#' @inheritSection Sim_Jaccard method
#' @inheritSection Sim_Dice method
#' @inheritSection Sim_Overlap method
#' 
#' @return A numeric symmetric matrix.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = 1:3,
#'     "b" = 3:4,
#'     "c" = 5,
#'     "d" = 7,
#'     "e" = 4:7,
#'     "f" = 8
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' term_sim(dag, dag_all_terms(dag), "Sim_Ling_1998")
term_sim = function(dag, terms, method, control = list()) {
	sim_fun = get_term_sim_method(method)

	sim_fun(dag, terms)
}

get_term_sim_method = function(method) {
	if(!method %in% ALL_TERM_SIM_METHODS) {
		stop("Supported similarity method should be in `ALL_TERM_SIM_METHODS`")
	}
	get(method, envir = topenv(), inherits = FALSE)
}

# cat(paste("#' @inheritSection",ALL_GROUP_SIM_METHODS, "method") , sep = "\n")

#' Semantic similarity between two groups of terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param group1 A vector of term names.
#' @param group2 A vector of term names.
#' @param method A group similarity method. All available methods are in [`ALL_GROUP_SIM_METHODS`].
#' @param sim_method A Term similarity method. All available methods are in [`ALL_TERM_SIM_METHODS`].
#' @param control A list of parameters passing to individual methods. See the subsections.
#'
#' @inheritSection GroupSim_pairwise_avg method
#' @inheritSection GroupSim_pairwise_max method
#' @inheritSection GroupSim_pairwise_BMM method
#' @inheritSection GroupSim_pairwise_ABM method
#' @inheritSection GroupSim_pairwise_HDF method
#' @inheritSection GroupSim_pairwise_VHDF method
#' @inheritSection GroupSim_pairwise_Froehlich_2007 method
#' @inheritSection GroupSim_pairwise_Joeng_2014 method
#' @inheritSection GroupSim_SimALN method
#' @inheritSection GroupSim_SimINT method
#' @inheritSection GroupSim_spgk method
#' @inheritSection GroupSim_SimGIC method
#' @inheritSection GroupSim_SimDIC method
#' @inheritSection GroupSim_SimUIC method
#' @inheritSection GroupSim_SimUI method
#' @inheritSection GroupSim_SimDB method
#' @inheritSection GroupSim_SimUB method
#' @inheritSection GroupSim_SimNTO method
#' @inheritSection GroupSim_SimCOU method
#' @inheritSection GroupSim_SimCOT method
#' @inheritSection GroupSim_SimLP method
#' @inheritSection GroupSim_Ye_2005 method
#' @inheritSection GroupSim_Cho_2007 method
#' @inheritSection GroupSim_SimALD method
#' @inheritSection GroupSim_Jaccard method
#' @inheritSection GroupSim_Dice method
#' @inheritSection GroupSim_Overlap method
#' @inheritSection GroupSim_Kappa method
#' 
#' @details
#' If `annotation` is set in `create_ontology_DAG()` and you want to directly calculate semantic similarity between two
#' annotated items, you can first get the associated terms of the two items by [`annotated_terms()`]:
#' 
#' ```
#' group1 = annotated_terms(dag, item1)
#' group2 = annotated_terms(dag, item2)
#' group_sim(dag, group1, group2, ...)
#' ```
#' 
#' @return A numeric scalar.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = 1:3,
#'     "b" = 3:4,
#'     "c" = 5,
#'     "d" = 7,
#'     "e" = 4:7,
#'     "f" = 8
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' group_sim(dag, c("c", "e"), c("d", "f"), 
#'     method = "GroupSim_pairwise_avg", 
#'     sim_method = "Sim_Ling_1998"
#' )
group_sim = function(dag, group1, group2, method, sim_method, control = list()) {
	group_sim_fun = get_group_sim_method(method)

	group_sim_fun(dag, group1, group2, sim_method)
}

get_group_sim_method = function(method) {
	if(!method %in% ALL_GROUP_SIM_METHODS) {
		stop("Supported similarity method should be in `ALL_GROUP_SIM_METHODS`")
	}
	get(method, envir = topenv(), inherits = FALSE)
}
