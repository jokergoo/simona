


get_IC_method = function(method, control = list()) {
	param = .ALL_IC_METHODS[[method]]
	if(is.null(param)) {
		i = grep(method, names(.ALL_IC_METHODS), ignore.case = TRUE, value = TRUE)
		if(length(i) != 1) {
			stop("method should be in `all_ic_methods()`.")
		}
		param = .ALL_IC_METHODS[[i]]
	}

	f = get(method, envir = topenv(), inherits = FALSE)

	control = control[intersect(names(control), param)]
	
	function(dag) {
		do.call(f, c(list(dag = dag), control))
	}
}


get_term_sim_method = function(method, control = list()) {
	param = .ALL_TERM_SIM_METHODS[[method]]
	if(is.null(param)) {
		i = grep(method, names(.ALL_TERM_SIM_METHODS), ignore.case = TRUE, value = TRUE)
		if(length(i) != 1) {
			stop("method should be in `all_term_sim_methods()`.")
		}
		param = .ALL_TERM_SIM_METHODS[[i]]
	}

	f = get(method, envir = topenv(), inherits = FALSE)

	control = control[intersect(names(control), param)]
	
	function(dag, terms) {
		do.call(f, c(list(dag = dag, terms = terms), control))
	}
}


get_group_sim_method = function(method, control = list()) {
	param = .ALL_GROUP_SIM_METHODS[[method]]
	if(is.null(param)) {
		i = grep(method, names(.ALL_GROUP_SIM_METHODS), ignore.case = TRUE, value= TRUE)
		if(length(i) != 1) {
			stop("method should be in `all_group_sim_methods()`.")
		}
		param = .ALL_GROUP_SIM_METHODS[[i]]
	}

	f = get(method, envir = topenv(), inherits = FALSE)

	control = control[intersect(names(control), param)]
	
	function(dag, group1, group2) {
		do.call(f, c(list(dag = dag, group1 = group1, group2 = group2), control))
	}
}



# cat(paste("#' @inheritSection", all_ic_methods(), "Methods") , sep = "\n")

#' Information content
#' 
#' @param dag An `ontology_DAG` object.
#' @param method An IC method. All available methods are in [`all_ic_methods()`].
#' @param terms A vector of term names. If it is set, the returned vector will be subsetted to the terms that have been set here.
#' @param control A list of parameters passing to individual methods. See the subsections.
#' 
#' @inheritSection IC_offspring Methods
#' @inheritSection IC_depth Methods
#' @inheritSection IC_annotation Methods
#' @inheritSection IC_universal Methods
#' @inheritSection IC_Zhang_2006 Methods
#' @inheritSection IC_Seco_2004 Methods
#' @inheritSection IC_Zhou_2008 Methods
#' @inheritSection IC_Sanchez_2011 Methods
#' @inheritSection IC_Meng_2012 Methods
#' @inheritSection IC_Wang_2007 Methods
#' 
#' @return A numeric vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' term_IC(dag, "IC_annotation")
term_IC = function(dag, method, terms = NULL, control = list()) {
	IC_fun = get_IC_method(method, control)
	ic = IC_fun(dag)

	ic = structure(ic, names = dag@terms)

	if(!is.null(terms)) {
		i = term_to_node_id(dag, terms, strict = FALSE)
		ic[i]
	} else {
		ic
	}
}


# cat(paste("#' @inheritSection", all_term_sim_methods(), "Methods") , sep = "\n")

#' Semantic similarity
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param method A term similarity method. All available methods are in [`all_term_sim_methods()`].
#' @param control A list of parameters passing to individual methods. See the subsections.
#' 
#' @inheritSection Sim_Lin_1998 Methods
#' @inheritSection Sim_Resnik_1999 Methods
#' @inheritSection Sim_FaITH_2010 Methods
#' @inheritSection Sim_Relevance_2006 Methods
#' @inheritSection Sim_SimIC_2010 Methods
#' @inheritSection Sim_XGraSM_2013 Methods
#' @inheritSection Sim_GraSM_2005 Methods
#' @inheritSection Sim_EISI_2015 Methods
#' @inheritSection Sim_AIC_2014 Methods
#' @inheritSection Sim_Zhang_2006 Methods
#' @inheritSection Sim_universal Methods
#' @inheritSection Sim_Wang_2007 Methods
#' @inheritSection Sim_Rada_1989 Methods
#' @inheritSection Sim_Resnik_edge_2005 Methods
#' @inheritSection Sim_Leocock_1998 Methods
#' @inheritSection Sim_WP_1994 Methods
#' @inheritSection Sim_Slimani_2006 Methods
#' @inheritSection Sim_Shenoy_2012 Methods
#' @inheritSection Sim_Pekar_2002 Methods
#' @inheritSection Sim_Stojanovic_2001 Methods
#' @inheritSection Sim_Wang_edge_2012 Methods
#' @inheritSection Sim_Zhong_2002 Methods
#' @inheritSection Sim_AlMubaid_2006 Methods
#' @inheritSection Sim_Li_2003 Methods
#' @inheritSection Sim_RSS_2013 Methods
#' @inheritSection Sim_HRSS_2013 Methods
#' @inheritSection Sim_Shen_2010 Methods
#' @inheritSection Sim_SSDD_2013 Methods
#' @inheritSection Sim_Jiang_1997 Methods
#' @inheritSection Sim_Kappa Methods
#' @inheritSection Sim_Jaccard Methods
#' @inheritSection Sim_Dice Methods
#' @inheritSection Sim_Overlap Methods
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
#' term_sim(dag, dag_all_terms(dag), method = "Sim_Lin_1998")
term_sim = function(dag, terms, method, control = list()) {
	sim_fun = get_term_sim_method(method, control)

	sim_fun(dag, terms)
}


# cat(paste("#' @inheritSection", all_group_sim_methods(), "Methods") , sep = "\n")

#' Semantic similarity between two groups of terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param group1 A vector of term names.
#' @param group2 A vector of term names.
#' @param method A group similarity method. All available methods are in [`all_group_sim_methods()`].
#' @param control A list of parameters passing to individual methods. See the subsections.
#'
#' 
#' @section Methods:
#' 
#' @inheritSection GroupSim_pairwise_avg Methods
#' @inheritSection GroupSim_pairwise_max Methods
#' @inheritSection GroupSim_pairwise_BMA Methods
#' @inheritSection GroupSim_pairwise_BMM Methods
#' @inheritSection GroupSim_pairwise_ABM Methods
#' @inheritSection GroupSim_pairwise_HDF Methods
#' @inheritSection GroupSim_pairwise_MHDF Methods
#' @inheritSection GroupSim_pairwise_VHDF Methods
#' @inheritSection GroupSim_pairwise_Froehlich_2007 Methods
#' @inheritSection GroupSim_pairwise_Joeng_2014 Methods
#' @inheritSection GroupSim_SimALN Methods
#' @inheritSection GroupSim_SimINT Methods
#' @inheritSection GroupSim_spgk Methods
#' @inheritSection GroupSim_SimGIC Methods
#' @inheritSection GroupSim_SimDIC Methods
#' @inheritSection GroupSim_SimUIC Methods
#' @inheritSection GroupSim_SimUI Methods
#' @inheritSection GroupSim_SimDB Methods
#' @inheritSection GroupSim_SimUB Methods
#' @inheritSection GroupSim_SimNTO Methods
#' @inheritSection GroupSim_SimCOU Methods
#' @inheritSection GroupSim_SimCOT Methods
#' @inheritSection GroupSim_SimLP Methods
#' @inheritSection GroupSim_Ye_2005 Methods
#' @inheritSection GroupSim_SimCHO Methods
#' @inheritSection GroupSim_SimALD Methods
#' @inheritSection GroupSim_Jaccard Methods
#' @inheritSection GroupSim_Dice Methods
#' @inheritSection GroupSim_Overlap Methods
#' @inheritSection GroupSim_Kappa Methods
#' 
#' @details
#' If `annotation` is set in `create_ontology_DAG()` and you want to directly calculate semantic similarity between two
#' annotated items, you can first get the associated terms of the two items by [`annotated_terms()`]:
#' 
#' ```
#' group1 = annotated_terms(dag, item1)[[1]]
#' group2 = annotated_terms(dag, item2)[[1]]
#' group_sim(dag, group1, group2, ...)
#' ```
#' 
#' @return A numeric scalar.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' group_sim(dag, c("c", "e"), c("d", "f"), 
#'     method = "GroupSim_pairwise_avg", 
#'     control = list(term_sim_method = "Sim_Lin_1998")
#' )
group_sim = function(dag, group1, group2, method, control = list()) {
	group_sim_fun = get_group_sim_method(method, control)

	group_sim_fun(dag, group1, group2)
}
