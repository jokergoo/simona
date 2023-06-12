

#' Information content
#' 
#' @param dag An `ontology_DAG` object.
#' @param method An IC method. All available methods are in [ALL_IC_METHODS].
#' @param control A list of parameters passing to individual methods.
#' 
#' @details
#' When using `IC_annotation`, `annotation` should be already set in [create_ontology_DAG()].
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

#' Semantic similarity
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param method A similarity method. All available methods are in [ALL_TERM_SIM_METHODS].
#' @param control A list of parameters passing to individual methods.
#' 
#' @return A numeric matrix.
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


#' Semantic similarity between two groups of terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param group1 A vector of term names.
#' @param group2 A vector of term names.
#' @param method Group similarity method. All available methods are in [ALL_GROUP_SIM_METHODS].
#' @param sim_method Term similarity method. All available methods are in [ALL_TERM_SIM_METHODS].
#' @param control A list of parameters passing to individual methods.
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
