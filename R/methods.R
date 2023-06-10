

#' Information content
#' 
#' @param dag A `ontology_DAG` object.
#' @param IC_method An IC method. All available methods are in [ALL_IC_METHODS].
#' @param control A list of parameters passing to individual methods
#' 
#' @details
#' When using `IC_annotation`, `annotation` should be already set in [create_ontology_DAG()].
#' @return A numeric vector.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' for(method in ALL_IC_METHODS) {
#'     cat("===", method, "====\n")
#'     print(term_IC(dag, method)
#'     cat("\n")
#' }
term_IC = function(dag, IC_method, control = list()) {
	IC_fun = get_IC_method(IC_method)
	if(IC_method == "IC_Wang_2007") {
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
#' @param dag A `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param sim_method An IC method. All available methods are in [ALL_TERM_SIM_METHODS].
#' 
#' @return A numeric matrix.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' for(method in ALL_TERM_SIM_METHODS) {
#'     cat("===", method, "====\n")
#'     print(term_sim(dag, dag_all_terms(dag), method)
#'     cat("\n")
#' }
term_sim = function(dag, terms, sim_method) {
	sim_fun = get_term_sim_method(sim_method)

	sim_fun(dag, terms)
}

get_term_sim_method = function(method) {
	if(!method %in% ALL_TERM_SIM_METHODS) {
		stop("Supported similarity method should be in `ALL_TERM_SIM_METHODS`")
	}
	get(method, envir = topenv(), inherits = FALSE)
}


group_sim = function(dag, terms, group_sim_method) {
	sim_fun = get_group_sim_method(group_sim_method)

	sim_fun(dag, terms)
}

get_group_sim_method = function(method) {
	if(!method %in% ALL_GROUP_SIM_METHODS) {
		stop("Supported similarity method should be in `ALL_GROUP_SIM_METHODS`")
	}
	get(method, envir = topenv(), inherits = FALSE)
}
