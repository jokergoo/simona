
DEFAULT_RELATIONS = c("isa" = "is_a", 
	                  "is a" = "is_a", 
	                  "is_a" = "is_a", 
	                  "is-a" = "is_a", 
	                  "part_a" = "part_of", 
	                  "part a" = "part_of", 
	                  "part-a" = "part_of", 
	                  "part of" = "part_of", 
	                  "part_of" = "part_of",
	                  "part-of" = "part_of"
	                  )


.ALL_IC_METHODS = NULL
.ALL_TERM_SIM_METHODS = NULL
.ALL_GROUP_SIM_METHODS = NULL


ADD_IC_METHOD = function(method, param = character(0)) {
	.ALL_IC_METHODS[[method]] <<- param
}

ADD_TERM_SIM_METHOD = function(method, param = character(0)) {
	.ALL_TERM_SIM_METHODS[[method]] <<- param
}

ADD_GROUP_SIM_METHOD = function(method, param = character(0)) {
	.ALL_GROUP_SIM_METHODS[[method]] <<- param
}


#' Supported methods
#' 
#' @details
#' - `all_ic_methods()`: A vector all supported IC methods.
#' - `all_term_sim_methods()`: A vector all supported term similarity methods.
#' - `all_group_sim_methods()`: A vector all supported group similarity methods.
#' 
#' @rdname all_methods
#' @export
#' @return A vector of all supported methods.
#' @examples
#' all_ic_methods()
#' all_term_sim_methods()
#' all_group_sim_methods()
all_ic_methods = function() {
	names(.ALL_IC_METHODS)
}

#' @export
#' @rdname all_methods
all_term_sim_methods = function() {
	names(.ALL_TERM_SIM_METHODS)
}

#' @export
#' @rdname all_methods
all_group_sim_methods = function() {
	names(.ALL_GROUP_SIM_METHODS)
}


get_IC_method = function(method, control = list()) {
	param = .ALL_IC_METHODS[[method]]
	if(is.null(param)) {
		method = grep(method, names(.ALL_IC_METHODS), ignore.case = TRUE, value = TRUE)
		if(length(method) != 1) {
			stop("method should be in `all_ic_methods()`.")
		}
		param = .ALL_IC_METHODS[[method]]
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
		method = grep(method, names(.ALL_TERM_SIM_METHODS), ignore.case = TRUE, value = TRUE)
		if(length(method) != 1) {
			stop("method should be in `all_term_sim_methods()`.")
		}
		param = .ALL_TERM_SIM_METHODS[[method]]
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
		method = grep(method, names(.ALL_GROUP_SIM_METHODS), ignore.case = TRUE, value= TRUE)
		if(length(method) != 1) {
			stop("method should be in `all_group_sim_methods()`.")
		}
		param = .ALL_GROUP_SIM_METHODS[[method]]
	}

	f = get(method, envir = topenv(), inherits = FALSE)

	control = control[intersect(names(control), param)]
	
	function(dag, group1, group2) {
		do.call(f, c(list(dag = dag, group1 = group1, group2 = group2), control))
	}
}


#' Global options
#' 
#' @param RESET Ignore.
#' @param READ.ONLY Ignore.
#' @param LOCAL Ignore.
#' @param ADD Ignore.
#' @param ... Ignore.
#' 
#' @export
#' @import GlobalOptions
#' @examples
#' simone_opt
simone_opt = setGlobalOptions(
	use_cache = list(
		.value = TRUE,
		.class = "logical",
		.length = 1
	),
	verbose = list(
		.value = TRUE,
		.class = "logical",
		.length = 1
	),
	anno_uniquify = list(
		.value = TRUE,
		.class = "logical",
		.length = 1
	)
)
