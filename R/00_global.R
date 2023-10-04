

.ALL_IC_METHODS = NULL
.ALL_TERM_SIM_METHODS = NULL
.ALL_GROUP_SIM_METHODS = NULL


ADD_IC_METHOD = function(method, require_anno = FALSE) {
	fun = get(method, envir = topenv(), inherits = FALSE)
	param = names(formals(fun)[-1])
	attr(param, "require_anno") = require_anno
	.ALL_IC_METHODS[[method]] <<- param
}

ADD_TERM_SIM_METHOD = function(method, require_anno = FALSE) {
	fun = get(method, envir = topenv(), inherits = FALSE)
	param = names(formals(fun)[-(1:2)])
	attr(param, "require_anno") = require_anno
	.ALL_TERM_SIM_METHODS[[method]] <<- param
}

ADD_GROUP_SIM_METHOD = function(method, require_anno = FALSE) {
	fun = get(method, envir = topenv(), inherits = FALSE)
	param = names(formals(fun)[-(1:3)])
	attr(param, "require_anno") = require_anno
	.ALL_GROUP_SIM_METHODS[[method]] <<- param
}


#' Supported methods
#' 
#' @param require_anno If it is set to `TRUE`, methods that require external annotations are only returned. If 
#'    it is set to `FALSE`, methods that do not require annotations are returned. A value of `NULL` returns both.
#' 
#' @details
#' - `all_term_IC_methods()`: A vector of all supported IC methods.
#' - `all_term_sim_methods()`: A vector of all supported term similarity methods.
#' - `all_group_sim_methods()`: A vector of all supported group similarity methods.
#' 
#' @rdname all_methods
#' @export
#' @return A character vector of all supported methods.
#' @examples
#' all_term_IC_methods()
#' all_term_sim_methods()
#' all_group_sim_methods()
all_term_IC_methods = function(require_anno = NULL) {
	if(is.null(require_anno)) {
		names(.ALL_IC_METHODS)
	} else if(require_anno) {
		names(.ALL_IC_METHODS)[vapply(.ALL_IC_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	} else {
		names(.ALL_IC_METHODS)[!vapply(.ALL_IC_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	}
}

#' @export
#' @rdname all_methods
all_term_sim_methods = function(require_anno = NULL) {
	if(is.null(require_anno)) {
		names(.ALL_TERM_SIM_METHODS)
	} else if(require_anno) {
		names(.ALL_TERM_SIM_METHODS)[vapply(.ALL_TERM_SIM_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	} else {
		names(.ALL_TERM_SIM_METHODS)[!vapply(.ALL_TERM_SIM_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	}
}

#' @export
#' @rdname all_methods
all_group_sim_methods = function(require_anno = NULL) {
	if(is.null(require_anno)) {
		names(.ALL_GROUP_SIM_METHODS)
	} else if(require_anno) {
		names(.ALL_GROUP_SIM_METHODS)[vapply(.ALL_GROUP_SIM_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	} else {
		names(.ALL_GROUP_SIM_METHODS)[!vapply(.ALL_GROUP_SIM_METHODS, function(x) attr(x, "use_anno"), FUN.VALUE = character(1))]
	}
}

#' All Papameters of a given method
#' 
#' @param IC_method A single IC method name.
#' @param term_sim_method A single term similarity method name.
#' @param group_sim_method A single group similarity method name.
#' 
#' @return A vector of parameter names.
#' @export
#' @examples
#' method_param(IC_method = "IC_annotation")
#' method_param(term_sim_method = "Sim_Wang_2007")
method_param = function(IC_method = NULL, term_sim_method = NULL, group_sim_method = NULL) {
	if(!is.null(IC_method)) {
		p = .ALL_IC_METHODS[[IC_method]]
	} else if(!is.null(term_sim_method)) {
		p = .ALL_TERM_SIM_METHODS[[term_sim_method]]
	} else if(!is.null(group_sim_method)) {
		p = .ALL_GROUP_SIM_METHODS[[group_sim_method]]
	} else {
		return(NULL)
	}

	attributes(p) = NULL
	p
}


#' Global options
#' 
#' @param RESET Reset to default option values.
#' @param READ.ONLY Only return read only options.
#' @param LOCAL Only return local options.
#' @param ADD Add new options.
#' @param ... Name-value pairs for options.
#' 
#' @details
#' There are the following global options:
#' 
#' - `use_cache`: By default, information content of all terms is cached and reused. If `use_cache` is set to `FALSE`, IC will be re-calculated.
#' - `verbose`: Whether to print messages?
#' - `anno_uniquify`: In the annotation-based IC method, the union of items annotated to the term as well as all its offspring terms is used, which means
#'      the set of annotated items for the term is uniquified. If `anno_uniquify` is set to `FALSE`, the uniquification is not applied, we simply add the number
#'      of items annotated to the term and the numbers of items annotated to each of its offspring terms.
#' - `robot_jar`: Path of the `robot.jar` file. The file can be found from \url{https://github.com/ontodev/robot/releases}.
#' 
#' To set an option, you can use `$`:
#' 
#' ```
#' simona_opt$verbose = FALSE
#' ```
#' 
#' or use it as a function:
#' 
#' ```
#' simona_opt(verbose = FALSE)
#' ```
#' 
#' @export
#' @import GlobalOptions
#' @returns A single option value.
#' @examples
#' simona_opt
simona_opt = GlobalOptions::setGlobalOptions(
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
	),
	verbose = list(
		.value = TRUE,
		.class = "logical",
		.length = 1
	),
	robot_jar = NULL
)
