

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
#' - `all_ic_methods()`: A vector of all supported IC methods.
#' - `all_term_sim_methods()`: A vector of all supported term similarity methods.
#' - `all_group_sim_methods()`: A vector of all supported group similarity methods.
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
#' - `robot_jar`: Path of the `robot.jar` file. The file can be found from https://github.com/ontodev/robot/releases.
#' 
#' To set an option, you can use the `$`:
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
	robot_jar = NULL
)
