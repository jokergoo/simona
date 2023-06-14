
#' Supported methods
#' 
#' @details
#' - `ALL_IC_METHODS`: A vector all supported IC methods.
#' - `ALL_TERM_SIM_METHODS`: A vector all supported term similarity methods.
#' - `ALL_GROUP_SIM_METHODS`: A vector all supported group similarity methods.
#' 
#' @rdname all_methods
#' @export
#' @return A vector all supported methods.
#' @examples
#' ALL_IC_METHODS
#' ALL_TERM_SIM_METHODS
#' ALL_GROUP_SIM_METHODS
ALL_IC_METHODS = NULL


#' @export
#' @rdname all_methods
ALL_TERM_SIM_METHODS = NULL

#' @export
#' @rdname all_methods
ALL_GROUP_SIM_METHODS = NULL

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


ADD_IC_METHOD = function(method) {
	ALL_IC_METHODS <<- unique(c(ALL_IC_METHODS, method))
}

ADD_TERM_SIM_METHOD = function(method) {
	ALL_TERM_SIM_METHODS <<- unique(c(ALL_TERM_SIM_METHODS, method))
}

ADD_GROUP_SIM_METHOD = function(method) {
	ALL_GROUP_SIM_METHODS <<- unique(c(ALL_GROUP_SIM_METHODS, method))
}
