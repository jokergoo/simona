
#' All supported IC methods
#' @export
#' @return A vector all supported IC methods.
ALL_IC_METHODS = NULL


#' All supported similarity methods
#' @export
#' @return A vector all supported similarity methods.
ALL_TERM_SIM_METHODS = NULL

#' All supported group similarity methods
#' @export
#' @return A vector all supported grouped similarity methods.
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
