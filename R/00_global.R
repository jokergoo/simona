
#' All supported IC methods
ALL_IC_METHODS = NULL


#' All supported similarity methods
ALL_TERM_SIM_METHODS = NULL

#' All supported group similarity methods
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