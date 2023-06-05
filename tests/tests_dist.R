
library(testthat)


## export all functions
if(!identical(topenv(), .GlobalEnv)) {
	pkg_env = asNamespace("ontsim")
	all_objs = ls(envir = pkg_env)
	for(obj in all_objs) {
		assign(obj, get(obj, envir = pkg_env, inherits = FALSE))
	}
}

#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children)

shortest_path_length(dag, from  = letters[1:6])
longest_path_length(dag, from  = letters[1:6])

shortest_path_length(dag, from  = letters[1:3], to = letters[3:6])

test_that("test finding paths", {
	expect_equal(
		shortest_path(dag, 1, 3),
		c(1, 3)
	)
	expect_equal(
		shortest_path(dag, 1, 3, in_labels = TRUE),
		c("a", "c")
	)
	expect_equal(
		longest_path(dag, 1, 3),
		c(1, 2, 3)
	)
	expect_equal(
		longest_path(dag, 1, 3, in_labels = TRUE),
		c("a", "b", "c")
	)
})


parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children, rep("isa", 6))

shortest_path_length(dag, from  = letters[1:6], weight = c("isa" = 0.6))
shortest_path_length(dag, from  = letters[1:6], weight = 0.6)

dag = create_ontology_DAG(parents, children, c("small", "normal", "small", "normal", "normal", "normal"))
shortest_path_length(dag, from  = letters[1:6], weight = c("small" = 0.2, "normal" = 1))
shortest_path(dag, "a", "c", weight = c("small" = 0.2, "normal" = 1))