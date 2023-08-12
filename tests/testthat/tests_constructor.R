
library(testthat)

test_that("test cyclic_node", {
	parents = c("a", "b", "c", "d")
	children = c("b", "c", "d", "b")
	expect_error(
		create_ontology_DAG(parents, children),
		"find a cyclic"
	)

	parents = c("a", "b", "c", "g", "h")
	children = c("b", "c", "d", "h", "i")
	expect_message(
		dag <- create_ontology_DAG(parents, children),
		"more than one root"
	)
	expect_equal(
		length(dag@terms),
		length(unique(c(parents, children))) + 1
	)

	parents = c("a", "b", "c", "d")
	children = c("b", "c", "d", "a")
	expect_error(
		create_ontology_DAG(parents, children),
		"There might exist a cycle"
	)
})


#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children)

test_that("test sub-DAG", {
	expect_identical(
		dag[["c"]]@terms,
		c("c", "e")
	)
	expect_identical(
		dag[["b"]]@terms,
		c("b", "c", "d", "e", "f")
	)
	expect_error(
		dag["f"]
	)
})


test_that("test DAG filter", {
	expect_identical(
		dag_filter(dag, terms = c("b", "d", "f"))@terms,
		c("b", "d", "f")
	)
	expect_identical(
		dag_filter(dag, root = "b")@terms,
		c("b", "c", "d", "e", "f")
	)
	expect_identical(
		dag_filter(dag, root = c("b", "c"))@terms,
		c("b", "c", "d", "e", "f")
	)
	expect_identical(
		dag_filter(dag, leaves = c("c", "d"))@terms,
		c("a", "b", "c", "d")
	)
	expect_identical(
		dag_filter(dag, leaves = c("b", "c"))@terms,
		c("a", "b", "c")
	)
	expect_identical(
		dag_filter(dag, root = "b", leaves = "e")@terms,
		c("b", "c", "e")
	)
})


parents  = c("a", "b", "c", "d", "e")
children = c("b", "c", "d", "e", "b")

test_that("test cyclic path", {
	expect_error(
		create_ontology_DAG(parents, children), 
		"Found cyclic nodes"
	)

	expect_message(
		create_ontology_DAG(parents, children, remove_cyclic_paths = TRUE), 
		"Remove"
	)
})

parents  = c("a", "b", "c", "d", "f", "g", "h")
children = c("b", "c", "d", "e", "g", "h", "f")

test_that("test isolated rings", {
	expect_error(
		create_ontology_DAG(parents, children),
		"Found isolated rings"
	)

	expect_message(
		create_ontology_DAG(parents, children, remove_rings = TRUE), 
		"Remove"
	)
})
