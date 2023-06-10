
library(testthat)


#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children)

test_that("test transverse", {
	for(letter in c("a", "b", "c", "d", "e", "f")) {
		expect_true(all(dag_parents(dag, letter) %in% dag_ancestor(dag, letter)))
		expect_true(all(dag_children(dag, letter) %in% dag_offspring(dag, letter)))
	}
})


dag = create_ontology_DAG_from_GO_db()

all_go_id = dag@terms
test_that("test transverse on GO DAG", {
	for(term in sample(all_go_id, 10)) {
		expect_true(all(dag_parents(dag, term) %in% dag_ancestor(dag, term)))
		expect_true(all(dag_children(dag, term) %in% dag_offspring(dag, term)))
	}
})


parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children)

test_that("test ancestor/offspring within background", {
	expect_equal(
		cpp_offspring_within_background(dag, 2, 1:4),
		intersect(cpp_offspring(dag, 2), 1:4)
	)

	expect_equal(
		cpp_ancestor_within_background(dag, 5, 2:5),
		intersect(cpp_ancestor(dag, 5), 2:5)
	)
})

test_that("test ancestor/offspring of groups", {
	expect_equal(
		cpp_ancestor_of_two_groups(dag, 1, 2, 2),
		integer(0)
	)
	expect_equal(
		cpp_ancestor_of_two_groups(dag, 1, 2, 1),
		1
	)
	expect_equal(
		cpp_ancestor_of_a_group(dag, 1),
		integer(0)
	)
	expect_equal(
		cpp_ancestor_of_a_group(dag, 3),
		c(1, 2)
	)
	expect_equal(
		cpp_offspring_of_a_group(dag, 6),
		integer(0)
	)
	expect_equal(
		cpp_offspring_of_a_group(dag, 4),
		6
	)
})

test_that("test cpp_is_reachable", {
	lm = cpp_is_reachable(dag, 1:6)
	expect_false(lm[4, 5])
	expect_true(lm[2, 3])

	lm = cpp_is_reachable(dag, 1:6, TRUE)
	expect_true(lm[2, 3])
	expect_false(lm[3, 2])
})


test_that("test topological sorting", {
	expect_equal(
		dag_depth(dag)[dag@tpl_sorted],
		sort(dag_depth(dag))
	)
})

test_that("test dag_depth", {
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		c(0, 1, 2, 2, 3, 3)
	)
})

test_that("test dag_height", {
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		c(3, 2, 1, 1, 0, 0)
	)
})

test_that("test n_children/n_parents/n_leaves", {
	expect_equal(
		n_offspring(dag),
		c(5, 4, 1, 1, 0, 0)
	)
	expect_equal(
		n_ancestor(dag),
		c(0, 1, 2, 2, 3, 3)
	)
	expect_equal(
		n_leaves(dag),
		c(2, 2, 1, 1, 0, 0)
	)
})
