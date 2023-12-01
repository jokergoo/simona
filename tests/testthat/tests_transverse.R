
library(testthat)

## export all functions
if(!identical(topenv(), .GlobalEnv)) {
	pkg_env = asNamespace("simona")
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

test_that("test transverse", {
	for(letter in c("a", "b", "c", "d", "e", "f")) {
		expect_true(all(dag_parents(dag, letter) %in% dag_ancestors(dag, letter)))
		expect_true(all(dag_children(dag, letter) %in% dag_offspring(dag, letter)))
	}
})


test_that("test dag_siblings", {
	expect_equal(dag_siblings(dag, "b"), "c")
	expect_equal(dag_siblings(dag, "c"), c("b", "d"))
	expect_error(dag_siblings(dag, c("b", "c")), "only be length of one")
})

dag = create_ontology_DAG_from_GO_db()

all_go_id = dag@terms
test_that("test transverse on GO DAG", {
	for(term in sample(all_go_id, 10)) {
		expect_true(all(dag_parents(dag, term) %in% dag_ancestors(dag, term)))
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
		cpp_ancestors_within_background(dag, 5, 2:5),
		intersect(cpp_ancestors(dag, 5), 2:5)
	)
})

test_that("test ancestor/offspring of groups", {
	expect_equal(
		cpp_ancestors_of_two_groups(dag, 1, 2, 2),
		integer(0)
	)
	expect_equal(
		cpp_ancestors_of_two_groups(dag, 1, 2, 1),
		1
	)
	expect_equal(
		cpp_ancestors_of_a_group(dag, 1),
		integer(0)
	)
	expect_equal(
		cpp_ancestors_of_a_group(dag, 3),
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
		unname(dag_depth(dag)[dag@tpl_sorted]),
		unname(sort(dag_depth(dag)))
	)
})

test_that("test dag_depth", {
	expect_equal(
		unname(dag_depth(dag, use_cache = FALSE)),
		c(0, 1, 2, 2, 3, 3)
	)
})

test_that("test dag_height", {
	expect_equal(
		unname(dag_height(dag, use_cache = FALSE)),
		c(3, 2, 1, 1, 0, 0)
	)
})

test_that("test other distances", {
	expect_equal(
		unname(dag_shortest_dist_from_root(dag, use_cache = FALSE)),
		c(0, 1, 1, 2, 2, 3)
	)
	expect_equal(
		unname(dag_shortest_dist_to_leaves(dag, use_cache = FALSE)),
		c(2, 2, 1, 1, 0, 0)
	)

	expect_equal(
		unname(dag_longest_dist_to_offspring(dag, "b")),
		c(-1, 0, 1, 1, 2, 2)
	)

	expect_equal(
		dag_shortest_dist_to_offspring(dag, "a"),
		dag_shortest_dist_from_root(dag, use_cache = FALSE)
	)

	expect_equal(
		unname(dag_longest_dist_from_ancestors(dag, "e")),
		c(3, 2, 1, -1, 0, -1)
	)

	expect_equal(
		unname(dag_shortest_dist_from_ancestors(dag, "e")),
		c(2, 2, 1, -1, 0, -1)
	)
})

test_that("test n_children/n_parents/n_connected_leaves", {
	expect_equal(
		unname(n_offspring(dag)),
		c(5, 4, 1, 1, 0, 0)
	)
	expect_equal(
		unname(n_ancestors(dag)),
		c(0, 1, 2, 2, 3, 3)
	)
	expect_equal(
		unname(n_offspring(dag, include_self = TRUE)),
		c(5, 4, 1, 1, 0, 0) + 1
	)
	expect_equal(
		unname(n_ancestors(dag, include_self = TRUE)),
		c(0, 1, 2, 2, 3, 3) + 1
	)
	expect_equal(
		unname(n_connected_leaves(dag)),
		c(2, 2, 1, 1, 1, 1)
	)

	expect_equal(
		simona:::cpp_n_ancestors(dag_treelize(dag)),
		simona:::cpp_n_ancestors_on_tree(dag_treelize(dag))
	)

	expect_equal(
		simona:::cpp_n_offspring(dag_treelize(dag)),
		simona:::cpp_n_offspring_on_tree(dag_treelize(dag))
	)
})


tree = dag_treelize(dag)
test_that("test n_ancestors/n_offspring/n_connected_leaves on trees", {
	expect_equal(
		cpp_n_ancestors(tree),
		cpp_n_ancestors_on_tree(tree)
	)
	expect_equal(
		cpp_n_offspring(tree),
		cpp_n_offspring_on_tree(tree)
	)
	expect_equal(
		cpp_n_leaves(tree),
		cpp_n_leaves_on_tree(tree)
	)
})

test_that("test dag_distinct_ancestors", {
	expect_equal(
		dag_distinct_ancestors(dag, c("c", "d", "e", "f")), 
		c("c", "d")
	)
	expect_equal(
		dag_distinct_ancestors(dag, c("a", "b", "c", "d", "e", "f")), 
		c("a")
	)
})
