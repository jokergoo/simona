
library(testthat)


#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = createOntologyDAG(parents, children)

lca = LCA_term(dag)
lca = structure(dag@terms[lca], dim = dim(lca), dimnames = list(dag@terms, dag@terms))

test_that("test LCA_term", {
	expect_equal(lca["e", "f"], "b")
	expect_equal(lca["e", "d"], "b")
	expect_equal(lca["c", "f"], "b")
	expect_equal(lca["c", "c"], "c")
	expect_equal(lca["b", "c"], "b")
	expect_equal(lca["d", "f"], "d")
	expect_equal(lca["b", "b"], "b")
	expect_equal(lca["d", "d"], "d")
	expect_equal(lca["a", "a"], "a")
	expect_equal(lca["a", "b"], "a")
})

lca_terms = LCA_term(dag);lca_terms = structure(dag@terms[lca_terms], dim = dim(lca_terms), dimnames = list(dag@terms, dag@terms))
lca = LCA_depth(dag); dimnames(lca) = list(dag@terms, dag@terms)
depth = dag_depth(dag); names(depth) = dag@terms

test_that("test LCA_depth", {
	expect_equal(depth[[ lca_terms["e", "f"] ]], lca["e", "f"])
	expect_equal(depth[[ lca_terms["e", "d"] ]], lca["e", "d"])
	expect_equal(depth[[ lca_terms["c", "f"] ]], lca["c", "f"])
	expect_equal(depth[[ lca_terms["c", "c"] ]], lca["c", "c"])
	expect_equal(depth[[ lca_terms["b", "c"] ]], lca["b", "c"])
	expect_equal(depth[[ lca_terms["d", "f"] ]], lca["d", "f"])
	expect_equal(depth[[ lca_terms["b", "b"] ]], lca["b", "b"])
	expect_equal(depth[[ lca_terms["d", "d"] ]], lca["d", "d"])
	expect_equal(depth[[ lca_terms["a", "a"] ]], lca["a", "a"])
	expect_equal(depth[[ lca_terms["a", "b"] ]], lca["a", "b"])
})

