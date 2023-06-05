
library(testthat)


#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

annotation = list(
	"a" = 1:3,
	"b" = 3:4,
	"c" = 5,
	"d" = 7,
	"e" = 4:7,
	"f" = 8
)

dag = create_ontology_DAG(parents, children, annotation = annotation)

lca = LCA_term(dag, c("a", "b", "c", "d", "e", "f"))

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

lca = LCA_depth(dag, c("a", "b", "c", "d", "e", "f"))
lca_terms = LCA_term(dag, c("a", "b", "c", "d", "e", "f"))
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

mica = MICA_IC(dag, letters[1:6], IC_method = "IC_universal")
ic = IC_universal(dag)

test_that("test MICA_term", {
	for(i in 1:10) {
		nodes = sample(c("a", "b", "c", "d", "e", "f"), 2)
		expect_equal(
			mica[nodes[1], nodes[2]], 
			max(ic[dag_ancestor_from_two_groups(dag, nodes[1], nodes[2], "intersect", FALSE)])
		)
	}
})

