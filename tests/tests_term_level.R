
library(testthat)


test_that("test cyclic_node", {
	parents = c("a", "b", "c", "d")
	children = c("b", "c", "d", "a")
	expect_error(
		createOntologyDAG(parents, children),
		"not a DAG"
	)

	parents = c("a", "b", "c", "g", "h")
	children = c("b", "c", "d", "h", "i")
	expect_warning(
		createOntologyDAG(parents, children, add_super_root = FALSE),
		"more than one root"
	)

	dag = createOntologyDAG(parents, children, add_super_root = TRUE)
	expect_equal(
		dag@n_terms,
		length(unique(c(parents, children))) + 1
	)
})



#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = createOntologyDAG(parents, children)

test_that("test dag_depth", {
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		dag_depth_R(dag, use_cache = FALSE)
	)
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		c(0, 1, 2, 2, 3, 3)
	)
})

test_that("test dag_height", {
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		dag_height_R(dag, use_cache = FALSE)
	)
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		c(3, 2, 1, 1, 0, 0)
	)
})

test_that("test IC_universal", {
	expect_equal(
		IC_universal_recursive(dag, F),
		IC_universal_bfs(dag, F)
	)
	expect_equal(
		IC_universal_recursive(dag, F),
		-log(c(1, 1/2, 1/8, 1/4, 1/8, 1/4))
	)
})

test_that("test reachability", {
	expect_equal(
		reachability_bfs(dag, F),
		reachability_recursive(dag, F)
	)
	expect_equal(
		reachability_bfs(dag, F),
		c(3, 2, 1, 1, 1, 1)
	)
})

test_that("test totipotency", {
	expect_equal(
		totipotency_bfs(dag, F),
		totipotency_recursive(dag, F)
	)
	expect_equal(
		totipotency_bfs(dag, F),
		c(1, 5/6, 1/3, 1/3, 1/6, 1/6)
	)
})

IC_Meng_2012(dag)
IC_Zhou_2008(dag)
IC_Seco_2004(dag)
IC_Zhang_2006(dag)
IC_Seddiqui_2010(dag)
IC_Sanchez_2011(dag)

test_that("test IC_Wang_2007", {
	expect_error(
		IC_Wang_2007(dag),
		"was not set"
	)
})

annotation = list(
	"a" = 1:3,
	"b" = 3:4,
	"c" = 5,
	"d" = 7,
	"e" = 4:7,
	"f" = 8
)

annotation2 = list(
	"a" = c(1, 2, 3, 4, 5, 7, 6, 8), 
	"b" = c(3, 4, 5, 7, 6, 8), 
	"c" = c(5, 4, 6, 7), 
	"d" = 7:8, 
	"e" = 4:7, 
	"f" = 8
)

test_that("test totipotency", {
	expect_equal(
		IC_annotation(dag, annotation, merge_to_ancestor = TRUE),
		IC_annotation(dag, annotation2, merge_to_ancestor = FALSE)
	)
})


##### test GOBP

df = toTable(GOBPCHILDREN)
dag = createOntologyDAG(parents = df[, 2], children = df[, 1], relations = df[, 3])

test_that("test dag_depth", {
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		dag_depth_R(dag, use_cache = FALSE)
	)
})

test_that("test dag_height", {
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		dag_height_R(dag, use_cache = FALSE)
	)
})

test_that("test IC_universal", {
	expect_equal(
		IC_universal_recursive(dag, F),
		IC_universal_bfs(dag, F)
	)
})

test_that("test reachability", {
	expect_equal(
		reachability_bfs(dag, F),
		reachability_recursive(dag, F)
	)
})

test_that("test totipotency", {
	expect_equal(
		totipotency_bfs(dag, F),
		totipotency_recursive(dag, F)
	)
})

IC_Meng_2012(dag)
IC_Zhou_2008(dag)
IC_Seco_2004(dag)
IC_Zhang_2006(dag)
IC_Seddiqui_2010(dag)
IC_Sanchez_2011(dag)
IC_Wang_2007(dag)

library(org.Hs.eg.db)
annotation = lapply(as.list(org.Hs.egGO2ALLEGS), unique)
test_that("test IC_annotation", {
	expect_equal(
		IC_annotation(dag, annotation),
		IC_annotation(dag, as.list(org.Hs.egGO2ALLEGS))
	)

	expect_equal(
		IC_annotation(dag, annotation),
		IC_annotation(dag, as.list(org.Hs.egGO2EG), merge_to_ancestor = TRUE)
	)
})

