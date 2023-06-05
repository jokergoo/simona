
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

test_that("test topological sorting", {
	expect_equal(
		dag_depth(dag)[dag@tpl_sorted],
		sort(dag_depth(dag))
	)
})

test_that("test dag_depth", {
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		dag_depth_R(dag)
	)
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		c(0, 1, 2, 2, 3, 3)
	)
})

test_that("test dag_height", {
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		dag_height_R(dag)
	)
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

test_that("test IC_Meng_2012", {
	expect_equal(
		IC_Meng_2012(dag, FALSE),
		c(0, 0, log(2)/log(3)*(1-log(4/3)/log(6)), log(2)/log(3)*(1-log(4/3)/log(6)), 1, 1)
	)
	expect_equal(
		IC_Meng_2012(dag, FALSE, correct = TRUE),
		c(0, log(1+1)/log(3+1)*(1-log(8/3)/log(6)), log(2+1)/log(3+1)*(1-log(4/3)/log(6)), log(2+1)/log(3+1)*(1-log(4/3)/log(6)), 1, 1)
	)
})
IC_Zhou_2008(dag)
IC_Seco_2004(dag)
IC_Zhang_2006(dag)
IC_Seddiqui_2010(dag)
IC_Sanchez_2011(dag)

test_that("test IC_Wang_2007", {
	expect_error(
		IC_Wang_2007(dag),
		"not set"
	)
})



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

annotation2 = list(
	"a" = c(1, 2, 3, 4, 5, 7, 6, 8), 
	"b" = c(3, 4, 5, 7, 6, 8), 
	"c" = c(5, 4, 6, 7), 
	"d" = 7:8, 
	"e" = 4:7, 
	"f" = 8
)

dag1 = create_ontology_DAG(parents, children, annotation = annotation)
dag2 = create_ontology_DAG(parents, children, annotation = annotation2)

remove_attr = function(x) {
	attributes(x) = NULL
	x
}

test_that("test IC_annotation", {
	expect_equal(
		n_annotations(dag1),
		n_annotations(dag2)
	)
	expect_equal(
		remove_attr(IC_annotation(dag1)),
		-c(log(8/8), log(6/8), log(4/8), log(2/8), log(4/8), log(1/8))
	)
})

dag = create_ontology_DAG(parents, children, relations = rep("isa", length(parents)), 
	annotation = annotation)
IC_Wang_2007(dag)

##### test GOBP
dag = create_ontology_DAG_from_GO_db("BP")

test_that("test dag_depth", {
	expect_equal(
		dag_depth(dag, use_cache = FALSE),
		dag_depth_R(dag)
	)
})

test_that("test dag_height", {
	expect_equal(
		dag_height(dag, use_cache = FALSE),
		dag_height_R(dag)
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

invisible({
IC_Meng_2012(dag, F)
IC_Zhou_2008(dag, F)
IC_Seco_2004(dag, F)
IC_Zhang_2006(dag, F)
IC_Seddiqui_2010(dag, F)
IC_Sanchez_2011(dag, F)
IC_Wang_2007(dag, F)
})

if(FALSE) {

dag = create_ontology_DAG_from_GO_db("BP", org_db = "org.Hs.eg.db")
lt = lapply(ALL_IC_METHODS, function(method) {
	calc_IC(dag, method)
})
names(lt) = ALL_IC_METHODS

}
