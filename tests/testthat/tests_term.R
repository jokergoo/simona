
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


test_that("test IC_universal", {
	expect_equal(
		IC_universal(dag, F),
		-log(c(1, 1/2, 1/8, 1/4, 1/8, 1/4))
	)
})

test_that("test reachability", {
	expect_equal(
		reachability(dag, F),
		c(3, 2, 1, 1, 1, 1)
	)
})

test_that("test totipotency", {
	expect_equal(
		totipotency(dag, F),
		c(1, 5/6, 1/3, 1/3, 1/6, 1/6)
	)
})

test_that("test IC_Meng_2012", {
	expect_equal(
		IC_Meng_2012(dag, correct = FALSE, FALSE),
		c(0, 0, log(2)/log(3)*(1-log(4/3)/log(6)), log(2)/log(3)*(1-log(4/3)/log(6)), 1, 1)
	)
	expect_equal(
		IC_Meng_2012(dag, correct = TRUE, FALSE),
		c(0, log(1+1)/log(3+1)*(1-log(8/3)/log(6)), log(2+1)/log(3+1)*(1-log(4/3)/log(6)), log(2+1)/log(3+1)*(1-log(4/3)/log(6)), 1, 1)
	)
})
IC_Zhou_2008(dag, FALSE)
IC_Seco_2004(dag, FALSE)
IC_Zhang_2006(dag, FALSE)
IC_Seddiqui_2010(dag, FALSE)
IC_Sanchez_2011(dag, FALSE)

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

#####################
#   b--d--f
#  / \
# a---c--e
# upstream -> downstream

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")


dag = create_ontology_DAG(parents, children, relations = c("isa", "part of", "isa", "part of", "isa", "part of"), 
	annotation = annotation)
test_that("test IC_Wang_2007", {
	expect_equal(
		IC_Wang_2007(dag, c("isa" = 0.7, "part of" = 0.6)),
		c(1, 1.7, 2.3, 2.02, 2.61, 2.212)
	)
})


### test annotation
dag = create_ontology_DAG_from_GO_db("BP", org_db = "org.Hs.eg.db")
n = n_annotations(dag)
test_that("test n_annotations", {
	for(i in 1:10) {
		x = sample(dag@terms, 1)
		an = dag_ancestor(dag, x)
		expect_true(
			all(n[an] >= n[x])
		)
	}
})

if(FALSE) {

dag = create_ontology_DAG_from_GO_db("BP", org_db = "org.Hs.eg.db")
lt = lapply(ALL_IC_METHODS, function(method) {
	calc_IC(dag, method)
})
names(lt) = ALL_IC_METHODS

}
