
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


test_that("test IC_universal", {
	expect_equal(
		IC_universal(dag, use_cache = FALSE),
		-log(c(1, 1/2, 1/8, 1/4, 1/8, 1/4))
	)
})

test_that("test reachability", {
	expect_equal(
		reachability(dag, use_cache = FALSE),
		c(3, 2, 1, 1, 1, 1)
	)
})

test_that("test totipotency", {
	expect_equal(
		totipotency(dag, use_cache = FALSE),
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
IC_Zhou_2008(dag, use_cache = FALSE)
IC_Seco_2004(dag, use_cache = FALSE)
IC_Zhang_2006(dag, use_cache = FALSE)
# IC_Seddiqui_2010(dag, use_cache = FALSE)
IC_Sanchez_2011(dag, use_cache = FALSE)

test_that("test IC_Wang_2007", {
	expect_error(
		IC_Wang_2007(dag, use_cache = FALSE),
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
		n_annotations(dag1, use_cache = FALSE),
		n_annotations(dag2, use_cache = FALSE)
	)
	expect_equal(
		remove_attr(IC_annotation(dag1, use_cache = FALSE)),
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
		IC_Wang_2007(dag, c("is_a" = 0.7, "part of" = 0.6), use_cache = FALSE),
		c(1, 1.7, 2.3, 2.02, 2.61, 2.212)
	)
})

library(igraph)
g = dag_as_igraph(dag)
E(g)$weight = c("is_a" = 0.7, "part of" = 0.6)[E(g)$relation]
d = distances(g, mode = "out", weights = -log(E(g)$weight))
s = exp(-d)
test_that("test IC_Wang_2007 and shortest path weighted by 1/w", {
	expect_equal(
		IC_Wang_2007(dag, c("is_a" = 0.7, "part of" = 0.6), use_cache = FALSE),
		unname(colSums(s))
	)
})

### test annotation
dag = create_ontology_DAG_from_GO_db("BP", org_db = "org.Hs.eg.db")
n = n_annotations(dag)
test_that("test n_annotations", {
	for(i in 1:10) {
		x = sample(dag@terms, 1)
		an = dag_ancestors(dag, x)
		expect_true(
			all(n[an] >= n[x])
		)
	}
})

if(FALSE) {

dag = create_ontology_DAG_from_GO_db("BP", org_db = "org.Hs.eg.db")
lt = lapply(all_ic_methods(), function(method) {
	cat("=====", method, "=====\n")
	term_IC(dag, method)
})
names(lt) = all_ic_methods()

df = as.data.frame(lt)
pairs(df, pch = ".", col = dag_depth(dag))

}
