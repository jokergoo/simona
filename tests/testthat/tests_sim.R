
library(testthat)


#### test a small dag

#   b--d--f
#  / \
# a---c--e
# upstream -> downstream
annotation = list(
	"a" = 1:3,
	"b" = 3:4,
	"c" = 5,
	"d" = 7,
	"e" = 4:7,
	"f" = 8
)

parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children, relations = c("is_a", "part_of", "is_a", "part_of", "is_a", "part_of"), 
	annotation = annotation)

# IC_annotation:
# -c(log(8/8), log(6/8), log(4/8), log(2/8), log(4/8), log(1/8))
# test_that("test sim_XGraSM", {
# 	m = 
# })

# test_that("test sim_EISI", {

# })

# test_that("test sim_AIC", {

# })

# test_that("test sim_zhong", {

# })


### Sw:

#   | a      b       c    d    e   f
# -------------------------------------
# a | 1
# b | 0.7    1
# c | 0.6    0.7     1
# d | 0.42   0.6     x    1
# e | 0.42   0.49  0.7    x   1 
# f | 0.252  0.36    x  0.6   x    1

test_that("test sim_wang", {
	m = cpp_sim_wang(dag, 1:6, c(0.7, 0.6))
	expect_equal(
		m[, 1],
		c(1, 1.7/2.7, 1.6/3.3, 1.42/3.02, 1.42/3.61, 1.252/3.212)
	)
	expect_equal(
		m[5, 6],
		(0.36+0.49+0.252+0.42)/(2.61+2.212)
	)

	m2 = Sim_Wang_2007(dag, letters[1:6], contribution_factor = c("is_a" = 0.7, "part_of" = 0.6))

	expect_equal(as.vector(m), as.vector(m2))
})
