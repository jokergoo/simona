
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

test_that("test tpl paths", {
	expect_equal(
		cpp_tpl_shortest_path_length(dag, 1, 3),
		1
	)
	expect_equal(
		cpp_tpl_shortest_path_length(dag, 1, 5),
		2
	)
	expect_equal(
		cpp_tpl_longest_path_length(dag, 1, 3),
		2
	)
	expect_equal(
		cpp_tpl_longest_path_length(dag, 1, 5),
		3
	)
	expect_equal(
		cpp_tpl_shortest_path_length(dag, 1, 4),
		cpp_tpl_longest_path_length(dag, 1, 4)
	)

	## path
	expect_equal(
		cpp_tpl_shortest_path(dag, 1, 3),
		c(1, 3)
	)
	expect_equal(
		cpp_tpl_shortest_path(dag, 1, 5),
		c(1, 3, 5)
	)
	expect_equal(
		cpp_tpl_longest_path(dag, 1, 3),
		c(1, 2, 3)
	)
	expect_equal(
		cpp_tpl_longest_path(dag, 1, 5),
		c(1, 2, 3, 5)
	)
	expect_equal(
		cpp_tpl_shortest_path(dag, 1, 4),
		cpp_tpl_longest_path(dag, 1, 4)
	)

	## test the other distance method
	m = cpp_longest_distances_directed(dag, 1:6)
	for(i in 1:6) {
		for(j in 1:6) {
			expect_equal(
				m[i, j],
				cpp_tpl_longest_path_length(dag, i, j)
			)
		}
	}
	
	m = cpp_shortest_distances_directed(dag, 1:6)
	for(i in 1:6) {
		for(j in 1:6) {
			expect_equal(
				m[i, j],
				cpp_tpl_shortest_path_length(dag, i, j)
			)
		}
	}
})

### test on GO BP

dag = create_ontology_DAG_from_GO_db()
depth = dag_depth(dag)

test_that("test two dist methods with GO BP", {
	for(i in 1:10) {
		go_id_1 = sample(dag@terms[depth > 5], 1)
		go_id_2 = sample(dag_ancestor(dag, go_id_1), 1)
		j = which(dag@terms == go_id_1)
		i = which(dag@terms == go_id_2)

		expect_equal(
			cpp_tpl_shortest_path_length(dag, i, j),
			cpp_shortest_distances_directed(dag, c(i, j))[1, 2]
		)

		expect_equal(
			cpp_tpl_longest_path_length(dag, i, j),
			cpp_longest_distances_directed(dag, c(i, j))[1, 2]
		)
	}
})


