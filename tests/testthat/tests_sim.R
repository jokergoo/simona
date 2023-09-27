
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

dag = create_ontology_DAG(parents, children, 
	relations = c("is_a", "part_of", "is_a", "part_of", "is_a", "part_of"), 
	annotation = annotation)


IC_annotation = -c(log(8/8), log(6/8), log(4/8), log(2/8), log(4/8), log(1/8))
IC_annotation[1] = 0

test_that("test Sim_AIC_2014", {
	sw = 1/(1 + exp(-1/IC_annotation))
	sim = term_sim(dag, dag@terms, "Sim_AIC_2014")

	for(i in 1:10) {
		t = sample(dag@terms, 2)
		a1 = dag_ancestors(dag, t[1], in_labels = FALSE, include_self = TRUE)
		a2 = dag_ancestors(dag, t[2], in_labels = FALSE, include_self = TRUE)

		expect_equal(sim[t[1], t[2]], 2*sum(sw[intersect(a1, a2)])/(sum(sw[a1]) + sum(sw[a2])))
	}

})

### Sw:

#   | a      b       c    d    e   f
# -------------------------------------
# a | 1
# b | 0.7    1
# c | 0.6    0.7     1
# d | 0.42   0.6     x    1
# e | 0.42   0.49  0.7    x   1 
# f | 0.252  0.36    x  0.6   x    1

test_that("test Sim_Wang_2007", {
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


test_that("test Sim_Wang_edge_2012", {
	m = term_sim(dag, dag@terms, "Sim_Wang_edge_2012")

	for(i in 1:10) {
		t = sample(dag@terms, 2)
		lca = LCA_term(dag, t)[1,2]

		x1 = dag_depth(dag, lca)
		if(lca == t[1]) {
			x2 = dag_depth(dag, lca)
		} else {
			x2 = dag_depth(dag, lca) + longest_distances_via_LCA(dag, c(lca, t[1]))[2,1]
		}
		if(lca == t[2]) {
			x3 = dag_depth(dag, lca)
		} else {
			x3 = dag_depth(dag, lca) + longest_distances_via_LCA(dag, c(lca, t[2]))[2,1]
		}

		if(x1 == 0) {
			expect_equal(m[t[1], t[2]], 0)
		} else {
			expect_equal(m[t[1], t[2]], unname(x1^2/x2/x3))
		}
	}
})

test_that("test Sim_Zhong_2002", {
	m = term_sim(dag, dag@terms, "Sim_Zhong_2002", control = list(depth_via_LCA = FALSE))

	d = dag_depth(dag)

	for(i in 1:10) {
		t = sample(dag@terms, 2)
		lca = LCA_term(dag, t)[1,2]

		expect_equal(m[t[1], t[2]], unname(1 - ( 2^(-d[lca]) - 2^(-(1+d[t[1]])) - 2^(-(1+d[t[2]])) )))
	}

	m = term_sim(dag, dag@terms, "Sim_Zhong_2002", control = list(depth_via_LCA = TRUE))

	for(i in 1:10) {
		t = sample(dag@terms, 2)
		lca = LCA_term(dag, t)[1,2]

		x = 2^(-d[lca])
		if(lca == t[1]) {
			y = 1 - 2^(-(1 + 0)) - 
			        2^(-(1+longest_distances_directed(dag, c(lca, t[2]))[lca, t[2]]))
		} else if(lca == t[2]) {
			y = 1 - 2^(-(1+longest_distances_directed(dag, c(lca, t[1]))[lca, t[1]])) - 
			        2^(-(1 + 0))
		} else {
			y = 1 - 2^(-(1+longest_distances_directed(dag, c(lca, t[1]))[lca, t[1]])) - 
			        2^(-(1+longest_distances_directed(dag, c(lca, t[2]))[lca, t[2]]))
		}
		expect_equal(m[t[1], t[2]], unname(1 - x*y))
	}
})

test_that("test Sim_Shen_2010", {
	m = term_sim(dag, dag@terms, "Sim_Shen_2010")

	for(i in 1:10) {
		t = sample(dag@n_terms, 2)
		lca = LCA_term(dag, t, in_labels = FALSE)[1,2]
		path1 = cpp_tpl_shortest_path(dag, lca, t[1])
		path2 = cpp_tpl_shortest_path(dag, lca, t[2])
		
		expect_equal(m[t[1], t[2]], 1 - atan(sum(1/IC_annotation[union(path1, path2)]))/pi*2, tolerance = 1e-5)
	}	
})


test_that("test Sim_SSDD_2013", {
	m = term_sim(dag, dag@terms, "Sim_SSDD_2013")
	tot = totipotency(dag)

	for(i in 1:10) {
		t = sample(dag@n_terms, 2)
		lca = LCA_term(dag, t, in_labels = FALSE)[1,2]
		path1 = cpp_tpl_shortest_path(dag, lca, t[1])
		path2 = cpp_tpl_shortest_path(dag, lca, t[2])
		path = union(path1, path2)
		
		expect_equal(m[t[1], t[2]], 1 - atan(sum(tot[path]))/pi*2, tolerance = 1e-5)
	}
})


