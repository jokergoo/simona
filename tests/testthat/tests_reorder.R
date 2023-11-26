library(testthat)




parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")

dag = create_ontology_DAG(parents, children)

test_that("test dag_reorder", {
	dag2 = dag_reorder(dag, value = c(1, 1, 10, 1, 10, 1))
	expect_equal(dag2@lt_children[[2]], rev(dag@lt_children[[2]]))

	dag3 = dag_reorder(dag, value = c(10, 1))
	expect_equal(dag3@lt_children[[2]], rev(dag@lt_children[[2]]))
})




dag = create_ontology_DAG(c("a-h", "a-b", "a-c", "a-d", "b-e", "b-f", "c-g", "h-g", "d-e"))
tree = dag_treelize(dag)
lt = cpp_get_force_counterpart(dag@lt_children, dag@lt_parents, tree@lt_children, tree@lt_parents, dag@root)

test_that("test cpp_get_force_counterpart", {
	expect_equal(lt[[1]], integer(0))
	expect_equal(lt[[2]], 4)
	expect_equal(lt[[3]], 8)
	expect_equal(lt[[4]], 5)
	expect_equal(lt[[5]], 4)
	expect_equal(lt[[6]], integer(0))
	expect_equal(lt[[7]], 8)
	expect_equal(lt[[8]], 7)
})

test_that("test move_index", {
	x = c(2, 1, 5, 4, 3)
	od = order(-abs(x))
	expect_equal(x[move_index(x, od-1, 1) + 1], c(5, 2, 1, 4, 3))
	expect_equal(x[move_index(x, od-1, 2) + 1], c(5, 4, 2, 1, 3))
	expect_equal(x[move_index(x, od-1, 3) + 1], c(5, 4, 3, 2, 1))
	expect_equal(x[move_index(x, od-1, 4) + 1], c(5, 4, 3, 2, 1))
	expect_equal(x[move_index(x, od-1, 5) + 1], c(5, 4, 3, 2, 1))

	x = c(-2, -1, -5, -4, -3)
	od = order(-abs(x))
	expect_equal(x[move_index(x, od-1, 1) + 1], c(-2, -1, -4, -3, -5))
	expect_equal(x[move_index(x, od-1, 2) + 1], c(-2, -1, -3, -4, -5))
	expect_equal(x[move_index(x, od-1, 3) + 1], c(-2, -1, -3, -4, -5))
	expect_equal(x[move_index(x, od-1, 4) + 1], c(-1, -2, -3, -4, -5))
	expect_equal(x[move_index(x, od-1, 5) + 1], c(-1, -2, -3, -4, -5))

	x = c(-2, 1, 5, -4, 3)
	od = order(-abs(x))
	expect_equal(x[move_index(x, od-1, 1) + 1], c(5, -2, 1, -4, 3))
	expect_equal(x[move_index(x, od-1, 2) + 1], c(5, -2, 1, 3, -4))
	expect_equal(x[move_index(x, od-1, 3) + 1], c(5, 3, -2, 1, -4))
	expect_equal(x[move_index(x, od-1, 4) + 1], c(5, 3, 1, -2, -4))
	expect_equal(x[move_index(x, od-1, 5) + 1], c(5, 3, 1, -2, -4))

	x = c(3, 1, 5, -4, -2)
	od = order(-abs(x))
	expect_equal(x[move_index(x, od-1, 1, FALSE) + 1], c(3, 1, -4, -2, 5))
	expect_equal(x[move_index(x, od-1, 2, FALSE) + 1], c(-4, 3, 1, -2, 5))
	expect_equal(x[move_index(x, od-1, 3, FALSE) + 1], c(-4, 1, -2, 3, 5))
	expect_equal(x[move_index(x, od-1, 4, FALSE) + 1], c(-4, -2, 1, 3, 5))
	expect_equal(x[move_index(x, od-1, 5, FALSE) + 1], c(-4, -2, 1, 3, 5))
})


test_that("test calc_x_offset", {
	prev_od = 1:5
	new_od = 1:5
	expect_equal(calc_x_offset(1:5, prev_od - 1, new_od - 1, 1:5),
		         c(0, 0, 0, 0, 0))

	prev_od = 1:5
	new_od = c(2, 1, 3, 4, 5)
	expect_equal(calc_x_offset(1:5, prev_od - 1, new_od - 1, 1:5),
		         c(2, -1, 0, 0, 0))

	prev_od = 1:5
	new_od = 5:1
	expect_equal(calc_x_offset(1:5, prev_od - 1, new_od - 1, 1:5),
		         c(14, 11, 6, -1, -10))
})


# pos = cpp_node_pos_in_tree(tree, n_connected_leaves(tree))
# lt_counterpart = cpp_get_force_counterpart(dag@lt_children, dag@lt_parents, tree@lt_children, tree@lt_parents, dag@root)


# force = cpp_get_force(lt_counterpart, pos$x, dag_depth(tree))
# test_that("test cpp_get_force", {
# 	expect_equal(sign(force), c(0, 1, -1, -1, 1, 0, -1, 1))
# })


# n_cp = sapply(lt_counterpart, length)
# x = pos$x
# test_that("test reorder_children", {
# 	expect_equal(reorder_children(tree@lt_children[[1]], n_cp, force, pos$width, dag_depth(tree), x, tree@lt_children),
# 		         c(3, 4, 2, 8))
# 	expect_equal(order(x[c(3, 4, 2, 8)]), 1:4)
# })


# pos = cpp_node_pos_in_tree(tree, n_connected_leaves(tree))
# cpp_reorder_tree_x(tree, lt_counterpart, pos$x, pos$width)


