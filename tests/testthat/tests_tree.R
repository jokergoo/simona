


dag = create_ontology_DAG(c("a-h", "a-b", "a-c", "a-d", "b-e", "b-f", "c-g", "h-g", "d-e"))



lt_children2 = cpp_mark_tree_links(dag)
lt_children2 = lapply(lt_children2, function(x) {
	abs(x[x < 0])
})

n = length(lt_children2)
parents = rep(seq_len(n), times = vapply(lt_children2, length, FUN.VALUE = integer(1)))
children = unlist(lt_children2)

tree1 = create_ontology_DAG(dag@terms[parents], dag@terms[children])
	
tree2 = dag_treelize(dag)


test_that("test dag_treelize", {
	expect_identical(tree1@lt_children, tree2@lt_children)
	expect_identical(tree1@lt_parents, tree2@lt_parents)
	expect_identical(dag@terms, tree1@terms)
	expect_identical(dag@terms, tree2@terms)
})
