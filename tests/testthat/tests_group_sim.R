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


group1 = c("c", "e", "d")
group2 = c("b", "d", "f")


for(method in all_group_sim_methods()) {
	cat(method, "\n")
	print(group_sim(dag, group1, group2, method = method))
}
