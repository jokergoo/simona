## ----echo = FALSE, message = FALSE--------------------------------------------
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.width = 6, fig.height = 6,
    fig.align = "center")

## -----------------------------------------------------------------------------
library(simona)
parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")
dag = create_ontology_DAG(parents, children)

## -----------------------------------------------------------------------------
dag

## -----------------------------------------------------------------------------
dag_root(dag)
dag_leaves(dag)
dag_is_leaf(dag, letters[1:6])

## -----------------------------------------------------------------------------
dag_all_terms(dag)
dag_n_terms(dag)

## -----------------------------------------------------------------------------
n_children(dag)
n_parents(dag)
n_offspring(dag)
n_ancestors(dag)

## -----------------------------------------------------------------------------
n_connected_leaves(dag)

## -----------------------------------------------------------------------------
dag_parents(dag, "c")
dag_parents(dag, c("d", "e")) # union of d's parents and e's parents
dag_children(dag, "b")

## -----------------------------------------------------------------------------
dag_offspring(dag, "b")
dag_ancestors(dag, "e")

## -----------------------------------------------------------------------------
dag_offspring(dag, "b", include_self = TRUE)
dag_ancestors(dag, "e", include_self = TRUE)

## -----------------------------------------------------------------------------
dag_depth(dag)
dag_shortest_dist_from_root(dag)

## -----------------------------------------------------------------------------
dag_height(dag)
dag_shortest_dist_to_leaves(dag)

## -----------------------------------------------------------------------------
g = dag_as_igraph(dag)
g

## -----------------------------------------------------------------------------
library(igraph)
plot(g, layout = layout.sugiyama)

## -----------------------------------------------------------------------------
tree = dag_treelize(dag)
dag_depth(dag)
dag_depth(tree)

## -----------------------------------------------------------------------------
tree

## -----------------------------------------------------------------------------
dend = dag_as_dendrogram(tree)
dend = dendrapply(dend, function(d) {
    attr(d, "nodePar") = list(pch = attr(d, "label"))
    d
})
plot(dend, leaflab = "none")

## -----------------------------------------------------------------------------
relations = c("is_a", "is_a", "part_of", "part_of", "is_a", "is_a")
dag = create_ontology_DAG(parents, children, relations = relations)
dag

## -----------------------------------------------------------------------------
annotation = list(
    "a" = c("t1", "t2", "t3"),
    "b" = c("t3", "t4"),
    "c" = c("t5"),
    "d" = c("t7"),
    "e" = c("t4", "t5", "t6", "t7"),
    "f" = c("t8")
)
dag = create_ontology_DAG(parents, children, annotation = annotation)
dag

## -----------------------------------------------------------------------------
n_annotations(dag)

## -----------------------------------------------------------------------------
term_annotations(dag, letters[1:6])
annotated_terms(dag, c("t1", "t2", "t3"))

## -----------------------------------------------------------------------------
term_annotations(dag, letters[1:6], return = "matrix")
annotated_terms(dag, c("t1", "t2", "t3"), return = "matrix")

## ----message = TRUE-----------------------------------------------------------
parents  = c("a", "a", "b", "x", "x", "y")
children = c("b", "c", "c", "z", "y", "z")
create_ontology_DAG(parents, children)

## -----------------------------------------------------------------------------
# or with the double bracket: dag[["b"]]
dag["b"]

## -----------------------------------------------------------------------------
sessionInfo()

