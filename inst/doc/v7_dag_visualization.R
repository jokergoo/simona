## ----echo = FALSE-------------------------------------------------------------
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)

knitr::opts_chunk$set(
  message = FALSE,
  dev = "ragg_png",
  fig.align = "center",
  pngquant = "--speed=10 --quality=30"
)

## -----------------------------------------------------------------------------
library(simona)
parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")
dag_small = create_ontology_DAG(parents, children)

## -----------------------------------------------------------------------------
dag_graphviz(dag_small)

## -----------------------------------------------------------------------------
color = sample(colors(), 6)
shape = c("polygon", "box", "oval", "egg", "diamond", "parallelogram")
dag_graphviz(dag_small, color = color, shape = shape)

## ----comment = ''-------------------------------------------------------------
dag_as_DOT(dag_small, color = color, shape = shape) |> cat()

## -----------------------------------------------------------------------------
dag = create_ontology_DAG_from_GO_db()
dag_graphviz(dag[, "GO:0010228"])

## -----------------------------------------------------------------------------
dag_graphviz(dag[, "GO:0010228"], 
  edge_color = c("is_a" = "purple", "part_of" = "darkgreen"),
  edge_style = c("is_a" = "solid", "part_of" = "dashed"), width = 800, height = 800)

## ----fig.width = 10, fig.height = 7-------------------------------------------
dag_circular_viz(dag)

## ----fig.width = 10, fig.height = 7-------------------------------------------
head(mcols(dag))
dag_circular_viz(dag, legend_labels_from = "name")

## ----fig.width = 10, fig.height = 7-------------------------------------------
tree = dag_treelize(dag)
dag_circular_viz(tree, legend_labels_from = "name", edge_transparency = 0.95)

## ----fig.width = 10, fig.height = 7-------------------------------------------
sig_go_ids = readRDS(system.file("extdata", "sig_go_ids.rds", package = "simona"))
dag_circular_viz(dag, highlight = sig_go_ids, legend_labels_from = "name")

## -----------------------------------------------------------------------------
sessionInfo()

