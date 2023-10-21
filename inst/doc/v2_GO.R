## ----echo = FALSE, message = FALSE--------------------------------------------
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.align = "center")

## -----------------------------------------------------------------------------
library(simona)
dag = create_ontology_DAG_from_GO_db("BP")
dag

## -----------------------------------------------------------------------------
create_ontology_DAG_from_GO_db("BP", relations = c("part of", "regulates"))  # "part_of" is also OK

## -----------------------------------------------------------------------------
library(org.Hs.eg.db)
dag = create_ontology_DAG_from_GO_db("BP", org_db = org.Hs.eg.db)
dag

## -----------------------------------------------------------------------------
term_annotations(dag, c("GO:0000002", "GO:0000012"))

## -----------------------------------------------------------------------------
head(mcols(dag))

## -----------------------------------------------------------------------------
sessionInfo()

