## ----echo = FALSE, message = FALSE--------------------------------------------
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.align = "center")

## ----warning = FALSE----------------------------------------------------------
library(simona)
dag1 = import_obo("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo")
dag1

## -----------------------------------------------------------------------------
head(mcols(dag1))

## ----eval = Sys.info()["user"] == "guz"---------------------------------------
dag2 = import_ontology("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.owl", 
    robot_jar = "~/Downloads/robot.jar")

## ----eval = FALSE-------------------------------------------------------------
#  dag2

## ----echo = FALSE-------------------------------------------------------------
if(Sys.info()["user"] == "guz") {
    print(dag2)
} else {
    cat(
"An ontology_DAG object:
  Source: po, releases/2021-08-13
  1654 terms / 2510 relations
  Root: _all_
  Terms: PO:0000001, PO:0000002, PO:0000003, PO:0000004, ...
  Max depth: 13
  Aspect ratio: 24.85:1 (based on the longest distance to root)
                39.6:1 (based on the shortest distance to root)
  Relations: is_a, part_of

With the following columns in the metadata data frame:
  id, short_id, name, namespace, definition
")
}

## -----------------------------------------------------------------------------
dag3 = import_owl("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.owl")
dag3

## ----eval = FALSE-------------------------------------------------------------
#  # https://bioportal.bioontology.org/ontologies/MSTDE
#  dag4 = import_ttl("https://jokergoo.github.io/simona/MSTDE.ttl")
#  dag4

## -----------------------------------------------------------------------------
sessionInfo()

