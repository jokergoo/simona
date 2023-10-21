## ----echo = FALSE, message = FALSE--------------------------------------------
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.align = "center")

## ----echo = FALSE-------------------------------------------------------------
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)

knitr::opts_chunk$set(
  dev = "ragg_png",
  fig.align = "center",
  pngquant = "--speed=10 --quality=30"
)

## -----------------------------------------------------------------------------
library(simona)
dag = create_ontology_DAG_from_GO_db(namespace = "BP", org_db = "org.Hs.eg.db")
dag

## -----------------------------------------------------------------------------
all_term_IC_methods()

## ----fig.width = 6, fig.height = 6--------------------------------------------
ic1 = term_IC(dag, method = "IC_annotation")
ic2 = term_IC(dag, method = "IC_annotation", control = list(uniquify = FALSE))
# ranges on both x- and y-axes
rg = c(0, max(ic1, ic2, na.rm = TRUE))
plot(ic1, ic2, xlim = rg, ylim = rg,
    xlab = "uniquified (first method)", ylab = "not uniquified (second method)", 
    pch = 16, col = "#00000020", main = "compare IC_annotation")
abline(a = 0, b = 1, col = "red")

## -----------------------------------------------------------------------------
sessionInfo()

