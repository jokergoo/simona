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
all_term_sim_methods()

## -----------------------------------------------------------------------------
sessionInfo()

