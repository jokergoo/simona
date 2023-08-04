


## Comparisons

Let's make a comparison of various IC metrics, using the GO BP ontology.

```r
dag = create_ontology_DAG_from_GO_db(namespace = "BP", org_db = "org.Hs.eg.db")
lt = lapply(all_ic_methods(), function(method) {
    cat("=====", method, "=====\n")
    term_IC(dag, method)
})
names(lt) = all_ic_methods()

df = as.data.frame(lt)
pairs(df, pch = ".", col = dag_depth(dag), gap = 0)
```

```{r, echo = FALSE, fig.width = 10, fig.height = 10}
df = readRDS(system.file("extdata", "go_ic.rds", package = "simone"))
pairs(df, pch = ".", col = dag_depth(dag), gap = 0)
```

And the heatmap of the correlations of IC metrics. 

```{r, fig.width = 6, fig.height = 6}
cor = cor(df, method = "spearman", use = "pairwise.complete.obs")
library(ComplexHeatmap)
Heatmap(cor, name = "correlation", column_title = "Spearman correlation")
```



# Comparisons

Let's make a comparison of various similarity methods, using the GO BP ontology.

```r
dag = create_ontology_DAG_from_GO_db(namespace = "BP", org_db = "org.Hs.eg.db")
set.seed(123)
ic = term_IC(dag, method = "IC_annotation")
ic = ic[!is.na(ic)]
go_id = sample(names(ic), 100)

lt = lapply(all_term_sim_methods(), function(method) {
    cat("=====", method, "=====\n")
    m = term_sim(dag, go_id, method)
    m[lower.tri(m)]
})
names(lt) = all_term_sim_methods()
df = as.data.frame(lt)

LCA_depth = LCA_depth(dag, go_id)
LCA_depth = LCA_depth[lower.tri(LCA_depth)]

pairs(df[, 1:29], pch = ".", gap = 0, col = LCA_depth+1)
```

```{r, echo = FALSE, fig.width = 12, fig.height = 12, out.width = "100%"}
df = readRDS(system.file("extdata", "go_sim.rds", package = "simone"))
LCA_depth = readRDS(system.file("extdata", "LCA_depth.rds", package = "simone"))
pairs(df[, 1:29], pch = ".", gap = 0, col = LCA_depth+1)
```

And the heatmap of the correlations of IC metrics. 

```{r, fig.width = 8, fig.height = 8}
cor = cor(df[, 1:29], method = "spearman", use = "pairwise.complete.obs")
library(ComplexHeatmap)
Heatmap(cor, name = "correlation", column_title = "Spearman correlation")
```

