

#' Enrichment analysis on offspring terms
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' 
#' @details
#' Given a list of terms in `terms`, the function tests whether they are enriched in some terms's offspring terms.
#' The test is based on the hypergeometric distribution. In the following 2x2 contigency table, `S` is the set of `terms`,
#' for a term `t` in the DAG, `T` is the set of its offspring plus the `t` itself, the aim is to test whether `S` is over-represented
#' in `T`.
#' 
#' ```
#' +----------+------+----------+-----+
#' |          | in S | not in S | all |
#' +----------+------+----------+-----+
#' | in T     |  x11 |    x12   | x10 |
#' | not in T |  x21 |    x22   | x20 |
#' +----------+------+----------+-----+
#' | all      |  x01 |    x02   |  x  |
#' +----------+------+----------+-----+
#' ``` 
#' 
#' @return A data frame with the following columns:
#' 
#' - `term`: Term names.
#' - `n_hits`: Number of terms in `terms` intersecting to `t`'s offspring terms. 
#' - `n_offspring`: Number of offspring terms of `t` (including `t` itself).
#' - `n_terms`: Number of terms in `term` intersecting to all terms in the DAG.
#' - `n_all`: Number of all terms in the DAG.
#' - `p_value`: P-values from hypergeometric test.
#' - `p_adjust`: Adjusted p-values from the BH method.
#' 
#' The number of rows in the data frame is the same as the number of terms in the DAG.
#' 
#' @importFrom stats p.adjust phyper
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db() 
#' terms = random_terms(dag, 100)
#' df = dag_enrich_terms(dag, terms)
#' }
#' 1
dag_enrich_terms = function(dag, terms) {
	n = dag@n_terms
	ind = which(dag@terms %in% terms)
	m = length(ind)

	n_offspring = n_offspring(dag, include_self = TRUE)
	n_hits = cpp_n_offspring_with_intersect(dag, ind, include_self = TRUE)

	p = phyper(n_hits-1, m, n - m, n_offspring, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	df = data.frame(term = dag@terms, n_hits = n_hits, n_offspring = n_offspring, n_terms = m, n_all = n, p_value = p, p_adjust = padj)
	df
}

#' Enrichment analysis on offspring terms by permutation test
#' 
#' @param dag An `ontology_DAG` object.
#' @param value A numeric value. The value should correspond to terms in `dag@terms`.
#' @param n Number of permutations.
#' 
#' @details
#' The function tests whether a term `t`'s offspring terms have an over-represented pattern on values in `value`.
#' Denote `T` as the set of `t`'s offspring terms plus `t` itself, and `v` as the numeric vector of `value`, we first
#' calculate a score `s` based on values in `T`:
#' 
#' ```
#' s = mean_{terms in T}(v)
#' ```
#' 
#' To construct a random version of `s`, we randomly sample `n_T` terms from the DAG where `n_T` is the size of set `T`:
#' 
#' ```
#' sr_i = mean_{n_T randomly sampled terms}(v)
#' ```
#' 
#' where index `i` represents the i^th sampling. If we sample `k` times, the p-value is calculated as:
#' 
#' ```
#' p = sum_{go over k}(I(sr > s))/k
#' ```
#' 
#' @return A data frame with the following columns:
#' 
#' - `term`: Term names.
#' - `stats`: The statistics of terms.
#' - `n_offspring`: Number of offspring terms of `t` (including `t` itself).
#' - `z_score`: Defined as `(s - mean)/sd` where `mean` and `sd` are calculated from random permutation.
#' - `log2_fold_enrichment`: defined as `log2(s/mean)` where `mean` is calculated from random permutation.
#' - `p_value`: P-values from permutation test.
#' - `p_adjust`: Adjusted p-values from the BH method.
#' 
#' The number of rows in the data frame is the same as the number of terms in the DAG.
#' 
#' @importFrom stats sd
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db() 
#' value = runif(dag_n_terms(dag))
#' df = dag_enrich_terms_by_permutation(dag, value)
#' }
#' 1
dag_enrich_terms_by_permutation = function(dag, value, n = 1000) {
	n = dag@n_terms

	lt_offspring = dag_all_offspring(dag, include_self = TRUE, in_labels = FALSE)
	n_offspring = vapply(lt_offspring, length, FUN.VALUE = integer(1))

	s = vapply(lt_offspring, function(x) mean(value[x], na.rm = TRUE), FUN.VALUE = numeric(1))
	sr = matrix(nrow = length(s), ncol = n)
	p = numeric(n)
	
	for(i in seq_len(n)) {
		message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\brandom sampling", i, "/", n, appendLF = FALSE)
		sr[, i] = vapply(n_offspring, function(x) {
			mean(sample(value, x), na.rm = TRUE)
		}, FUN.VALUE = numeric(1))
	}
	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\brandom sampling", n, "/", n, "\n")

	p = vapply(seq_len(n), function(i) sum(sr[, i] > s[i])/length(s), FUN.VALUE = numeric(1))
	padj = p.adjust(p, "BH")
	z = vapply(seq_len(n), function(i) (s[i] - mean(sr[, i]))/sd(sr[, i]), FUN.VALUE = numeric(1))
	log2fe = vapply(seq_len(n), function(i) log2(s[i]/mean(sr[, i])), FUN.VALUE = numeric(1))

	data.frame(term = dag@terms, stat = s, n_offspring = unname(n_offspring), z_score = z, log2_fold_enrichment = log2fe, p_value = p, p_adjust = padj)
}

#' Enrichment analysis on numbers of annotated items
#' 
#' @param dag An `ontology_DAG` object.
#' @param items A vector of item names.
#' 
#' @details
#' The function tests whether the list of items are enriched in some terms.
#' The test is based on the hypergeometric distribution. In the following 2x2 contigency table, `S` is the set of `items`,
#' for a term `t` in the DAG, `T` is the set of items annotated to `t` (by automatically merging from its offspring terms), 
#' the aim is to test whether `S` is over-represented in `T`.
#' 
#' ```
#' +----------+------+----------+-----+
#' |          | in S | not in S | all |
#' +----------+------+----------+-----+
#' | in T     |  x11 |    x12   | x10 |
#' | not in T |  x21 |    x22   | x20 |
#' +----------+------+----------+-----+
#' | all      |  x01 |    x02   |  x  |
#' +----------+------+----------+-----+
#' ``` 
#' 
#' @return A data frame with the following columns:
#' 
#' - `term`: Term names.
#' - `n_hits`: Number of items in `items` intersecting to `t`'s annotated items.
#' - `n_anno`: Number of annotated items of `t`.
#' - `n_items`: Number of items in `items` intersecting to all annotated items in the DAG.
#' - `n_all`: Number of all annotated items in the DAG.
#' - `p_value`: P-values from hypergeometric test.
#' - `p_adjust`: Adjusted p-values from the BH method.
#' 
#' The number of rows in the data frame is the same as the number of terms in the DAG.
#' 
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db(org_db = "org.Hs.eg.db") 
#' items = random_items(dag, 1000)
#' df = dag_enrich_items(dag, items)
#' }
#' 1
dag_enrich_items = function(dag, items) {
	validate_dag_has_annotation(dag)

	n = length(dag@annotation$names)
	ind = which(dag@annotation$names %in% items)
	m = length(ind)

	n_anno = n_annotations(dag, uniquify = TRUE)
	n_hits = cpp_n_annotations_with_intersect(dag, ind)

	p = phyper(n_hits-1, m, n - m, n_anno, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	df = data.frame(term = dag@terms, n_hits = n_hits, n_anno = n_anno, n_items = m, n_all = n, p_value = p, p_adjust = padj)
	df$depth = dag_depth(dag)
	df$height = dag_height(dag)
	df
}

#' Randomly sample terms/items
#' 
#' @param dag An `ontology_DAG` object.
#' @param n Number of terms or items.
#' 
#' @export
#' @returns A character vector of terms or items.
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' annotation = list(
#'     "a" = c("t1", "t2", "t3"),
#'     "b" = c("t3", "t4"),
#'     "c" = "t5",
#'     "d" = "t7",
#'     "e" = c("t4", "t5", "t6", "t7"),
#'     "f" = "t8"
#' )
#' dag = create_ontology_DAG(parents, children, annotation = annotation)
#' random_terms(dag, 3)
#' random_items(dag, 3)
random_terms = function(dag, n) {
	if(n > dag@n_terms) {
		stop("`n` should not be larger than the total number of terms in the DAG.")
	}
	sample(dag@terms, n)
}

#' @rdname random_terms
#' @export
random_items = function(dag, n) {
	n_anno = length(dag@annotation$names)
	if(n_anno == 0) {
		stop("annotation has not been set in the DAG.")
	}
	if(n > n_anno) {
		stop("`n` should not be larger than the total number of annotated items in the DAG.")
	}
	sample(dag@annotation$names, n)
}
