

#' Enrichment analysis on offspring terms
#' 
#' The analysis task is to evaluate how significant a term includes `terms`. 
#' 
#' @param dag An `ontology_DAG` object.
#' @param terms A vector of term names.
#' @param min_hits Minimal number of terms in an offspring set.
#' @param min_offspring Minimal size of the offspring set.
#' 
#' @details
#' Given a list of terms in `terms`, the function tests whether they are enriched in a term's offspring terms.
#' The test is based on the hypergeometric distribution. In the following 2x2 contigency table, `S` is the set of `terms`,
#' for a term `t` in the DAG, `T` is the set of its offspring plus the `t` itself, the aim is to test whether `S` is over-represented
#' in `T`.
#' 
#' If there is a significant p-value, we can say the term `t` preferably includes terms in `term`.
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
#' - `log2_fold_enrichment`: Defined as log2(observation/expected).
#' - `z_score`: Defined as (observed-expected)/sd.
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
#' df = dag_enrich_on_offsprings(dag, terms)
#' }
#' 1
dag_enrich_on_offsprings = function(dag, terms, min_hits = 3, min_offspring = 10) {
	n = dag@n_terms
	ind = which(dag@terms %in% terms)
	m = length(ind)

	n_offspring = n_offspring(dag, include_self = TRUE)
	n_hits = cpp_n_offspring_with_intersect(dag, ind, include_self = TRUE)
	all_terms = dag@terms

	l = n_offspring >= min_offspring & n_hits >= min_hits
	n_offspring = n_offspring[l]
	n_hits = n_hits[l]
	all_terms = all_terms[l]

	p = phyper(n_hits-1, m, n - m, n_offspring, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	if("name" %in% colnames(mcols(dag))) {
		df = data.frame(term = all_terms, name = mcols(dag)[l, "name"],
			n_hits = n_hits, n_offspring = n_offspring, n_terms = m, n_all = n, 
			log2_fold_enrichment = log2(n_hits/(m*n_offspring/n)), 
			z_score = (n_hits - (m*n_offspring/n)) / sqrt(m*n_offspring/n * (n-n_offspring)/n * (n-m)/(n-1)),
			p_value = p, p_adjust = padj)
	} else {
		df = data.frame(term = all_terms, n_hits = n_hits, n_offspring = n_offspring, n_terms = m, n_all = n, 
			log2_fold_enrichment = log2(n_hits/(m*n_offspring/n)), 
			z_score = (n_hits - (m*n_offspring/n)) / sqrt(m*n_offspring/n * (n-n_offspring)/n * (n-m)/(n-1)),
			p_value = p, p_adjust = padj)
	}
	df$depth = dag_depth(dag, all_terms)
	df
}

#' Enrichment analysis on offspring terms by permutation test
#' 
#' @param dag An `ontology_DAG` object.
#' @param value A numeric value. The value should correspond to terms in `dag@terms`.
#' @param min_offspring Minimal size of the offspring set.
#' @param perm Number of permutations.
#' @param verbose Whether to print messages.
#' 
#' @details
#' In the function [`dag_enrich_on_offsprings()`], the statistic for testing is the number of terms in each category. Here
#' this funtion makes the testing procedure more general
#' 
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
#' p = sum_{i in 1..k}(I(sr_i > s))/k
#' ```
#' 
#' @return A data frame with the following columns:
#' 
#' - `term`: Term names.
#' - `stats`: The statistics of terms.
#' - `n_offspring`: Number of offspring terms of `t` (including `t` itself).
#' - `log2_fold_enrichment`: defined as `log2(s/mean)` where `mean` is calculated from random permutation.
#' - `z_score`: Defined as `(s - mean)/sd` where `mean` and `sd` are calculated from random permutation.
#' - `p_value`: P-values from permutation test.
#' - `p_adjust`: Adjusted p-values from the BH method.
#' 
#' The number of rows in the data frame is the same as the number of terms in the DAG.
#' 
#' @importFrom stats sd
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db() 
#' value = runif(dag_n_terms(dag)) # a set of random values
#' df = dag_enrich_on_offsprings_by_permutation(dag, value)
#' }
#' 1
dag_enrich_on_offsprings_by_permutation = function(dag, value, perm = 1000, min_offspring = 10, verbose = simona_opt$verbose) {
	n = dag@n_terms

	if(length(value) != n) {
		stop("Length of `value` should be the same as the total number of terms.")
	}

	n_offspring = n_offspring(dag)
	s = exec_under_message_condition({
		cpp_offspring_aggregate(dag, value)
	}, verbose = verbose)
	all_terms = dag@terms

	l = n_offspring >= min_offspring
	n_offspring = n_offspring[l]
	s = s[l]
	all_terms = all_terms[l]

	sr = exec_under_message_condition({
		cpp_random_aggregatioin(n_offspring, value, perm)  # columns are permutations
	}, verbose = verbose)

	p = vapply(seq_len(nrow(sr)), function(i) sum(sr[i, ] > s[i])/perm, FUN.VALUE = numeric(1))
	padj = p.adjust(p, "BH")
	z = vapply(seq_len(nrow(sr)), function(i) (s[i] - mean(sr[i, ]))/sd(sr[i, ]), FUN.VALUE = numeric(1))
	log2fe = vapply(seq_len(nrow(sr)), function(i) log2(s[i]/mean(sr[i, ])), FUN.VALUE = numeric(1))

	if("name" %in% colnames(mcols(dag))) {
		df = data.frame(term = all_terms, name = mcols(dag)[l, "name"],
			stat = s, n_offspring = unname(n_offspring), 
			log2_fold_enrichment = log2fe, z_score = z, 
			p_value = p, p_adjust = padj)
	} else {
		df = data.frame(term = all_terms, stat = s, n_offspring = unname(n_offspring), 
			log2_fold_enrichment = log2fe, z_score = z, 
			p_value = p, p_adjust = padj)
	}
	df$depth = dag_depth(dag, all_terms)
	df
}

#' Enrichment analysis on the number of annotated items
#' 
#' The analysis task is to evaluate which terms the given items are enriched to.
#' 
#' @param dag An `ontology_DAG` object.
#' @param items A vector of item names.
#' @param min_hits Minimal number of items in the term set.
#' @param min_items Minimal size of the term set.
#' 
#' @details
#' The function tests whether the list of items are enriched in terms on the DAG.
#' The test is based on the hypergeometric distribution. In the following 2x2 contigency table, `S` is the set of `items`,
#' for a term `t` in the DAG, `T` is the set of items annotated to `t` (by automatically merging from its offspring terms), 
#' the aim is to test whether `S` is over-represented in `T`.
#' 
#' The universal set `all` correspond to the full set of items annotated to the DAG.
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
#' - `n_anno`: Number of annotated items of `t`. Specifically for `dag_enrich_on_genes()`, this column
#'             is renamed to `n_gs`.
#' - `n_items`: Number of items in `items` intersecting to all annotated items in the DAG. Specifically
#'             for `dag_enrich_on_genes()`, this column is renamed to `n_genes`.
#' - `n_all`: Number of all annotated items in the DAG.
#' - `log2_fold_enrichment`: Defined as log2(observation/expected).
#' - `z_score`: Defined as (observed-expected)/sd.
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
#' df = dag_enrich_on_items(dag, items)
#' }
#' 1
dag_enrich_on_items = function(dag, items, min_hits = 5, min_items = 10) {
	validate_dag_has_annotation(dag)

	n = length(dag@annotation$names)
	ind = which(dag@annotation$names %in% items)
	m = length(ind)

	n_anno = n_annotations(dag, uniquify = TRUE)
	n_hits = cpp_n_annotations_with_intersect(dag, ind)
	all_terms = dag@terms

	l = n_anno >= min_items & n_hits >= min_hits
	n_anno = n_anno[l]
	n_hits = n_hits[l]
	all_terms = all_terms[l]

	p = phyper(n_hits-1, m, n - m, n_anno, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	if("name" %in% colnames(mcols(dag))) {
		df = data.frame(term = all_terms, name = mcols(dag)[l, "name"],
			n_hits = n_hits, n_anno = n_anno, n_items = m, n_all = n, 
			log2_fold_enrichment = log2(n_hits/(m*n_anno/n)), 
			z_score = (n_hits - (m*n_anno/n)) / sqrt(m*n_anno/n * (n-n_anno)/n * (n-m)/(n-1)),
			p_value = p, p_adjust = padj)
	} else {
		df = data.frame(term = all_terms, n_hits = n_hits, n_anno = n_anno, n_items = m, n_all = n, 
			log2_fold_enrichment = log2(n_hits/(m*n_anno/n)), 
			z_score = (n_hits - (m*n_anno/n)) / sqrt(m*n_anno/n * (n-n_anno)/n * (n-m)/(n-1)),
			p_value = p, p_adjust = padj)
	}
	df$depth = dag_depth(dag, all_terms)
	df
}

#' @param genes A vector of gene IDs. The gene ID type can be found by directly printing the `ontology_DAG` object.
#' @param min_genes Minimal number of genes.
#' 
#' @details
#' `dag_enrich_on_genes()` is the same as `dag_enrich_on_items()` which only changes the argument `item` to `gene`.
#' 
#' @rdname dag_enrich_on_items
#' @export
dag_enrich_on_genes = function(dag, genes, min_hits = 5, min_genes = 10) {
	df = dag_enrich_on_items(dag, genes, min_hits = min_hits, min_items = min_genes)
	colnames(df)[colnames(df) == "n_items"]  = "n_genes"
	colnames(df)[colnames(df) == "n_anno"]  = "n_gs"
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
