

#' @importFrom stats p.adjust
dag_enrich_terms = function(dag, terms) {
	n = dag@n_terms
	ind = which(dag@terms %in% terms)
	m = length(ind)

	n_offspring = n_offspring(dag, include_self = TRUE)
	n_hits = cpp_n_offspring_with_intersect(dag, ind, include_self = TRUE)

	p = phyper(n_hits-1, m, n - m, n_offspring, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	df = data.frame(term = dag@terms, n_hits = n_hits, n_offspring = n_offspring, n_terms = m, n_all = n, p_value = p, p_adjust = padj)
	df$depth = dag_depth(dag)
	df$height = dag_height(dag)
	df
}

dag_enrich_items = function(dag, p_value) {
	n = dag@n_terms

	lt_offspring = dag_all_offspring(dag, in_labels = FALSE)
	n_offspring = sapply(lt_offspring, length)

	s = sapply(lt_offspring, function(x) mean(-log(p_value[x])))
	sr = matrix(nrow = length(s), ncol = 1000)
	p = numeric(n)
	
	for(i in 1:1000) {
		cat(i, "/", 1000, "\n")
		sr[, i] = sapply(n_offspring, function(x) {
			mean(-log(p_value[sample(p_value, x)]))
		})
	}

	p = sapply(1:n, function(i) sum(sr[, i] - s[i])/length(s))

	data.frame(term = dag@terms, s = s, p = p)
}

#' @importFrom stats phyper
dag_enrich_items = function(dag, items) {
	validate_dag_has_annotation(dag)

	n = length(dag@annotation$names)
	ind = which(dag@annotation$names %in% items)
	m = length(ind)

	n_anno = n_annotations(dag)
	n_hits = cpp_n_annotations_with_intersect(dag, ind)

	p = phyper(n_hits-1, m, n - m, n_anno, lower.tail = FALSE)
	padj = p.adjust(p, "BH")
	df = data.frame(term = dag@terms, n_hits = n_hits, n_anno = n_anno, n_items = m, n_all = n, p_value = p, p_adjust = padj)
	df$depth = dag_depth(dag)
	df$height = dag_height(dag)
	df
}


random_terms = function(dag, n) {
	sample(dag@terms, n)
}

random_items = function(dag, n) {
	sample(dag@annotation$names, n)
}
