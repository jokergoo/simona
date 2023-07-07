
#' Partition the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param level Depth in the DAG to cut. The DAG is cut below terms (or cut the links to their child terms) with `depth == level`.
#' @param from A list of terms to cut. If it is set, `level` is ignored.
#' 
#' @details
#' Let's call the terms below the `from` term as "top terms" because they will be on top of the sub-DAGs after the partitioning.
#' It is possible that a term in the middle of the DAG can be traced back to more than one top terms.
#' To partition all terms exclusively, a term partitioned to the sub-DAG from the top term with the largest distance to the term.
#' If a term has the same largest distances to several top terms, a random top term is selected.
#' 
#' In `partition_by_k()`, the DAG is first reduced to a tree where a child term only has one parent. 
#' The partition is done recursively by cutting into its child-trees. 
#' The splitting stops when all the child-trees have size less than `k`.
#' 
#' `NA` is assigned to the `from` terms, their ancestor terms, and terms having infinite directed distance to `from` terms.
#' 
#' @export
#' @examples
#' \dontrun{
#' dag = create_ontology_DAG_from_GO_db()
#' pa = partition_by_level(dag)
#' table(pa)
#' pa = partition_by_k(dag, k = 1000)
#' table(pa)
#' }
partition_by_level = function(dag, level = 0, from = NULL) {

	if(is.null(from)) {
		depth = dag_depth(dag)
		max_depth = max(depth)
		if(level < 0 && level >= max_depth) {
			stop("wrong value of `level`.")
		}
		from = which(depth == level)
	} else {
		from = term_to_node_id(dag, from, strict = FALSE)
	}

	children = dag_children(dag, from, in_labels = FALSE)
	m = matrix(nrow = length(children), ncol = dag@n_terms)
	for(i in seq_along(children)) {
		m[i, ] = cpp_dag_longest_dist_to_offspring(dag, children[i])
	}

	ind = apply(m, 2, which.max)
	partition = dag@terms[ children[ind] ]
	partition[apply(m < 0, 2, all)] = NA
	structure(partition, names = dag@terms)
}

#' @param k Number of terms in a cluster. The splitting stops on a term if all its child-tree are smaller than `k`.
#' @rdname partition_by_level
#' @importFrom stats dendrapply
partition_by_k = function(dag, k) {
	
	tree = dag_treelize(dag)
	dend = dag_as_dendrogram(tree)

	pa = rep(NA_integer_, dag@n_terms)

	get_nodes_attr = function(dend, attr) {
		v = vector("list", attr(dend, "n_nodes"))
		i = 1
		dendrapply(dend, function(d) {
			v[[i]] <<- attr(d, "attr")
			i <<- i + 1
		})

		unlist(v)
	}

	scan_downstream = function(dend) {

		group = attr(dend, "term_id")
		if(attr(dend, "leaf")) {
			pa[group] <<- group
			return(NULL)
		}

		nn = sapply(dend, function(x) attr(x, "n_nodes"))
		for(i in seq_along(nn)) {
			if(nn[i] > k) {
				scan_downstream(dend[[i]])
			} else {
				pa[get_nodes_attr(dend[[i]], "term_id")] <<- group
			}
		}
	}

	scan_downstream(dend)

	structure(tree@terms[pa], names = tree@terms)
}


