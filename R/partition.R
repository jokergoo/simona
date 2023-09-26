
#' Partition the DAG
#' 
#' @param dag An `ontology_DAG` object.
#' @param level Depth in the DAG to cut. The DAG is cut below terms (or cut the links to their child terms) with `depth == level`.
#' @param from A list of terms to cut. If it is set, `level` is ignored.
#' @param term_pos Internally used.
#' 
#' @details
#' Let's call the terms below the `from` term as "top terms" because they will be on top of the sub-DAGs after the partitioning.
#' It is possible that a term in the middle of the DAG can be traced back to more than one top terms.
#' To partition all terms exclusively, a term partitioned to the sub-DAG from the top term with the largest distance to the term.
#' If a term has the same largest distances to several top terms, a random top term is selected.
#' 
#' In `partition_by_size()`, the DAG is first reduced to a tree where a child term only has one parent. 
#' The partition is done recursively by cutting into its child-trees. 
#' The splitting stops when all the child-trees have size less than `size`.
#' 
#' `NA` is assigned to the `from` terms, their ancestor terms, and terms having infinite directed distance to `from` terms.
#' 
#' @export
#' @returns A character vector of top terms in each partition.
#' @examples
#' \donttest{
#' dag = create_ontology_DAG_from_GO_db()
#' pa = partition_by_level(dag)
#' table(pa)
#' pa = partition_by_size(dag, size = 1000)
#' table(pa)
#' }
#' 1
partition_by_level = function(dag, level = 0, from = NULL, term_pos = NULL) {

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

	if(is.null(term_pos)) {
		if(dag_is_tree(dag)) {
			tree = dag
		} else {
			tree = dag_treelize(dag)
		}
		term_pos = cpp_term_pos_on_circle(tree, n_offspring(dag), 0, 360) ## in polar coordinate
	}

	from = from[order(term_pos[from, "rho"])]
	range = data.frame(left = term_pos[from, "theta"] - term_pos[from, "width"]/2,
		               right = term_pos[from, "theta"] + term_pos[from, "width"]/2)

	partition = rep(NA_character_, dag@n_terms)
	all_offspring = setdiff(seq_len(dag@n_terms), dag_ancestors(dag, from, in_labels = FALSE))
	l_offspring = rep(FALSE, dag@n_terms)
	l_offspring[all_offspring] =  TRUE
	for(i in seq_along(from)) {
		l = term_pos$theta >= range$left[i] & term_pos$theta <= range$right[i] & l_offspring
		partition[l] =  dag@terms[ from[i] ]
	}

	partition
}

#' @param size Number of terms in a cluster. The splitting stops on a term if all its child-tree are smaller than `size`.
#' @rdname partition_by_level
#' @importFrom stats dendrapply
#' @export
partition_by_size = function(dag, size = ceiling(dag_n_terms(dag)/10)) {
	
	tree = dag_treelize(dag)
	dend = dag_as_dendrogram(tree)

	pa = rep(NA_integer_, dag@n_terms)

	get_nodes_attr = function(dend, attr) {
		v = vector("list", attr(dend, "n_nodes"))
		i = 1
		dendrapply(dend, function(d) {
			v[[i]] <<- attr(d, attr)
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

		nn = vapply(dend, function(x) attr(x, "n_nodes"), FUN.VALUE = double(1))
		for(i in seq_along(nn)) {
			if(nn[i] > size) {
				scan_downstream(dend[[i]])
			} else {
				pa[get_nodes_attr(dend[[i]], "term_id")] <<- group
			}
		}
	}

	scan_downstream(dend)

	structure(tree@terms[pa], names = tree@terms)
}


