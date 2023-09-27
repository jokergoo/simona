

#' Reduce the DAG to a tree
#' 
#' @param dag An `ontology_DAG` object.
#' 
#' @details
#' A tree is a reduced DAG where a child only has one parent. The reducing is applied by a breadth-first searching
#' 
#' Starting from the root and on a certain depth (the depth is the maximal distance to root), for every term `t` on this depth,
#' its child term `c` and parent-child relation are kept only when `depth(c) == depth(t) + 1`. If `c` is selected, it is
#' marked as visited and will not be checked again.
#' 
#' In this way, depths of all terms in the orignal DAG are still identical to the depths in the tree (see the Examples section).
#' 
#' @returns A tree is also an `ontology_DAG` object.
#' @export
#' @examples
#' parents  = c("a", "a", "b", "b", "c", "d")
#' children = c("b", "c", "c", "d", "e", "f")
#' dag = create_ontology_DAG(parents, children)
#' tree = dag_treelize(dag)
#' d1 = dag_depth(dag)
#' d2 = dag_depth(tree)
#' identical(d1, d2)
#' 
#' dend = dag_as_dendrogram(tree)
#' dend
dag_treelize = function(dag) {

	if(dag_is_tree(dag)) {
		return(dag)
	}

	lt_children2 = cpp_mark_tree_links(dag)
	lt_children2 = lapply(lt_children2, function(x) {
		abs(x[x < 0])
	})

	n = length(lt_children2)
	parents = rep(seq_len(n), times = vapply(lt_children2, length, FUN.VALUE = integer(1)))
	children = unlist(lt_children2)

	tree = create_ontology_DAG(dag@terms[parents], dag@terms[children])
	tree@annotation = dag@annotation

	if(!is.null(mcols(dag))) {
		mcols(tree) = mcols(dag)[tree@terms, , drop = FALSE]
	}
	tree
}

#' @rdname dag_treelize
#' @export
#' @details
#' `dag_as_dendrogram()` coverts the tree to a `dendrogram` object.
dag_as_dendrogram = function(dag) {

	if(!dag_is_tree(dag)) {
		stop("dag is not a tree. run `dag_treelize()` first.")
	}

	add_dend = function(node, lt_children, height) {
		if(length(lt_children[[node]]) == 0) {			
			dend = node
			attributes(dend) = list(members = 1, 
				                    height = height, 
				                    label = dag@terms[node], 
				                    term_id = node, 
				                    n_nodes = 1,
				                    max_dist_to_leaf = 0,
				                    leaf = TRUE, 
				                    class = "dendrogram")
		} else {
			dend = lapply(lt_children[[node]], add_dend, lt_children, height - 1)
			attributes(dend) = list(members = length(unlist(dend)), 
				                    height = height, 
				                    label = dag@terms[node], 
				                    term_id = node, 
				                    n_nodes = 1 + sum(vapply(dend, function(x) attr(x, "n_nodes"), FUN.VALUE = numeric(1))),
				                    max_dist_to_leaf = 1 + max(vapply(dend, function(x) attr(x, "max_dist_to_leaf"), FUN.VALUE = numeric(1))),
				                    leaf = FALSE, 
				                    class = "dendrogram")
		}
		dend
	}

	max_depth = max(dag_depth(dag))
	dend = add_dend(dag@root, dag@lt_children, max_depth)

	dend = dend_set_midpoint(dend)
	class(dend) = c("ontology_tree", "dendrogram")

	dend
}


#' @param x An `ontology_DAG` object.
#' @param ... Ignored.
#' @rdname dag_treelize
#' @export
print.ontology_tree = function(x, ...) {
	cat("This is a normal `dendrogram` object, but with several more attributes attached to each node:\n")
	cat("  - label: term name\n")
	cat("  - term_id: the integer index of the term\n")
	cat("  - n_nodes: number of offspring terms, including itself\n")
	cat("  - members: this is the standard attribute, which is the number of leaves the term can reach to\n")
	cat("  - max_dist_to_leaf: max distance to the term's leaves\n")
	cat(" \n")
	getFromNamespace("print.dendrogram", "stats")(x)
}

dend_add_x = function(dend) {

	env = new.env(parent = emptyenv())
	i_leaf = 0

	add_x = function(dend) { # current sub-dend
		if(attr(dend, "leaf")) {
			i_leaf <<- i_leaf + 1
			x = i_leaf
			attr(dend, "x") = x
			attr(dend, "width") = 0
		} else {
			k = length(dend) # number of sub-dend
			for(i in seq_len(k)) {
				dend[[i]] = add_x(dend[[i]])
			}
			x = (attr(dend[[1]], "x") + attr(dend[[k]], "x"))/2
			attr(dend, "x") = x
			attr(dend, "width") = attr(dend[[k]], "x") - attr(dend[[1]], "x")
		}
		dend
	}

	add_x(dend)
}

dend_set_midpoint = function(dend) {

	dend = dend_add_x(dend)
	
	add_midpoint = function(dend) {
		if(attr(dend, "leaf")) {
			attr(dend, "midpoint") = 0
		} else {
			k = length(dend)
			for(i in seq_len(k)) {
				dend[[i]] = add_midpoint(dend[[i]])
			}
			midpoint = attr(dend[[1]], "midpoint") + attr(dend, "width")/2
			
			attr(dend, "midpoint") = midpoint
		}
		dend
	}

	add_midpoint(dend)
}


dend_sort = function(dend, by = "n_nodes", decreasing = TRUE) {

	if(!by %in% names(attributes(dend))) {
		stop("cannot find attribute:", by)
	}

	reorder = function(dend) {
		k = length(dend)
		if(k > 1) {
			v = vapply(dend, function(x) attr(x, by), FUN.VALUE = numeric(1))
			attr = attributes(dend)

			dend = dend[order(v, decreasing = decreasing)]
			attributes(dend) = attr
			attr(dend, "width") = attr(dend[[k]], "x") - attr(dend[[1]], "x")

			for(i in seq_len(k)) {
				dend[[i]] = reorder(dend[[i]])
			}
		}
		dend
	}

	dend = reorder(dend)
	dend_set_midpoint(dend)
}

