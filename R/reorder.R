
# reorder children based on the cross-cluster links in each set of offsprings

#' @importFrom TSP TSP solve_TSP
#' @importFrom stats as.dist
reorder_children = function(dag, node) {

	children = dag_children(dag, node, in_labels = FALSE)
	n_children = length(children)
	if(n_children == 0) {
		return(integer(0))
	} else if(n_children <= 2) {
		return(children)
	}

	xl = lapply(children, function(x) {
		dag_offspring(dag, x, in_labels = FALSE, include_self = TRUE)
	})

	m = matrix(nrow = n_children, ncol = n_children)
	for(i in seq(1, n_children-1)) {
		for(j in seq(i+1, n_children)) {
			m[i, j] = m[j, i] = n_links_from_two_groups_of_nodes(dag, xl[[i]], xl[[j]])
		}
	}
	n_relations = sum(sapply(dag@lt_children, length))
	m = n_relations - m
	d = as.dist(m)
	tsp = TSP(d)
	fit = solve_TSP(tsp)

	children[as.vector(fit)]

}


#' Reorder child terms
#' 
#' @param dag An `ontology_Dag` object.
#' @param max_level Maximal depth of terms in DAG to apply reordering. The value
#'      should be set to a small integer because normally the effect of reordering
#'      too deep in the DAG is not noticable.
#' 
#' @details 
#' For a given term, its child terms are reordered based on the numbers
#' of cross-links of their offspring terms.
#' 
#' Denote c1 and c2 are two child terms, S1 and S2 are two sets of offspring terms of c1 and c2.
#' We calculate k as the number of links between S1 and S2, which are the links linking parents in S1
#' and children in S2, or parents in S2 and children in S1. Denote n as the total number of links
#' in the DAG, n - k is the distance bewteen c1 and c2 regarding how close the offspring sub-trees are connected.
#' 
#' With the distances from the set of child terms, [`TSP::TSP()`] is applied to reorder child terms.
#' 
#' @return An `ontology_Dag` object.
#' @export
#' @importFrom GetoptLong qqcat
#' @examples
#' \donttest{
#' dag = create_ontology_DAG_from_GO_db()
#' dag_reorder(dag)
#' }
#' 1
dag_reorder = function(dag, max_level = 2) {
	max_adjust_level = max_level
	current_terms = dag_root(dag, in_labels = FALSE)
	current_level = 0
	while(current_level < max_adjust_level) {
		message(qq("adjust child terms on level @{current_level}, @{length(current_terms)} terms"))
		current_terms2 = integer(0)
		for(t in current_terms) {
			dag@lt_children[[t]] = reorder_children(dag, t)
			current_terms2 = c(current_terms2, dag@lt_children[[t]])
		}
		current_terms = current_terms2
		current_level = current_level + 1
	}
	dag
}