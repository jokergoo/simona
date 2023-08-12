

#' @importFrom matrixStats rowMaxs colMaxs
.GroupSim_pairwise = function(dag, group1, group2, term_sim_method, group_sim_method = "avg", ...) {

	group1 = lapply(group1, function(x) term_to_node_id(dag, x, strict = FALSE))
	group2 = lapply(group2, function(x) term_to_node_id(dag, x, strict = FALSE))

	id = unique(unlist(c(group1, group2)))
	sim = term_sim(dag, id, term_sim_method, control = list(...))
	rownames(sim) = colnames(sim) = dag@terms[id]

	.calc = function(g1, g2) {
		if(length(g1) == 0 || length(g2) == 0) {
			return(0)
		}
		sim = sim[g1, g2, drop = FALSE]

		group_sim_method = tolower(group_sim_method)

		best_match1 = rowMaxs(sim)
		best_match2 = colMaxs(sim)

		# avg
		if(group_sim_method == "avg") {
			mean(sim)
		} else if(group_sim_method == "max") {
			# max
			max(sim)
		} else if(group_sim_method == "bma") {
			#bma
			(mean(best_match1) + mean(best_match2))/2
		} else if(group_sim_method == "bmm") {
			# bmm
			max(mean(best_match1), mean(best_match2))
		} else if(group_sim_method == "abm") {
			# abm
			(sum(best_match1) + sum(best_match2))/(nrow(sim) + ncol(sim))
		} else if(group_sim_method == "hdf") {
			# hdf
			1 - max(1 - min(best_match1), 1 - min(best_match2))
		} else if(group_sim_method == "mhdf") {
			# MHDF
			1 - max(1 - mean(best_match1), 1 - mean(best_match2))
		} else if(group_sim_method == "vhdf") {
			1 - 0.5 * (sqrt( mean((1 - best_match1)^2)) + sqrt( mean((1 - best_match2)^2)))
		} else if(group_sim_method == "froehlich_2007") {
			exp(- max(1 - min(best_match1), 1 - min(best_match2)))
		} else if(group_sim_method == "joeng_2014") {
			0.5 * (sqrt( mean(best_match1^2)) + sqrt( mean(best_match2^2)))
		} else {
			stop("method not supported.")
		}
	}

	n1 = length(group1)
	n2 = length(group2)
	m = matrix(nrow = n1, ncol = n2)
	rownames(m) = names(group1)
	colnames(m) = names(group2)

	for(i in seq(1, length(group1))) {
		for(j in seq(1, length(group2))) {
			m[i, j] = .calc(dag@terms[ group1[[i]] ], dag@terms[ group2[[j]] ])
		}
	}
	m
}

#' GroupSim_pairwise_avg
#' 
#' @section Methods:
#' ## GroupSim_pairwise_avg
#' 
#' Denote `S(a, b)` as the semantic similarity between terms `a` and `b` where `a` is from `group1` and `b` is from `group2`, 
#' The similarity between `group1` and `group2` is the average similarity of every pair of individual terms in the two groups:
#' 
#' ```
#' group_sim = mean_{a in group1, b in group2}(S(a, b))
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_avg
GroupSim_pairwise_avg = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "avg", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_avg", c("term_sim_method", "..."))


#' GroupSim_pairwise_max
#' 
#' @section Methods:
#' ## GroupSim_pairwise_max
#' 
#' This is the maximal `S(a, b)` among all pairs of terms in `group1` and `group2`:
#' 
#' ```
#' group_sim = max_{a in group1, b in group2}(S(a, b))
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_max
GroupSim_pairwise_max = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "max", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_max", c("term_sim_method", "..."))

#' GroupSim_pairwise_BMA
#' 
#' @section Methods:
#' ## GroupSim_pairwise_BMA
#' 
#' BMA stands for "best-match average". First define similarity of a term to a group of terms as
#' 
#' ```
#' S(x, group) = max_{y in group}(x, y)
#' ```
#' 
#' which is the most similar terms in `group` to `x`.
#' 
#' Then the BMA similarity is calculated as:
#' 
#' ```
#' group_sim = 0.5*(mean_{a in group1}(S(a, group2)) + mean_{b in group2}(S(b, group1)))
#' ```
#' 
#' So it is the average of the similarity of every term in `group1` to the whole `group2` and every term in `group2` to the whole `group1`.
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_BMA
GroupSim_pairwise_BMA = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "BMA", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_BMA", c("term_sim_method", "..."))


#' GroupSim_pairwise_BMM
#' 
#' @section Methods:
#' ## GroupSim_pairwise_BMM
#' 
#' BMM stands for "best-match max". It is defined as:
#' 
#' ```
#' group_sim = max(mean_{a in group1}(S(a, group2)), mean_{b in group2}(S(b, group1)))
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_BMM
GroupSim_pairwise_BMM = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "BMM", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_BMM", c("term_sim_method", "..."))


#' GroupSim_pairwise_ABM
#' 
#' @section Methods:
#' ## GroupSim_pairwise_ABM
#' 
#' ABM stands for "average best-match". It is defined as:
#' 
#' ```
#' group_sim = (sum_{a in group1}(S(a, group2)) + sum_{b in group2}(S(b, group1)))/(n1 + n2)
#' ```
#' 
#' where `n1` and `n2` are the number of terms in `group1` and `group2`.
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_ABM
GroupSim_pairwise_ABM = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "ABM", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_ABM", c("term_sim_method", "..."))


#' GroupSim_pairwise_HDF
#' 
#' @section Methods:
#' ## GroupSim_pairwise_HDF
#' 
#' First define the distance of a term to a group of terms:
#' 
#' ```
#' D(x, group) = 1 - S(x, group)
#' ```
#' 
#' Then the Hausdorff distance between two groups are:
#' 
#' ```
#' HDF = max(max_{a in group1}(D(a, group2)), max_{b in group2}(D(b, group1)))
#' ```
#' 
#' This final similarity is:
#' 
#' ```
#' group_sim = 1 - HDF
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_HDF
GroupSim_pairwise_HDF = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "HDF", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_HDF", c("term_sim_method", "..."))


#' GroupSim_pairwise_MHDF
#' 
#' @section Methods:
#' ## GroupSim_pairwise_MHDF
#' 
#' Instead of using the maximal distance from a group to the other group, MHDF uses mean distance:
#' 
#' ```
#' MHDF = max(mean_{a in group1}(D(a, group2)), mean_{b in group2}(D(b, group1)))
#' ```
#' 
#' This final similarity is:
#' 
#' ```
#' group_sim = 1 - MHDF
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_MHDF
GroupSim_pairwise_MHDF = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "MHDF", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_MHDF", c("term_sim_method", "..."))



#' GroupSim_pairwise_VHDF
#' 
#' @section Methods:
#' ## GroupSim_pairwise_VHDF
#' 
#' It is defined as:
#' 
#' ```
#' VHDF = 0.5*(sqrt(mean_{a in group1}(D(a, group2)^2)) + sqrt(mean_{b in group2}(D(b, group1)^2)))
#' group_sim = 1 - VHDF
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_VHDF
GroupSim_pairwise_VHDF = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "VHDF", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_VHDF", c("term_sim_method", "..."))


#' GroupSim_pairwise_Froehlich_2007
#' 
#' @section Methods:
#' ## GroupSim_pairwise_Froehlich_2007
#' 
#' The similarity is:
#' 
#' ```
#' group_sim = exp(-HDF(group1, group2))
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_Froehlich_2007
GroupSim_pairwise_Froehlich_2007 = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "froehlich_2007", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_Froehlich_2007", c("term_sim_method", "..."))


#' GroupSim_pairwise_Joeng_2014
#' 
#' @section Methods:
#' ## GroupSim_pairwise_Joeng_2014
#' 
#' Similar to *VHDF*, but it directly uses the similarity:
#' 
#' ```
#' group_sim = 0.5*(sqrt(mean_{a in group1}(S(a, group2)^2)) + sqrt(mean_{b in group2}(S(b, group1)^2)))
#' ```
#' 
#' The term semantic similarity method can be set via `control = list(term_sim_method = ...)`.
#' 
#' @rdname temp__GroupSim_pairwise_Joeng_2014
GroupSim_pairwise_Joeng_2014 = function(dag, group1, group2, term_sim_method, ...) {
	.GroupSim_pairwise(dag, group1, group2, term_sim_method, "joeng_2014", ...)
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_Joeng_2014", c("term_sim_method", "..."))


#' GroupSim_SimALN
#' 
#' @section Methods:
#' ## GroupSim_SimALN
#' 
#' It is based on the average distances between every pair of terms in the two groups:
#' 
#' ```
#' exp(-mean_{a in group1, b in group2}(d(a, b)))
#' ```
#' 
#' `d(a, b)` is the distance between `a` and `b`, which can be the shortest distance between the two terms or
#' the longest distnace via LCA.
#' 
#' Distance method can be set via `control = list(distance = ...)`.
#' 
#' @rdname temp__GroupSim_SimALN
GroupSim_SimALN = function(dag, group1, group2, distance = "longest_distances_via_LCA") {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	id = union(id1, id2)
	if(distance == "shortest_distances_via_CA") {
		d = cpp_shortest_distances_via_CA(dag, id)
	} else if(distance == "longest_distances_via_LCA") {
		d = cpp_longest_distances_via_LCA(dag, id)
	} else {
		stop("`distance` can only be in 'shortest_distances_via_CA' or 'longest_distances_via_LCA'.")
	}
	dimnames(d) = list(dag@terms[id], dag@terms[id])

	d = d[dag@terms[id1], dag@terms[id2], drop = FALSE]
	exp(-mean(d))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimALN", "distance")


#' GroupSim_SimINT
#' 
#' @section Methods:
#' ## GroupSim_SimINT
#' 
#' To be added.
#' 
#' @rdname temp__GroupSim_SimINT
GroupSim_SimINT = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)
	ic_anno = IC_annotation(dag)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[id1]^2)/sum(ic_anno[id2]^2)
}
ADD_GROUP_SIM_METHOD("GroupSim_SimINT")


#' GroupSim_spgk
#' 
#' @section Methods:
#' ## GroupSim_spgk
#' 
#' To be added.
#' 
#' @rdname temp__GroupSim_spgk
GroupSim_spgk = function(dag, group1, group2) {
	
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	d1 = cpp_shortest_distances_directed(dag, ancestors1)
	dimnames(d1) = list(dag@terms[ancestors1], dag@terms[ancestors1])
	d2 = cpp_shortest_distances_directed(dag, ancestors2)
	dimnames(d2) = list(dag@terms[ancestors2], dag@terms[ancestors2])

	common_ancestors = intersect(ancestors1, ancestors2)

	if(length(common_ancestors) == 0) {
		0
	} else {
		common_ancestors = dag@terms[common_ancestors]
		sum(pmax(0, 2 - abs(d1[common_ancestors, common_ancestors, drop = FALSE] - d2[common_ancestors, common_ancestors, drop = FALSE])))
	}
}
ADD_GROUP_SIM_METHOD("GroupSim_spgk")


GROUPSIM_ANCESTORS_UNION = 1
GROUPSIM_ANCESTORS_INTERSECT = 2

#' GroupSim_SimGIC
#' 
#' @section Methods:
#' ## GroupSim_SimGIC
#' 
#' Denote `A` and `B` as the two sets of ancestors terms of terms in `group1` and `group2` respectively,
#' the SimGIC is:
#' 
#' ```
#' group_sim = sum_{x in intersect(A, B)}(IC(x))/sum_{x in union(A, B)}(IC(x))
#' ```
#' 
#' IC method can be set via `control = list(IC_method = ...)`.
#' 
#' @rdname temp__GroupSim_SimGIC
GroupSim_SimGIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_UNION, TRUE)

	ic = term_IC(dag, IC_method)
	sum(ic[ancestors1])/sum(ic[ancestors2])
}
ADD_GROUP_SIM_METHOD("GroupSim_SimGIC", "IC_method")


#' GroupSim_SimDIC
#' 
#' @section Methods:
#' ## GroupSim_SimDIC
#' 
#' Similar as *GroupSim_SimGIC*, it calculates the Dice coeffcient:
#' 
#' ```
#' group_sim = 2*sum_{x in intersect(A, B)}(IC(x))/(sum_{x in A}(IC(x)) + sum_{x in B}(IC(x)))
#' ```
#' 
#' IC method can be set via `control = list(IC_method = ...)`.
#' 
#' @rdname temp__GroupSim_SimDIC
GroupSim_SimDIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors3 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	ic = term_IC(dag, IC_method)
	2*sum(ic[ancestors1])/(sum(ic[ancestors2]) + sum(ic[ancestors3]))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimDIC", "IC_method")


#' GroupSim_SimUIC
#' 
#' @section Methods:
#' ## GroupSim_SimUIC
#' 
#' Similar as *GroupSim_SimGIC*, it is calculated as:
#' 
#' ```
#' group_sim = sum_{x in intersect(A, B)}(IC(x))/max(sum_{x in A}(IC(x)), sum_{x in B}(IC(x)))
#' ```
#' 
#' IC method can be set via `control = list(IC_method = ...)`.
#' 
#' @rdname temp__GroupSim_SimUIC
GroupSim_SimUIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors3 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	ic = term_IC(dag, IC_method)
	sum(ic[ancestors1])/max(sum(ic[ancestors2]), sum(ic[ancestors3]))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUIC", "IC_method")


#' GroupSim_SimUI
#' 
#' @section Methods:
#' ## GroupSim_SimUI
#' 
#' It is only based on the number of terms. `A` is the set of all ancestors of group1 terms and `B`
#' is the set of all ancestors of group2 terms.
#' 
#' ```
#' group_sim = length(intersect(A, B))/length(union(A, B))
#' ```
#' 
#' @rdname temp__GroupSim_SimUI
GroupSim_SimUI = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_UNION, TRUE)

	length(ancestors1)/length(ancestors2)
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUI")


#' GroupSim_SimDB
#' 
#' @section Methods:
#' ## GroupSim_SimDB
#' 
#' It is:
#' 
#' ```
#' group_sim = 2*length(intersect(A, B))/(length(A) + length(B))
#' ```
#' 
#' @rdname temp__GroupSim_SimDB
GroupSim_SimDB = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors3 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	2*length(ancestors1)/(length(ancestors2) + length(ancestors3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimDB")


#' GroupSim_SimUB
#' 
#' @section Methods:
#' ## GroupSim_SimUB
#' 
#' It is:
#' 
#' ```
#' group_sim = length(intersect(A, B))/max(length(A), length(B))
#' ```
#' 
#' @rdname temp__GroupSim_SimUB
GroupSim_SimUB = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors3 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	length(ancestors1)/max(length(ancestors2), length(ancestors3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUB")


#' GroupSim_SimNTO
#' 
#' @section Methods:
#' ## GroupSim_SimNTO
#' 
#' It is:
#' 
#' ```
#' group_sim = length(intersect(A, B))/min(length(A), length(B))
#' ```
#' 
#' @rdname temp__GroupSim_SimNTO
GroupSim_SimNTO = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestors1 = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	ancestors2 = cpp_ancestors_of_a_group(dag, id1, TRUE)
	ancestors3 = cpp_ancestors_of_a_group(dag, id2, TRUE)

	length(ancestors1)/min(length(ancestors2), length(ancestors3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimNTO")


#' GroupSim_SimCOU
#' 
#' @section Methods:
#' ## GroupSim_SimCOU
#' 
#' It is based on the dot product of two vectors `p` and `q` which correspond to terms in `group1` and `group2`. 
#' `p` and `q` have the same length as the total number of terms. Value of position i in `p` or `q` corresponds to term `t`. The value
#' takes `IC(t)` if `t` is an ancestor of any term in `p` or `q`, and the value takes zero if `t` is not. The
#' similarity betweem `group1` terms and `group2` terms is calculated as:
#' 
#' ```
#' <p,q>/||p||/||q||
#' ``` 
#' 
#' where `<p,q>` is the dot product between the two, and `||p||` or `||q||` is the norm of the vector.
#' The equation can be written as:
#' 
#' ```
#' group_sim = sum_{x in intersect(A, B)}(IC(x)^2) / 
#'               sqrt(sum_{x in A}(IC(x)^2)) / 
#'               sqrt(sum_{x in B}(IC(x)^2))
#' ```
#' 
#' IC method can be set via `control = list(IC_method = ...)`.
#' 
#' @rdname temp__GroupSim_SimCOU
GroupSim_SimCOU = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	ic = term_IC(dag, IC_method)

	sum(ic[i_intersect]^2)/sqrt(sum(ic[id1]^2))/sqrt(sum(ic[id2]^2))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimCOU")


#' GroupSim_SimCOT
#' 
#' @section Methods:
#' ## GroupSim_SimCOT
#' 
#' Similar as *GroupSim_SimCOU*, the similarity is:
#' 
#' ```
#' <p,q>/(||p||^2 + ||q||^2 - <p,q>)
#' ``` 
#' 
#' And it can be rewritten as:
#' 
#' ```
#' group_sim = sum_{x in intersect(A, B)}(IC(x)^2) /
#'     (sum_{x in A}(IC(x)^2) + sum_{x in B}(IC(x)^2) - sum_{x in intersect(A, B)}(IC(x)^2))
#' ```
#' 
#' IC method can be set via `control = list(IC_method = ...)`.
#' 
#' @rdname temp__GroupSim_SimCOT
GroupSim_SimCOT = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	ic = term_IC(dag, IC_method)

	sum(ic[i_intersect]^2)/( sum(ic[id1]^2) + sum(ic[id2]^2) - sum(ic[i_intersect]^2))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimCOT")


#' GroupSim_SimLP
#' 
#' @section Methods:
#' ## GroupSim_SimLP
#' 
#' It is the longest depth for the terms in `intersect(A, B)`.
#' 
#' 
#' ```
#' group_sim = max(depth(intersect(A, B)))
#' ```
#' 
#' @rdname temp__GroupSim_SimLP
GroupSim_SimLP = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	depth = dag_depth(dag)
	ancestors = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)
	max(depth[ancestors])
}
ADD_GROUP_SIM_METHOD("GroupSim_SimLP")


#' GroupSim_Ye_2005
#' 
#' @section Methods:
#' ## GroupSim_Ye_2005
#' 
#' It is a normalized version of *GroupSim_SimLP*:
#' 
#' ```
#' group_sim = max(depth(intersect(A, B)))/max_depth
#' ```
#' 
#' Since the minimal depth is zero for root.
#' 
#' @rdname temp__GroupSim_Ye_2005
GroupSim_Ye_2005 = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	depth = dag_depth(dag)
	ancestors = cpp_ancestors_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTORS_INTERSECT, TRUE)

	max(depth[ancestors])/max(depth)
}
ADD_GROUP_SIM_METHOD("GroupSim_Ye_2005")


#' GroupSim_SimCHO
#' 
#' @section Methods:
#' ## GroupSim_SimCHO
#' 
#' It is based on the annotated items. Denote `sigma(t)` as the total annotated items of `t`. The similarity is calculated as
#' 
#' ```
#' group_sim = log(C/sigma_max)/log(sigma_min/sigma_max)
#' ```
#' 
#' where `C` is `min(sigma_{x in intersect(A, B)}(x))`, i.e., the minimal sigma in the intersection of group1 and group2. Note
#' Now `A` and `B` are just two sets of terms in `group1` and `group2`.
#' `sigma_max` is the total number of items annotated to the DAG, `sig_min` is the minimal number of items annotated to a term, which
#' is mostly 1.
#' 
#' @rdname temp__GroupSim_SimCHO
GroupSim_SimCHO = function(dag, group1, group2) {

	n = n_annotations(dag)

	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)
	if(length(i_intersect) == 0) {
		0
	} else {
		log(min(n[i_intersect])/max(n))/log(min(n)/max(n))
	}
}
ADD_GROUP_SIM_METHOD("GroupSim_SimCHO")


#' GroupSim_SimALD
#' 
#' @section Methods:
#' ## GroupSim_SimALD
#' 
#' `A` and `B` are just two sets of terms in `group1` and `group2`. The similarity is calculated as:
#' 
#' ```
#' group_sim = max_{t in intersect(A, B)}(1 - sigma(t)/N)
#' ```
#' 
#' @rdname temp__GroupSim_SimALD
GroupSim_SimALD = function(dag, group1, group2) {

	n = n_annotations(dag)

	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	if(length(i_intersect) == 0) {
		0
	} else {
		max(1 - n[i_intersect]/max(n))
	}
}
ADD_GROUP_SIM_METHOD("GroupSim_SimALD")


#' GroupSim_Jaccard
#' 
#' @section Methods:
#' ## GroupSim_Jaccard
#' 
#' Say `A` is the set of items annotated to terms in `group1` and `B` is the set of items annotated to `group2`.
#' This is the Jaccard coeffcient between two sets.
#' 
#' 
#' The universe/background can be set via `control = list(universe = ...)`.
#' 
#' @rdname temp__GroupSim_Jaccard
GroupSim_Jaccard = function(dag, group1, group2, universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, universe, method = "jaccard")
}
ADD_GROUP_SIM_METHOD("GroupSim_Jaccard", "universe")


#' GroupSim_Dice
#' 
#' @section Methods:
#' ## GroupSim_Dice
#' 
#' It is the Dice coeffcient between `A` and `B`.
#' 
#' The universe/background can be set via `control = list(universe = ...)`.
#' 
#' @rdname temp__GroupSim_Dice
GroupSim_Dice = function(dag, group1, group2, universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, universe, method = "dice")
}
ADD_GROUP_SIM_METHOD("GroupSim_Dice", "universe")


#' GroupSim_Overlap
#' 
#' @section Methods:
#' ## GroupSim_Overlap
#' 
#' It is the Overlap coeffcient between `A` and `B`.
#' 
#' The universe/background can be set via `control = list(universe = ...)`.
#' 
#' @rdname temp__GroupSim_Overlap
GroupSim_Overlap = function(dag, group1, group2, universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, universe, method = "overlap")
}
ADD_GROUP_SIM_METHOD("GroupSim_Overlap", "universe")


#' GroupSim_Kappa
#' 
#' @section Methods:
#' ## GroupSim_Kappa
#' 
# It is the kappa coeffcient between `A` and `B`.
#' 
#' The universe/background can be set via `control = list(universe = ...)`.
#' 
#' @rdname temp__GroupSim_Kappa
GroupSim_Kappa = function(dag, group1, group2, universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, universe, method = "kappa")
}
ADD_GROUP_SIM_METHOD("GroupSim_Kappa", "universe")


.sim_overlap_from_two_groups = function(dag, id1, id2, universe = NULL, method = c("kappa", "jaccard", "dice", "overlap")) {

	check_pkg("proxyC", bioc = FALSE)

	if(!is.null(universe)) {
		universe = which(dag@terms %in% universe)
	} else {
		universe = seq_along(dag@terms)
	}

	if(length(dag@annotation$list) == 0) {
		stop("`annotation` should be set in `create_ontology_DAG()`.")
	}

	mg1 = cpp_get_term_annotations(dag, id1)
	mg1 = mg1[, universe, drop = FALSE]

	mg2 = cpp_get_term_annotations(dag, id2)
	mg2 = mg2[, universe, drop = FALSE]

	mg = rbind( (colSums(mg1) > 0) + 0, (colSums(mg2) > 0) + 0)
	mg = as(mg, "sparseMatrix")

	method = match.arg(method)[1]
	if(method == "kappa") {
		mat = kappa_dist(mg)
	} else if(method == "overlap") {
		mat = overlap_dist(mg)
	} else {
		mat = proxyC::simil(mg, method = method)
	}

	mat = as.matrix(mat)
	return(mat[1, 2])
}
