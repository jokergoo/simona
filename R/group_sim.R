

#' @importFrom matrixStats rowMaxs colMaxs
.GroupSim_pairwise = function(dag, group1, group2, sim_method, group_sim_method = "avg") {

	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)
	group1 = dag@terms[id1]
	group2 = dag@terms[id2]

	group = union(group1, group2)
	sim = term_sim(dag, group, sim_method)
	sim = sim/max(sim)
	group1 = intersect(group1, rownames(sim))
	group2 = intersect(group2, rownames(sim))
	sim = sim[group1, group2, drop = FALSE]

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

GroupSim_pairwise_avg = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "avg")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_avg")


GroupSim_pairwise_max = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "max")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_max")


GroupSim_pairwise_BMM = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "BMM")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_BMM")


GroupSim_pairwise_ABM = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "ABM")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_ABM")


GroupSim_pairwise_HDF = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "HDF")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_HDF")


GroupSim_pairwise_VHDF = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "VHDF")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_VHDF")


GroupSim_pairwise_Froehlich_2007 = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "froehlich_2007")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_Froehlich_2007")


GroupSim_pairwise_Joeng_2014 = function(dag, group1, group2, sim_method) {
	.GroupSim_pairwise(dag, group1, group2, sim_method, "joeng_2014")
}
ADD_GROUP_SIM_METHOD("GroupSim_pairwise_Joeng_2014")


GroupSim_SimALN = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	id = union(id1, id2)
	d = cpp_longest_distances_via_LCA(dag, id)
	dimnames(d) = list(dag@terms[id], dag@terms[id])

	d = d[dag@terms[id1], dag@terms[id2], drop = FALSE]
	exp(-mean(d))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimALN")


GroupSim_SimINT = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)
	ic_anno = IC_annotation(dag)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[id1]^2)/sum(ic_anno[id2]^2)
}
ADD_GROUP_SIM_METHOD("GroupSim_SimINT")


GroupSim_spgk = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	d1 = cpp_shortest_distances_directed(dag, ancestor1)
	dimnames(d1) = list(dag@terms[ancestor1], dag@terms[ancestor1])
	d2 = cpp_shortest_distances_directed(dag, ancestor2)
	dimnames(d2) = list(dag@terms[ancestor2], dag@terms[ancestor2])

	common_ancestor = intersect(ancestor1, ancestor2)

	if(length(common_ancestor) == 0) {
		0
	} else {
		common_ancestor = dag@terms[common_ancestor]
		sum(pmax(0, 2 - abs(d1[common_ancestor, common_ancestor, drop = FALSE] - d2[common_ancestor, common_ancestor, drop = FALSE])))
	}
}
ADD_GROUP_SIM_METHOD("GroupSim_spgk")


GROUPSIM_ANCESTOR_UNION = 1
GROUPSIM_ANCESTOR_INTERSECT = 2

GroupSim_SimGIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_UNION, TRUE)

	ic = term_IC(dag, IC_method)
	sum(ic[ancestor1])/sum(ic[ancestor2])
}
ADD_GROUP_SIM_METHOD("GroupSim_SimGIC")


GroupSim_SimDIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	ic = term_IC(dag, IC_method)
	2*sum(ic[ancestor1])/(sum(ic[ancestor2]) + sum(ic[ancestor3]))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimDIC")


GroupSim_SimUIC = function(dag, group1, group2, IC_method) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	ic = term_IC(dag, IC_method)
	sum(ic[ancestor1])/max(sum(ic[ancestor2]), sum(ic[ancestor3]))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUIC")


GroupSim_SimUI = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_UNION, TRUE)

	length(ancestor1)/length(ancestor2)
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUI")


GroupSim_SimDB = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	2*length(ancestor1)/(length(ancestor2) + length(ancestor3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimDB")


GroupSim_SimUB = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	length(ancestor1)/max(length(ancestor2), length(ancestor3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimUB")


GroupSim_SimNTO = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	ancestor1 = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, id1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, id2, TRUE)

	length(ancestor1)/min(length(ancestor2), length(ancestor3))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimNTO")


GroupSim_SimCOU = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[id1]^2)/sum(ic_anno[id2]^2)
}
ADD_GROUP_SIM_METHOD("GroupSim_SimCOU")


GroupSim_SimCOT = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/( sum(ic_anno[id1]^2)/sum(ic_anno[id2]^2) - sum(ic_anno[i_intersect]^2))
}
ADD_GROUP_SIM_METHOD("GroupSim_SimCOT")


GroupSim_SimLP = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	depth = dag_depth(dag)
	ancestor = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	max(depth[ancestor])
}
ADD_GROUP_SIM_METHOD("GroupSim_SimLP")


GroupSim_Ye_2005 = function(dag, group1, group2) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	depth = dag_depth(dag)
	ancestor = cpp_ancestor_of_two_groups(dag, id1, id2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)

	max(depth[ancestor])/max(depth)
}
ADD_GROUP_SIM_METHOD("GroupSim_Ye_2005")


GroupSim_Cho_2007 = function(dag, group1, group2) {

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
ADD_GROUP_SIM_METHOD("GroupSim_Cho_2007")


GroupSim_SimALD = function(dag, group1, group2) {

	n = n_annotations(dag)

	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	i_intersect = intersect(id1, id2)

	if(length(i_intersect) == 0) {
		0
	} else {
		max(1 - n[i_intersect]/sum(n))
	}
}
ADD_GROUP_SIM_METHOD("GroupSim_SimALD")


GroupSim_Jaccard = function(dag, group1, group2, anno_universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, anno_universe, method = "jaccard")
}
ADD_GROUP_SIM_METHOD("GroupSim_Jaccard")


GroupSim_Dice = function(dag, group1, group2, anno_universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, anno_universe, method = "dice")
}
ADD_GROUP_SIM_METHOD("GroupSim_Dice")


GroupSim_Overlap = function(dag, group1, group2, anno_universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, anno_universe, method = "overlap")
}
ADD_GROUP_SIM_METHOD("GroupSim_Overlap")


GroupSim_Kappa = function(dag, group1, group2, anno_universe = NULL) {
	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	l1 = validate_annotated_terms(dag, id1)
	id1 = id1[l1]
	l2 = validate_annotated_terms(dag, id2)
	id2 = id2[l2]
	.sim_overlap_from_two_groups(dag, id1, id2, anno_universe, method = "kappa")
}
ADD_GROUP_SIM_METHOD("GroupSim_Kappa")


.sim_overlap_from_two_groups = function(dag, id1, id2, anno_universe = NULL, method = c("kappa", "jaccard", "dice", "overlap")) {

	check_pkg("proxyC", bioc = FALSE)

	if(!is.null(anno_universe)) {
		anno_universe = which(dag@annotation$names %in% anno_universe)
	} else {
		anno_universe = seq_along(dag@annotation$names)
	}

	if(length(dag@annotation$list) == 0) {
		stop("`annotation` should be set in `create_ontology_DAG()`.")
	}

	mg1 = cpp_get_term_annotations(dag, id1)
	mg1 = mg1[, anno_universe, drop = FALSE]

	mg2 = cpp_get_term_annotations(dag, id2)
	mg2 = mg2[, anno_universe, drop = FALSE]

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
