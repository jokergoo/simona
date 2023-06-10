
GroupSim_pairwise = function(dag, group1, group2, sim_method, group_sim_method = "avg") {
	group = union(group1, group2)
	sim = term_sim(dag, group, sim_method)
	sim = sim[group1, group2, drop = FALSE]

	# avg
	if(group_sim_method == "avg") {
		mean(sim)
	} else if(group_sim_method == "max") {
		# max
		max(sim)
	} else if(group_sim_method == "bma") {
		#bma
		(rowMeans(sim) + colMeans(sim))/2
	} else if(group_sim_method == "bmm") {
		# bmm
		max(rowMeans(sim), colMeans(sim))
	} else if(group_sim_method == "abm") {
		# abm
		(rowSums(sim) + colSums(sim))/(nrow(sim) + ncol(sim))
	} else if(group_sim_method == "hdf") {
		# hdf
		1 - max(1 - rowMins(sim), 1 - colMins(sim))
	} else if(group_sim_method == "mhdf") {
		# MHDF
		1 - max(1 - rowMeans(sim), 1 - colMeans(sim))
	} else if(group_sim_method == "vhdf") {
		D = (1-sim)^2
		(sqrt(rowMeans(D)) + sqrt(colMeans(D)))/2
	} else if(group_sim_method == "froehlich_2007") {
		exp(-(min(rowMins(sim), colMins(sim))))
	} else if(group_sim_method == "joeng_2014") {
		S = sim^2
		(sqrt(rowMeans(S)) + sqrt(colMeans(S)))/2
	} else {
		stop("method not supported.")
	}
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_pairwise")

GroupSim_SimALN = function(dag, group1, group2) {
	group = union(group1, group2)
	d = cpp_longest_distances_via_LCA(dag, group)
	d = d[group1, group2, drop = FALSE]
	exp(-mean(d))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimALN")


GroupSim_SimINT = function(dag, group1, group2) {
	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_union = union(i1, i2)
	i_intersect = intersect(i1, i2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[i1]^2)/sum(ic_anno[i2]^2)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimINT")


GroupSim_spgk = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	d1 = cpp_shortest_distances_directed(dag, ancestor1)
	d2 = cpp_shortest_distances_directed(dag, ancestor2)

	common_ancestor = intersect(ancestor1, ancestor2)

	if(length(common_ancestor) == 0) {
		0
	} else {
		sum(pmax(0, 2 - abs(d1[common_ancestor, common_ancestor, drop = FALSE] - d2[common_ancestor, common_ancestor, drop = FALSE])))
	}
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_spgk")


GROUPSIM_ANCESTOR_UNION = 1
GROUPSIM_ANCESTOR_INTERSECT = 2

GroupSim_SimGIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_UNION, TRUE)

	sum(ic[ancestor1])/sum(ic[ancestor2])
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimGIC")


GroupSim_SimDIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	2*sum(ic[ancestor1])/(sum(ic[ancestor2]) + sum(ic[ancestor3]))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimDIC")


GroupSim_SimUIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	sum(ic[ancestor1])/max(sum(ic[ancestor2]), sum(ic[ancestor3]))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUIC")


GroupSim_SimUI = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_UNION, TRUE)

	length(ancestor1)/length(ancestor2)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUI")


GroupSim_SimDB = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	2*length(ancestor1)/(length(ancestor2) + length(ancestor3))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimDB")


GroupSim_SimUB = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	length(ancestor1)/max(length(ancestor2), length(ancestor3))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUB")


GroupSim_SimNTO = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1, TRUE)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2, TRUE)

	length(ancestor1)/min(length(ancestor2), length(ancestor3))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimNTO")


GroupSim_SimCOU = function(dag, group1, group2) {
	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[i1]^2)/sum(ic_anno[i2]^2)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimCOU")


GroupSim_SimCOT = function(dag, group1, group2) {
	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/( sum(ic_anno[i1]^2)/sum(ic_anno[i2]^2) - sum(ic_anno[i_intersect]^2))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimCOT")


GroupSim_SimLP = function(dag, group1, group2) {
	depth = dag_depth(dag)
	ancestor = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)
	max(depth[ancestor])
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimLP")


GroupSim_Ye_2005 = function(dag, group1, group2) {
	depth = dag_depth(dag)
	ancestor = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT, TRUE)

	max(depth[ancestor])/max(depth)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Ye_2005")


GroupSim_Cho_2007 = function(dag, group1, group2) {

	n = n_annotations(dag)

	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	log(min(n[i_intersect])/max(n))/log(min(n)/max(n))

}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Cho_2007")


GroupSim_SimALD = function(dag, group1, group2) {

	n = n_annotations(dag)

	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	max(1 - n[i_intersect]/sum(n))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimALD")


GroupSim_Jaccard = function(dag, group1, group2, universe = NULL) {
	.sim_overlap_from_two_groups(dag, group1, group2, universe, method = "jaccard")
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Jaccard")


GroupSim_Dice = function(dag, group1, group2, universe = NULL) {
	.sim_overlap_from_two_groups(dag, group1, group2, universe, method = "dice")
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Dice")


GroupSim_Overlap = function(dag, group1, group2, universe = NULL) {
	.sim_overlap_from_two_groups(dag, group1, group2, universe, method = "overlap")
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Overlap")


GroupSim_Kappa = function(dag, group1, group2, universe = NULL) {
	.sim_overlap_from_two_groups(dag, group1, group2, universe, method = "kappa")
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Kappa")


.sim_overlap_from_two_groups = function(dag, group1, group2, universe = NULL, method = c("kappa", "jaccard", "dice", "overlap")) {

	check_pkg("proxyC", bioc = FALSE)

	id1 = term_to_node_id(dag, group1, strict = FALSE)
	id2 = term_to_node_id(dag, group2, strict = FALSE)

	if(!is.null(universe)) {
		universe = term_to_node_id(dag, universe, strict = FALSE)
		id = intersect(id, universe)
	} else {
		universe = seq_along(dag@annotation$names)
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
