
getAnnotatedTerms = function(dag, anno, in_labels = TRUE) {
	ind = which(dag@annotation$names == anno)
	l = sapply(dag@annotation$list, function(x) ind %in% x)

	if(in_labels) {
		names(dag@terms)[l]
	} else {
		which(l)
	}
}

GroupSim_pairwise = function(group1, group2, sim_method, group_sim_method = "avg") {
	group = union(group1, group2)
	sim = get_Sim_method(sim_method)(group)

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

GroupSim_pairwise_edge = function(dag, group1, group2) {
	d = shortest_path_length(dag, group1, group2)
	exp(-mean(d))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_pairwise_edge")


GroupSim_SimINT = function(dag, group1, group2) {
	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_union = union(i1, i2)
	i_intersect = intersect(i1, i2)

	ic_anno = IC_annotation(dag, use_cache = TRUE)

	sum(ic_anno[i_intersect]^2)/sum(ic_anno[i1]^2)/sum(ic_anno[i2]^2)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimINT")


GROUPSIM_ANCESTOR_UNION = 1
GROUPSIM_ANCESTOR_INTERSECT = 2

GroupSim_SimGIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_UNION)

	sum(ic[ancestor1])/sum(ic[ancestor2])
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimGIC")


GroupSim_SimDIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2)

	2*sum(ic[ancestor1])/(sum(ic[ancestor2]) + sum(ic[ancestor3]))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimDIC")


GroupSim_SimUIC = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2)

	sum(ic[ancestor1])/max(sum(ic[ancestor2]), sum(ic[ancestor3]))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUIC")


GroupSim_SimUI = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_UNION)

	length(ancestor1)/length(ancestor2)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUI")


GroupSim_SimDB = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2)

	2*length(ancestor1)/(length(ancestor2) + length(ancestor3))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimDB")


GroupSim_SimUB = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2)

	length(ancestor1)/max(length(ancestor2), length(ancestor3))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimUB")


GroupSim_SimNTO = function(dag, group1, group2) {
	ancestor1 = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	ancestor2 = cpp_ancestor_of_a_group(dag, group1)
	ancestor3 = cpp_ancestor_of_a_group(dag, group2)

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
	ancestor = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)
	max(depth[ancestor])
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimLP")


GroupSim_Ye_2005 = function(dag, group1, group2) {
	depth = dag_depth(dag)
	ancestor = cpp_ancestor_of_two_groups(dag, group1, group2, GROUPSIM_ANCESTOR_INTERSECT)

	max(depth[ancestor])/max(depth)
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Ye_2005")


GroupSim_Cho_2007 = function(dag, group1, group2) {

	n = get_annotated_terms(dag)

	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	log(min(n[i_intersect])/max(n))/log(min(n)/max(n))

}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Cho_2007")


GroupSim_SimALD = function(dag, group1, group2) {

	n = get_annotated_terms(dag)

	i1 = which(dag@terms == group1)
	i2 = which(dag@terms == group2)

	i_intersect = intersect(i1, i2)

	max(1 - n[i_intersect]/sum(n))
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_SimALD")


GroupSim_Jaccard = function(dag, group1, group2) {
	.term_similarity(list(group1, group2), method = "jaccard")[1, 1]
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Jaccard")


GroupSim_Dice = function(dag, group1, group2) {
	.term_similarity(list(group1, group2), method = "dice")[1, 1]
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Dice")


GroupSim_Overlap = function(dag, group1, group2) {
	.term_similarity(list(group1, group2), method = "overlap")[1, 1]
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Overlap")


GroupSim_Kappa = function(dag, group1, group2) {
	.term_similarity(list(group1, group2), method = "kappa", all = dag@terms)[1, 1]
}
ALL_GROUP_SIM_METHODS = c(ALL_GROUP_SIM_METHODS, "GroupSim_Kappa")


.term_similarity = function (gl, method = c("kappa", "jaccard", "dice", "overlap"), all = NULL) {
    if (is.null(all)) {
        all = unique(unlist(gl))
    }
    else {
        gl = lapply(gl, intersect, all)
    }
    gl = lapply(gl, function(x) as.numeric(factor(x, levels = all)))
    n = length(gl)
    mg = matrix(0, ncol = length(all), nrow = n)
    for (i in seq_len(n)) {
        mg[i, gl[[i]]] = 1
    }
    mg = as(mg, "sparseMatrix")
    method = match.arg(method)[1]
    if (method == "kappa") {
        mat = kappa_dist(mg)
    }
    else if (method == "overlap") {
        mat = overlap_dist(mg)
    }
    else {
        mat = proxyC::simil(mg, method = method)
    }
    mat = as.matrix(mat)
    diag(mat) = 1
    rownames(mat) = colnames(mat) = names(gl)
    return(mat)
}
