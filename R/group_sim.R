
GroupSim_pairwise = function(group1, group2) {
	group = union(group1, group2)
	sim = ...

	sim = sim[group1, group2, drop = FALSE]

	# avg
	mean(sim)

	# max
	max(sim)

	#bma
	(rowMeans(sim) + colMeans(sim))/2

	# bmm
	max(rowMeans(sim), colMeans(sim))

	# abm
	(rowSums(sim) + colSums(sim))/(nrow(sim) + ncol(sim))

	# hdf
	1 - max(1 - rowMins(sim), 1 - colMins(sim))

	# MHDF
	max(1 - rowMeans(sim), 1)
}

GroupSim_pairwise_edge = function(group1, group2) {
	group = union(group1, group2)
	d = distances(dag@graph, group1, group2, mode = "all")
	exp(-mean(d))
}

GroupSim_IC_based = function() {
	lt_ancestor = dag@lt_ancestor
	an1 = unique(unlist(lt_ancestor[group1]))
	an2 = unique(unlist(lt_ancestor[group2]))
	set1 = intersect(an1, an2)
	set2 = union(an1, an2)

	sum(ic[set1])/sum(ic[set2])
}
