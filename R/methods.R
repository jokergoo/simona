
check_cyclic_node = function(root, nv, lt_children) {
	
	e$visited = rep(FALSE, nv)
	e$cyclic = integer(0)

	lt_children = dag@lt_children
	transverse = function(x) {
		if(length(lt_children[[x]]) == 0) {
			return(NULL)
		}
		if(length(e$cyclic)) {
			return(NULL)
		}
		if(e$visited[x]) {
			e$cyclic = x
			return(NULL)
		}
		e$visited[x] = TRUE
		
		for(t in lt_children[[x]]) {
			transverse(t)
			if(length(e$cyclic) > 0) {
				return(NULL)
			}
		}
	}

	transverse(root)

	e$cyclic

}