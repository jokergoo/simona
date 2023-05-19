
df = toTable(GOMFCHILDREN)

dag = createOntologyDAG(parents = df[, 2], children = df[, 1], relations = df[, 3])


annotation = lapply(as.list(org.Hs.egGO2ALLEGS), unique)
IC_annotation(dag, annotation)

root = dag@root

depth = dag_depth(dag)
lt_children = dag@lt_children
lt_offspring = dag@lt_offspring

term_pos = matrix(nrow = dag@n_terms, ncol = 2)

get_circular_range = function(parent_range, nodes) {
	size = sapply(lt_offspring[nodes], length)
	p = cumsum(size)/sum(size)
	right = p*(parent_range[2] - parent_range[1]) + parent_range[1]
	left = c(0, p[-length(p)])*(parent_range[2] - parent_range[1]) + parent_range[1]
	
	cbind(left = left, right = right)
}

reorder_node = function(node, lt_offspring) {
	lt = lt_offspring[node]
	d = 1 - term_similarity(lt, method = "jaccard")
	node[ hclust(as.dist(d))$order ]
}


parents = root
parent_range = cbind(left = 0, right = 360)
parent_pos = cbind(x = 0, y = 0)

term_pos[parents, ] = c(0, 0)
col = rep("black", dag@n_terms)

level1 = lt_children[[root]]
library(Polychrome)

col[level1] = Polychrome::alphabet.colors(length(level1))

current_depth = 0
while(1) {
	parent_range2 = NULL
	parents2 = NULL
	current_depth = current_depth + 1
	for(ip in seq_along(parents)) {
		children = lt_children[[ parents[ip] ]]
		children = children[depth[children] == current_depth]

		if(length(children)) {
			# reorder children
			if(length(children) > 2) {
				children = reorder_node(children, lt_offspring)
			}

			df = get_circular_range(parent_range[ip, ], children)

			mid = (df[, 1] + df[, 2])/2
			pos = circlize:::polar2Cartesian(cbind(mid, current_depth))
			term_pos[children, ] = pos

			parent_range2 = rbind(parent_range2, df)
			parents2 = c(parents2, children)

			if(current_depth > 1) {
				col[children] = col[parents[ip]]
			}
		}
	}

	parents = parents2
	parent_range = parent_range2

	if(length(parents2) == 0) {
		break
	}
}


x1 = x2 = y1 = y2 = numeric(0)
edge_id = integer(0)
for(i in seq_along(lt_children)) {
	children = lt_children[[i]]
	if(length(children)) {
		x1 = c(x1, rep(term_pos[i, 1], length(children)))
		y1 = c(y1, rep(term_pos[i, 2], length(children)))
		x2 = c(x2, term_pos[children, 1])
		y2 = c(y2, term_pos[children, 2])
		for(c in children) {
			edge_id = c(edge_id, get.edge.ids(dag@graph, c(i, c)))
		}
	}
}
type = E(dag@graph)$relation[edge_id]

max_depth = max(depth) + 1
grid.newpage()
pushViewport(viewport(x = unit(0, "npc"), just = "left", width = unit(1, "snpc"), height = unit(1, "snpc"), 
	xscale = c(-max_depth, max_depth), yscale = c(-max_depth, max_depth)))
grid.points(term_pos[, 1], term_pos[, 2], default.units = "native", pch = 16, size = unit(4, "pt"), 
	gp = gpar(col = add_transparency(col, 0.5)))

grid.segments(x1, y1, x2, y2, default.units = "native", 
	gp = gpar(col = ifelse(type == "isa", "#FF000010", ifelse(type == "part of", "#0000FF20", "#00000020"))))

level1_go = dag@terms[level1]
level1_go = paste0(level1_go, "~", Term(GOTERM[level1_go]))
level1_go = ifelse(nchar(level1_go) > 50, paste0(substr(level1_go, 0, 50), "..."), level1_go)
lgd = packLegend(
	Legend(title = "Terms", labels = level1_go, 
		type = "points", legend_gp = gpar(col = col[level1])),
	Legend(title = "Relations", labels = c("isa", "part of", "regulate"),
		type = "lines", legend_gp = gpar(col = c("red", "blue", "grey")))
)

draw(lgd, x = unit(1, "npc") - unit(1.5, "cm"), y = unit(0.5, "npc"), just = "left")

grid.text("GO BP ontology", x = unit(1, "npc") - unit(1.5, "cm"), y = unit(0.92, "npc"), 
	just = "left", gp = gpar(fontsize = 20))





