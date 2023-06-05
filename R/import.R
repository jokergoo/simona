
#' Import .obo file to an ontology_DAG object
#' 
#' @param file Path of the file or an URL. Pass to [base::readLines()].
#' 
#' @export
import_obo = function(file) {
	ln = readLines(file)
	ind1 = grep("^\\[Term\\]$", ln)
	ind2 = grep("^$", ln)
	ind2 = cpp_match_index(ind1, ind2)
	n = length(ind1)

	id = character(n)
	lt = vector("list", n)
	for(i in seq_len(n)) {
		parsed = .process_obo_term(ln[seq(ind1[i], ind2[i])])
		if(length(parsed)) {
			id[i] = parsed$id
			lt[[i]] = parsed$relations
		}
	}

	children = rep(id, sapply(lt, length))
	parents = unlist(lt)
	relations = unlist(lapply(lt, names))

	createOntologyDAG(parents = parents, children = children, relations = relations)
}

.process_obo_term = function(ln) {

	i0 = grep("^is_obsolete:", ln)
	if(length(i0)) {
		return(list())
	}

	i1 = grep("^id:", ln)
	id = gsub("^id: (.*)(\\s*.*)$", "\\1", ln[i1])

	i2 = grep("^is_a:", ln)
	is_a = gsub("^is_a: (.*) !.*$", "\\1", ln[i2])
	names(is_a) = rep("is_a", length(is_a))
	
	i3 = grep("^relationship: part_of", ln)
	part_of = gsub("^relationship: part_of (.*) !.*$", "\\1", ln[i3])
	names(part_of) = rep("part_of", length(part_of))

	relations = c(is_a, part_of)
	relations = setdiff(relations, id)
	list(id = id, relations = relations)
}
