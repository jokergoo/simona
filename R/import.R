
	      
#' Import .obo file to an ontology_DAG object
#' 
#' @param file Path of the ontology file or an URL.
#' @param relation_type Semantic relation types to include. Note `is_a` relation is always included.
#' 
#' @details 
#' 
#' Public bio-ontologies can be obtained from [Ontology Foundry](http://obofoundry.org/) or [BioPortal](https://bioportal.bioontology.org/). 
#' 
#' The `import_obo()` function parses the ontology file in `.obo` format. To parse other formats, external tool `robot.jar` is required.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @examples
#' \dontrun{
#' # The plant ontology: http://obofoundry.org/ontology/po.html 
#' dag = import_obo("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo")
#' }
import_obo = function(file, relation_type = "part_of") {
	
	if(grepl("^(http|https|ftp)://.*\\.gz$", file)) {
		ln = readLines(gzcon(url(file)))
	} else if(grepl("\\.gz$", file)) {
		ln = readLines(gzfile(file))
	} else if(grepl("^(http|https|ftp)://", file)) {
		ln = readLines(url(file))
	} else {
		ln = readLines(file)
	}

	ind_stanza = c(grep("^\\[.*\\]$", ln), length(ln))

	### relations ###
	ind1 = grep("^\\[Typedef\\]$", ln)
	ind2 = cpp_match_index(ind1, ind_stanza[-1]-1)
	n = length(ind1)

	if(n) {

		lt_relations = vector("list", n)
		for(i in seq_len(n)) {
			if(i %% 1000 == 0) {
				message(strrep("\b", 100), "parsing [Typedef] sections in the obo file [", i, "/", n, "]", appendLF = FALSE)
			}
			lt_relations[[i]] = process_obo_stanza(ln[seq(ind1[i], ind2[i])])
		}
		message(strrep("\b", 100), "parsing [Typedef] sections in the obo file [", n, "/", n, "]", appendLF = TRUE)

		## relations
		lt = .wrap_relations(lt_relations, type = "relation")
		relation_meta = lt$meta
		relation_relations = lt$relations
	} else {
		relation_meta = data.frame()
		relation_relations = data.frame()
	}

	relations_id_to_name = structure(ifelse(is.na(relation_meta$name), relation_meta$id, relation_meta$name), names = relation_meta$id)
	relations_id_to_name["is_a"] = "is_a"
	
	if(nrow(relation_relations)) {
		
		relation_relations$parent = unname(relations_id_to_name[relation_relations$parent])
		relation_relations$child = unname(relations_id_to_name[relation_relations$child])

		suppressWarnings(suppressMessages(relations_DAG <- create_ontology_DAG(relation_relations$parent, relation_relations$child)))

		relation_type = merge_offspring_relation_types(relations_DAG, relation_type)
	} else {
		relations_DAG = NULL
	}


	### term ###
	ind1 = grep("^\\[Term\\]$", ln)
	ind2 = cpp_match_index(ind1, ind_stanza[-1]-1)
	n = length(ind1)

	if(n == 0) {
		stop("Cannot find any [Term].")
	}

	lt_terms = vector("list", n)
	for(i in seq_len(n)) {
		if(i %% 1000 == 0) {
			message(strrep("\b", 100), "parsing [Term] sections in the obo file [", i, "/", n, "]", appendLF = FALSE)
		}
		lt_terms[[i]] = process_obo_stanza(ln[seq(ind1[i], ind2[i])], relation_type)
	}
	message(strrep("\b", 100), "parsing [Term] sections in the obo file [", n, "/", n, "]", appendLF = TRUE)

	lt = .wrap_relations(lt_terms, type = "term")
	term_meta = lt$meta
	term_relations = lt$relations
	
	term_relations$relation = unname(relations_id_to_name[term_relations$relation])


	## some meta for the whole ontology
	for(i in seq_along(ln)) {
		if(grepl("^\\[", ln[i])) {
			break
		}
	}
	version = gsub("data-version: ", "", grep("^data-version:", ln[1:(i-1)], value = TRUE))
	ontology = gsub("ontology: ", "", grep("^ontology:", ln[1:(i-1)], value = TRUE))
	default_namespace = gsub("default-namespace: ", "", grep("^default-namespace:", ln[1:(i-1)], value = TRUE))

	if(length(default_namespace)) {
		term_meta$namespace[is.na(term_meta$namespace)] = default_namespace
	}

	dag = create_ontology_DAG(parents = term_relations$parent, children = term_relations$child, relations = term_relations$relation,
		source = paste0(ontology, ", ", version), relations_DAG = relations_DAG)

	rownames(term_meta) = term_meta$id
	term_meta = term_meta[dag@terms, , drop = FALSE]
	
	if(dag_root(dag) == "_all_") {
		nr = nrow(term_meta)
		term_meta$id[nr] = "_all_"
		term_meta$short_id[nr] = "_all_"
	}

	mcols(dag) = term_meta

	if(!any(duplicated(term_meta$short_id))) {
		dag@terms = term_meta$short_id
		rownames(dag@elementMetadata) = dag@terms
	}

	dag
}

# validate relations
.validate_relations = function(lt_data) {
	is_obsolete = sapply(lt_data, "[[", "is_obsolete") == "true"
	if(any(is_obsolete)) {
		message("remove ", sum(is_obsolete), " obsolete terms")
	}
	lt_data = lt_data[!is_obsolete]

	all_terms = sapply(lt_data, "[[", "id")
	n_terms = length(all_terms)
	
	lt = lapply(lt_data, "[[", "relationship")

	lt2 = intersectToList_logical(lt, all_terms)
	for(i in seq_len(n_terms)) {
		if(any(lt2[[i]])) {
			lt_data[[i]]$relationship = lt_data[[i]]$relationship[ lt2[[i]] ]
		} else {
			lt_data[[i]]$relationship = character(0)
		}
	}

	lt_data
}

.wrap_relations = function(lt_data, type = "term") {
	## terms
	lt_data = .validate_relations(lt_data)
	all_elements = sapply(lt_data, "[[", "id")

	if(type == "term") {
		meta = data.frame(
			id = all_elements,
			short_id = sapply(lt_data, "[[", "short_id"),
			name = sapply(lt_data, "[[", "name"),
			namespace = sapply(lt_data, "[[", "namespace"),
			definition = sapply(lt_data, "[[", "def")
		)
	} else {
		meta = data.frame(
			id = all_elements,
			short_id = sapply(lt_data, "[[", "short_id"),
			name = sapply(lt_data, "[[", "name"),
			namespace = sapply(lt_data, "[[", "namespace"),
			definition = sapply(lt_data, "[[", "def")
		)
	}
	
	rl = lapply(lt_data, "[[", "relationship")
	nr = sapply(rl, length)
	child = rep(all_elements, times = nr)
	parent = unlist(rl)
	
	relations = data.frame(child = child, parent = unname(parent), relation = names(parent))
	relations = relations[relations$child != relations$parent, , drop = FALSE]
	
	list(meta = meta, relations = relations)
}

process_obo_stanza = function(ln, relation_type = "part_of") {

	lt = list()
	i = grep("^id:", ln)
	lt$id = gsub("^id: (\\S+)(\\s*.*)$", "\\1", ln[i])[1]
	lt$short_id = basename(lt$id)

	i = grep("^name:", ln)
	if(length(i)) {
		lt$name = gsub("^name: (.*)$", "\\1", ln[i])[1]
	} else {
		i = grep("^property_value: prefLabel", ln)
		if(length(i)) {
			lt$name = gsub('^property_value: prefLabel "(.*)".*$', "\\1", ln[i])[1]
		} else {
			lt$name = NA
		}
	}

	i = grep("^namespace:", ln)
	if(length(i)) {
		lt$namespace = gsub("^namespace: (.*)$", "\\1", ln[i])[1]
	} else {
		lt$namespace = NA
	}

	i = grep("^def:", ln)
	if(length(i)) {
		lt$def = gsub('^def: "(.*)" ?\\[.*\\].*$', "\\1", ln[i])[1]
	} else {
		i = grep("^property_value: definition", ln)
		if(length(i)) {
			lt$def = gsub('^property_value: definition "(.*)" ?\\[.*\\].*$', "\\1", ln[i])[1]
		} else {
			lt$def = NA
		}
	}

	lt$relationship = character(0)
	i = grep("^is_a:", ln)
	if(length(i)) {
		l = grepl("^is_a: (\\S+)\\s+\\{.*\\}\\s*", ln[i])
		i = i[!l]
		if(length(i)) {
			is_a = gsub("^is_a: (\\S+)(\\s*.*)$", "\\1", ln[i])
			is_a = unique(is_a)
			lt$relationship = lt$relationship = c(lt$relationship, structure(is_a, names = rep("is_a", length(is_a))))
		}
	}

	i = grep("^relationship:", ln)
	if(length(i)) {

		rl = strsplit(ln[i], " ")
		rl = lapply(rl, function(x) c(x, ""))
		rl_type = sapply(rl, "[[", 2)
		rl_term = sapply(rl, "[[", 3)

		l1 = !grepl("^\\{", sapply(rl, "[[", 4))
		l2 = rl_type %in% relation_type
		l = l1 & l2
		rl_type = rl_type[l]
		rl_term = rl_term[l]

		lt$relationship = c(lt$relationship, structure(rl_term, names = rl_type))
	}

	i = grep("^is_obsolete:", ln)
	if(length(i)) {
		lt$is_obsolete = gsub("^is_obsolete: (.*)$", "\\1", ln[i])
	} else {
		lt$is_obsolete = "false"
	}

	# specific for Typedef
	i = grep("^is_transitive:", ln)
	if(length(i)) {
		lt$is_transitive = gsub("^is_transitive: (.*)$", "\\1", ln[i])[1]
	} else {
		lt$is_transitive = "false"
	}

	i = grep("^inverse_of:", ln)
	if(length(i)) {
		lt$inverse_of = gsub("^inverse_of: (\\S+)(\\s*.*)$", "\\1", ln[i])[1]
	} else {
		lt$inverse_of = ""
	}

	lt
}


.owl_get_text = function(nodes, xpath, default = NA, return_list = FALSE) {
	if(return_list) {
		lapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
			if(length(x) == 0) {
				character(0)
			} else {
				xml_text(x)
			}
		})
	} else {
		sapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
			if(length(x) == 0) {
				default
			} else {
				xml_text(x)[1]
			}
		})
	}
}

.owl_get_attr = function(nodes, xpath, attr, default = NA) {
	lapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
		xml_attr(x, attr)
	})
}

#' @rdname import_obo
#' @export
#' @import xml2
#' @export
import_owl = function(file, relation_type = "part_of") {
	
	owl = read_xml(file, options = "HUGE")

	####### relation / ObjectProperty ########
	message("parsing <owl:ObjectProperty> ...")
	ObjectProperty = xml_find_all(owl, ".//owl:ObjectProperty")
	
	id = xml_attr(ObjectProperty, "about")
	short_id = .owl_get_text(ObjectProperty, ".//*[local-name()='id']", NA)
	short_id = ifelse(is.na(short_id), gsub("^.*#", "", basename(id)), short_id)
	name = .owl_get_text(ObjectProperty, ".//rdfs:label[@xml:lang='en'] | .//rdfs:label[not(@xml:lang)]", NA); name = gsub(" ", "_", name);
	def = .owl_get_text(ObjectProperty, ".//*[local-name()='IAO_0000115']", NA)
	namespace = .owl_get_text(ObjectProperty, ".//*[local-name()='hasOBONamespace']", NA)
	is_obsolete = .owl_get_text(ObjectProperty, ".//owl:deprecated", "false")
	is_a = .owl_get_text(ObjectProperty, ".//rdfs:subPropertyOf[@rdf:resource]/@rdf:resource", character(0), return_list = TRUE)
	is_a = lapply(is_a, function(x) {
		structure(x, names = rep("is_a", length(x)))
	})
	
	lt_relations = vector("list", length(id))
	for(i in seq_along(id)) {
		lt_relations[[i]]$id = id[i]
		lt_relations[[i]]$short_id = short_id[i]
		lt_relations[[i]]$name = name[i]
		lt_relations[[i]]$def = def[i]
		lt_relations[[i]]$namespace = namespace[i]
		lt_relations[[i]]$is_obsolete = is_obsolete[i]
		lt_relations[[i]]$relationship = is_a[[i]]
	}

	## relations
	if(length(lt_relations)) {
		lt = .wrap_relations(lt_relations, type = "relation")
		relation_meta = lt$meta
		relation_relations = lt$relations
	} else {
		relation_meta = data.frame()
		relation_relations = data.frame()
	}

	relations_id_to_name = structure(ifelse(is.na(relation_meta$name), relation_meta$id, relation_meta$name), names = relation_meta$id)
	relations_id_to_name["is_a"] = "is_a"
	
	if(nrow(relation_relations)) {
		
		relation_relations$parent = unname(relations_id_to_name[relation_relations$parent])
		relation_relations$child = unname(relations_id_to_name[relation_relations$child])

		suppressWarnings(suppressMessages(relations_DAG <- create_ontology_DAG(relation_relations$parent, relation_relations$child)))

		relation_type = merge_offspring_relation_types(relations_DAG, relation_type)
	} else {
		relations_DAG = NULL
	}


	### Class #####
	Class = xml_find_all(owl, ".//owl:Class")
	id = xml_attr(Class, "about")
	l = !is.na(id)

	id = id[l]
	Class = Class[l]

	n = length(Class)

	if(n == 0) {
		stop("Cannot find any owl:Class.")
	}

	message("parsing <owl:Class> ...")
	id = xml_attr(Class, "about")
	short_id = .owl_get_text(Class, ".//*[local-name()='id']", NA)
	short_id = ifelse(is.na(short_id), gsub("^.*#", "", basename(id)), short_id)
	name = .owl_get_text(Class, ".//rdfs:label[@xml:lang='en'] | .//rdfs:label[not(@xml:lang)]", NA)
	def = .owl_get_text(Class, ".//*[local-name()='IAO_0000115']", NA)
	namespace = .owl_get_text(Class, ".//*[local-name()='hasOBONamespace']", NA)
	is_obsolete = .owl_get_text(Class, ".//owl:deprecated", "false")
	is_a = .owl_get_text(Class, ".//rdfs:subClassOf[@rdf:resource]/@rdf:resource", character(0), return_list = TRUE)
	is_a = lapply(is_a, function(x) {
		structure(x, names = rep("is_a", length(x)))
	})
	value = .owl_get_attr(Class, ".//rdfs:subClassOf/owl:Restriction/*[self::owl:someValuesFrom or self::owl:allValuesFrom or self::owl:onClass]", "resource")
	property = .owl_get_attr(Class, ".//rdfs:subClassOf/owl:Restriction/*[self::owl:someValuesFrom or self::owl:allValuesFrom or self::owl:onClass]/preceding-sibling::owl:onProperty", "resource")

	lt_terms = vector("list", length(id))
	for(i in seq_along(id)) {
		lt_terms[[i]]$id = id[i]
		lt_terms[[i]]$short_id = short_id[i]
		lt_terms[[i]]$name = name[i]
		lt_terms[[i]]$def = def[i]
		lt_terms[[i]]$namespace = namespace[i]
		lt_terms[[i]]$is_obsolete = is_obsolete[i]
		
		rlv = value[[i]]
		rlp = relations_id_to_name[ property[[i]] ]
		l = !is.na(rlv)
		rlv = rlv[l]
		rlp = rlp[l]

		l = rlp %in% relation_type
		rlv = rlv[l]
		rlp = rlp[l]

		lt_terms[[i]]$relationship = c(is_a[[i]], structure(rlv, names = rlp))
	}

	###### Description ######
	Description = xml_find_all(owl, ".//rdf:Description")

	if(length(Description)) {
		message("parsing <rdf:Description> ...")
		id = xml_attr(Description, "about")
		name = .owl_get_text(Description, ".//*[local-name()='prefLabel']", NA)
		def = .owl_get_text(Description, ".//*[local-name()='definition']", NA)

		lt_description = vector("list", length(id))
		for(i in seq_along(id)) {
			lt_description[[i]]$id = id[i]
			lt_description[[i]]$name = name[i]
			lt_description[[i]]$def = def[i]
		}
		names(lt_description) = id
	} else {
		lt_description = list()
	}

	###### fill name and def in `lt_terms` #######
	if(length(lt_description)) {
		lt_terms = lapply(lt_terms, function(x) {
			if(is.na(x$name)) {
				v = lt_description[[x$id]]$name
				if(length(v) > 0) {
					x$name = v
				}
			}
			if(is.na(x$def)) {
				v = lt_description[[x$id]]$def
				if(length(v) > 0) {
					x$def = v
				}
			}
			x
		})

		lt_relations = lapply(lt_relations, function(x) {
			if(is.na(x$name)) {
				v = lt_description[[x$id]]$name
				if(length(v) > 0) {
					x$name = v
				}
			}
			if(is.na(x$def)) {
				v = lt_description[[x$id]]$def
				if(length(v) > 0) {
					x$def = v
				}
			}
			x
		})
	}

	## terms
	lt = .wrap_relations(lt_terms, "term")
	term_meta = lt$meta
	term_relations = lt$relations

	## some meta for the whole ontology
	version = xml_text(xml_find_all(owl, ".//owl:Ontology/owl:versionInfo"))
	if(length(version) == 0) {
		version = xml_attr(xml_find_all(owl, ".//owl:Ontology/owl:versionIRI"), "resource")
	}
	ontology = xml_text(xml_find_all(owl, ".//owl:Ontology/*[local-name()='title']"))
	if(length(ontology) == 0) {
		ontology = xml_attr(xml_find_all(owl, ".//owl:Ontology"), "about")
	}

	dag = create_ontology_DAG(parents = term_relations$parent, children = term_relations$child, relations = term_relations$relation,
		source = paste0(ontology, ", ", version), relations_DAG = relations_DAG)


	rownames(term_meta) = term_meta$id
	term_meta = term_meta[dag@terms, , drop = FALSE]
	
	if(dag_root(dag) == "_all_") {
		nr = nrow(term_meta)
		term_meta$id[nr] = "_all_"
		term_meta$short_id[nr] = "_all_"
	}

	mcols(dag) = term_meta


	if(!any(duplicated(term_meta$short_id))) {
		dag@terms = term_meta$short_id
		rownames(dag@elementMetadata) = dag@terms
	}

	dag
	
}

#' @param robot_jar The path of the `robot.jar` file. It can be downloaded from https://github.com/ontodev/robot/releases.
#'         Internally, the file is converted to the obo format and parsed by `import_obo()`. The value of `robot_jar` can be
#'         set as a global option `simone_opt$robot_jar = ...`.
#' @param JAVA_ARGS Options for `java`. For example you can set `-Xmx20G` if you want to increase the memory to 20G for java.
#' @details
#' `robot.jar` can automatically recognize the following formats:
#' 
#' - `json`: OBO Graphs JSON
#' - `obo`: OBO Format
#' - `ofn`: OWL Functional
#' - `omn`: Manchester
#' - `owl`: RDF/XML
#' - `owx`: OWL/XML
#' - `ttl`: Turtle 
#' 
#' @rdname import_obo
#' @export
#' @importFrom utils download.file
#' @examples
#' \dontrun{
#' # The plant ontology: http://obofoundry.org/ontology/po.html 
#' dag = import_ontology("http://purl.obolibrary.org/obo/po.owl", robot_jar = ...)
#' }
import_ontology = function(file, robot_jar = simone_opt$robot_jar, JAVA_ARGS = "") {

	if(grepl("\\.(obo|obo.gz)$", file, ignore.case = TRUE)) {
		return(import_obo(file))
	}

	# if(grepl("\\.(owl|owl.gz)$", file, ignore.case = TRUE)) {
	# 	oe = try(obj <- import_owl(file), silent = TRUE)
	# 	if(!inherits(oe, "try-error")) {
	# 		return(obj)
	# 	}
	# }

	if(Sys.which("java") == "") {
		stop("Java is not available.")
	}

	if(!file.exists(robot_jar)) {
		stop("Cannot find robot.jar. It can be downloaded from https://github.com/ontodev/robot/releases.")
	}
	robot_jar = normalizePath(robot_jar)

	if(grep("^(http|ftp)", file)) {
		file2 = tempfile(fileext = paste0("_", basename(file)))
		download.file(file, destfile = file2, quiet = TRUE)
		on.exit(file.remove(file2))

		file = file2
	}

	message("converting ", basename(file), " to the obo format.")
	output = tempfile(fileext = ".obo.gz")

	file = normalizePath(file)
	cmd = qq("java @{JAVA_ARGS} -jar '@{robot_jar}' convert --input '@{file}' --format obo --output '@{output}' --check false")
	message("  ", cmd)

	code = system(cmd)
	if(code != 0) {
		stop("There is an error when executing robot.jar.")
	}

	lt = import_obo(output)

	if(file.exists(output)) {
		file.remove(output)
	}

	lt
}
