
#' Import .obo file to an ontology_DAG object
#' 
#' @param file Path of the .obo/.owl file or an URL.
#' 
#' @details Public bio-ontologies can be obtained from [Ontology Foundry](http://obofoundry.org/) or [BioPortal](https://bioportal.bioontology.org/) where
#'    .obo files for the corresponding ontologies can be easily found. If there is no such .obo format file, the .owl file
#'    can be converted to the .obo file using [the **ROBOT** tool](http://robot.obolibrary.org/). An example command is:
#' 
#' ```   
#' robot convert --input ontology.ttl --output ontology.obo --check false
#' ```
#' 
#' @return A list.
#' @export
#' @examples
#' \dontrun{
#' # The plant ontology: http://obofoundry.org/ontology/po.html 
#' dag = import_obo("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo")
#' }
import_obo = function(file) {
	
	ln = readLines(file)

	ind1 = grep("^\\[Term\\]$", ln)
	ind_emptyline = grep("^$", ln)
	ind2 = cpp_match_index(ind1, ind_emptyline)
	n = length(ind1)

	if(n == 0) {
		stop("Cannot find any [Term].")
	}

	lt_terms = vector("list", n)
	for(i in seq_len(n)) {
		if(i %% 1000 == 0) {
			message(strrep("\b", 100), "parsing [Term] sections in the obo file [", i, "/", n, "]", appendLF = FALSE)
		}
		lt_terms[[i]] = process_obo_stanza(ln[seq(ind1[i], ind2[i])])
	}
	message(strrep("\b", 100), "parsing [Term] sections in the obo file [", n, "/", n, "]", appendLF = TRUE)

	## terms
	lt = .wrap_relations(lt_terms)
	term_meta = lt$meta
	term_relations = lt$relations
	

	ind1 = grep("^\\[Typedef\\]$", ln)
	ind2 = cpp_match_index(ind1, ind_emptyline)
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
		lt = .wrap_relations(lt_relations)
		relation_meta = lt$meta
		relation_relations = lt$relations
	} else {
		relation_meta = NULL
		relation_relations = NULL
	}

	## some meta for the whole ontology
	version = gsub("data-version: ", "", grep("^data-version:", ln[1:ind_emptyline[1]], value = TRUE))
	ontology = gsub("ontology: ", "", grep("^ontology:", ln[1:ind_emptyline[1]], value = TRUE))

	list(ontology = ontology, version = version,
		 term_meta = term_meta, term_relations = term_relations, 
		 relation_meta = relation_meta, relation_relations = relation_relations)
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

	# for(i in seq_len(n_terms)) {
	# 	if(i %% 1000 == 0) {
	# 		message(strrep("\b", 100), "validating all terms [", i, "/", n_terms, "]", appendLF = FALSE)
	# 	}

	# 	lt = lt_data[[i]]

	# 	if(length(lt$relationship)) {
	# 		ind = which(lt$relationship %in% all_terms)
	# 		if(length(ind) == 0) {
	# 			lt$relationship = NULL
	# 		} else {
	# 			lt$relationship = lt$relationship[ind]
	# 		}
	# 	}
	# 	lt_data[[i]] = lt
	# }
	# message(strrep("\b", 100), "validating all terms [", n_terms, "/", n_terms, "]", appendLF = TRUE)
	
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

.wrap_relations = function(lt_data) {
	## terms
	lt_data = .validate_relations(lt_data)
	all_elements = sapply(lt_data, "[[", "id")

	meta = data.frame(
		id = all_elements,
		short_id = sapply(lt_data, "[[", "short_id"),
		name = sapply(lt_data, "[[", "name"),
		namespace = sapply(lt_data, "[[", "namespace"),
		def = sapply(lt_data, "[[", "def")
	)
	l_col = sapply(meta, function(x) all(is.na(x)))
	meta = meta[, !l_col, drop = FALSE]


	rl = lapply(lt_data, "[[", "relationship")
	nr = sapply(rl, length)
	child = rep(all_elements, times = nr)
	parent = unlist(rl)
	
	relations = data.frame(child = child, parent = unname(parent), relation = names(parent))
	
	list(meta = meta, relations = relations)
}

process_obo_stanza = function(ln) {

	lt = list()
	i = grep("^id:", ln)
	lt$id = gsub("^id: (\\S+)(\\s*.*)$", "\\1", ln[i])
	lt$short_id = NA_character_

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
		lt$def = gsub('^def: "(.*)".*$', "\\1", ln[i])[1]
	} else {
		i = grep("^property_value: definition", ln)
		if(length(i)) {
			lt$def = gsub('^property_value: definition "(.*)".*$', "\\1", ln[i])[1]
		} else {
			lt$def = NA
		}
	}

	lt$relationship = character(0)
	i = grep("^is_a:", ln)
	if(length(i)) {
		lt$relationship = lt$relationship = c(lt$relationship, structure(gsub("^is_a: (\\S+)(\\s*.*)$", "\\1", ln[i]), names = rep("isa", length(i))))
	}
	
	i = grep("^relationship:", ln)
	if(length(i)) {
		rl = strsplit(ln[i], " ")
		rl_type = sapply(rl, "[[", 2)
		rl_term = sapply(rl, "[[", 3)

		lt$relationship = c(lt$relationship, structure(rl_term, names = rl_type))
	}
	
	i = grep("^is_obsolete:", ln)
	if(length(i)) {
		lt$is_obsolete = gsub("^is_obsolete: (.*)$", "\\1", ln[i])
	} else {
		lt$is_obsolete = "false"
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
import_owl = function(file) {
	
	owl = read_xml(file)

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
	isa = .owl_get_text(Class, ".//rdfs:subClassOf[@rdf:resource]/@rdf:resource", character(0), return_list = TRUE)
	isa = lapply(isa, function(x) {
		structure(x, names = rep("isa", length(x)))
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
		rlp = property[[i]]
		l = !is.na(rlv)
		rlv = rlv[l]
		rlp = rlp[l]

		lt_terms[[i]]$relationship = c(isa[[i]], structure(rlv, names = rlp))
	}


	####### relation / ObjectProperty ########
	message("parsing <owl:ObjectProperty> ...")
	ObjectProperty = xml_find_all(owl, ".//owl:ObjectProperty")
	
	id = xml_attr(ObjectProperty, "about")
	short_id = .owl_get_text(ObjectProperty, ".//*[local-name()='id']", NA)
	short_id = ifelse(is.na(short_id), gsub("^.*#", "", basename(id)), short_id)
	name = .owl_get_text(ObjectProperty, ".//rdfs:label[@xml:lang='en'] | .//rdfs:label[not(@xml:lang)]", NA)
	def = .owl_get_text(ObjectProperty, ".//*[local-name()='IAO_0000115']", NA)
	namespace = .owl_get_text(ObjectProperty, ".//*[local-name()='hasOBONamespace']", NA)
	is_obsolete = .owl_get_text(ObjectProperty, ".//owl:deprecated", "false")
	isa = .owl_get_text(ObjectProperty, ".//rdfs:subPropertyOf[@rdf:resource]/@rdf:resource", character(0), return_list = TRUE)
	isa = lapply(isa, function(x) {
		structure(x, names = rep("isa", length(x)))
	})
	
	lt_relations = vector("list", length(id))
	for(i in seq_along(id)) {
		lt_relations[[i]]$id = id[i]
		lt_relations[[i]]$short_id = short_id[i]
		lt_relations[[i]]$name = name[i]
		lt_relations[[i]]$def = def[i]
		lt_relations[[i]]$namespace = namespace[i]
		lt_relations[[i]]$is_obsolete = is_obsolete[i]
		lt_relations[[i]]$relationship = isa[[i]]
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
	lt = .wrap_relations(lt_terms)
	term_meta = lt$meta
	term_relations = lt$relations

	## relations
	if(length(lt_relations)) {
		lt = .wrap_relations(lt_relations)
		relation_meta = lt$meta
		relation_relations = lt$relations
	} else {
		relation_meta = NULL
		relation_relations = NULL
	}

	## some meta for the whole ontology
	version = xml_text(xml_find_all(owl, ".//owl:Ontology/owl:versionInfo"))
	if(length(version) == 0) {
		version = xml_attr(xml_find_all(owl, ".//owl:Ontology/owl:versionIRI"), "resource")
	}
	ontology = xml_text(xml_find_all(owl, ".//owl:Ontology/*[local-name()='title']"))
	if(length(ontology) == 0) {
		ontology = xml_attr(xml_find_all(owl, ".//owl:Ontology"), "about")
	}

	list(ontology = ontology, version = version,
		 term_meta = term_meta, term_relations = term_relations, 
		 relation_meta = relation_meta, relation_relations = relation_relations)
	
}

