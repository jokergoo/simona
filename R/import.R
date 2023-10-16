
	      
#' Import ontology file to an ontology_DAG object
#' 
#' @param file Path of the ontology file or an URL.
#' @param relation_type Semantic relation types to include. Note `is_a` relation is always included.
#' @param inherit_relations Relations may also be structured as a DAG. It controls whether to merge with a relations's offspring relations.
#' @param verbose Whether to print messages.
#' @param ... Pass to [`create_ontology_DAG()`].
#' 
#' @details 
#' 
#' Public bio-ontologies can be obtained from [Ontology Foundry](http://obofoundry.org/) or [BioPortal](https://bioportal.bioontology.org/). 
#' 
#' The `import_obo()` function parses the ontology file in `.obo` format. To parse other formats, external tool `robot.jar` is required.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom utils read.csv
#' @examples
#' \donttest{
#' # The plant ontology: http://obofoundry.org/ontology/po.html 
#' import_obo("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo")
#' }
import_obo = function(file, relation_type = character(0), inherit_relations = TRUE, verbose = simona_opt$verbose, ...) {
	
	if(grepl("^(http|https|ftp)://.*\\.gz$", file)) {
		con = url(file)
		con2 = gzcon(con)
		ln = readLines(con2)
		close(con)
		close(con2)
	} else if(grepl("\\.gz$", file)) {
		con2 = gzfile(file)
		ln = readLines(con2)
		close(con2)
	} else if(grepl("^(http|https|ftp)://", file)) {
		con = url(file)
		ln = readLines(con)
		close(con)
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
				if(verbose) message(strrep("\b", 100), "Parsing [Typedef] sections in the obo file [", i, "/", n, "]", appendLF = FALSE)
			}
			lt_relations[[i]] = process_obo_stanza(ln[seq(ind1[i], ind2[i])])
		}
		if(verbose) message(strrep("\b", 100), "Parsing [Typedef] sections in the obo file [", n, "/", n, "]", appendLF = TRUE)

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

		if(inherit_relations) {
			relation_type = merge_offspring_relation_types(relations_DAG, relation_type)
		}
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
			if(verbose) message(strrep("\b", 100), "Parsing [Term] sections in the obo file [", i, "/", n, "]", appendLF = FALSE)
		}
		lt_terms[[i]] = process_obo_stanza(ln[seq(ind1[i], ind2[i])], relation_type)
	}
	if(verbose) message(strrep("\b", 100), "Parsing [Term] sections in the obo file [", n, "/", n, "]", appendLF = TRUE)

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
	version = gsub("data-version: ", "", grep("^data-version:", ln[seq_len(i-1)], value = TRUE))
	ontology = gsub("ontology: ", "", grep("^ontology:", ln[seq_len(i-1)], value = TRUE))
	default_namespace = gsub("default-namespace: ", "", grep("^default-namespace:", ln[seq_len(i-1)], value = TRUE))

	if(length(default_namespace)) {
		term_meta$namespace[is.na(term_meta$namespace)] = default_namespace
	}

	term_meta = term_meta[term_meta$id %in% c(term_relations$parent, term_relations$child), , drop = FALSE]
	internal_id = term_meta$short_id
	dd = internal_id[duplicated(internal_id)]
	ldd = internal_id %in% dd
	internal_id[ldd] = term_meta$id[ldd]
	idmap = structure(internal_id, names = term_meta$id)
	
	dag = create_ontology_DAG(parents = idmap[term_relations$parent], children = idmap[term_relations$child], relations = term_relations$relation,
		source = paste0(ontology, ", ", version), relations_DAG = relations_DAG, verbose = verbose, ...)
	rownames(term_meta) = internal_id

	term_meta = term_meta[dag@terms, , drop = FALSE]
	
	if(dag_root(dag) == SUPER_ROOT) {
		term_meta$id[dag@root] = SUPER_ROOT
		term_meta$short_id[dag@root] = SUPER_ROOT
	}

	mcols(dag) = term_meta

	dag
}

# validate relations
.validate_relations = function(lt_data) {
	is_obsolete = vapply(lt_data, "[[", "is_obsolete", FUN.VALUE = character(1)) == "true"
	if(any(is_obsolete)) {
		message("remove ", sum(is_obsolete), " obsolete terms")
	}
	lt_data = lt_data[!is_obsolete]

	all_terms = vapply(lt_data, "[[", "id", FUN.VALUE = character(1))
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
	all_elements = vapply(lt_data, "[[", "id", FUN.VALUE = character(1))

	if(type == "term") {
		meta = data.frame(
			id = all_elements,
			short_id = vapply(lt_data, "[[", "short_id", FUN.VALUE = character(1)),
			name = vapply(lt_data, "[[", "name", FUN.VALUE = character(1)),
			namespace = vapply(lt_data, "[[", "namespace", FUN.VALUE = character(1)),
			definition = vapply(lt_data, "[[", "def", FUN.VALUE = character(1))
		)
	} else {
		meta = data.frame(
			id = all_elements,
			short_id = vapply(lt_data, "[[", "short_id", FUN.VALUE = character(1)),
			name = vapply(lt_data, "[[", "name", FUN.VALUE = character(1)),
			namespace = vapply(lt_data, "[[", "namespace", FUN.VALUE = character(1)),
			definition = vapply(lt_data, "[[", "def", FUN.VALUE = character(1))
		)
	}
	
	rl = lapply(lt_data, "[[", "relationship")
	nr = vapply(rl, length, FUN.VALUE = integer(1))
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
			lt$name = NA_character_
		}
	}

	i = grep("^namespace:", ln)
	if(length(i)) {
		lt$namespace = gsub("^namespace: (.*)$", "\\1", ln[i])[1]
	} else {
		lt$namespace = NA_character_
	}

	i = grep("^def:", ln)
	if(length(i)) {
		lt$def = gsub('^def: "(.*)" ?\\[.*\\].*$', "\\1", ln[i])[1]
	} else {
		i = grep("^property_value: definition", ln)
		if(length(i)) {
			lt$def = gsub('^property_value: definition "(.*)" ?\\[.*\\].*$', "\\1", ln[i])[1]
		} else {
			lt$def = NA_character_
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
		rl_type = vapply(rl, "[[", 2, FUN.VALUE = character(1))
		rl_term = vapply(rl, "[[", 3, FUN.VALUE = character(1))

		l1 = !grepl("^\\{", vapply(rl, "[[", 4, FUN.VALUE = character(1)))
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


.owl_get_text = function(nodes, xpath, default = NA_character_, return_list = FALSE) {
	if(return_list) {
		lapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
			if(length(x) == 0) {
				character(0)
			} else {
				xml_text(x)
			}
		})
	} else {
		vapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
			if(length(x) == 0) {
				default
			} else {
				xml_text(x)[1]
			}
		}, FUN.VALUE = character(1))
	}
}

.owl_get_attr = function(nodes, xpath, attr, default = NA_character_) {
	lapply(xml_find_all(nodes, xpath, flatten = FALSE), function(x) {
		xml_attr(x, attr)
	})
}

#' @rdname import_obo
#' @details `import_owl()` only recognizes `<owl:Class>` and `<owl:ObjectProperty>`. If the .owl file does not contain these tags,
#'     please use `import_ontology()` directly.
#' @export
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_text
#' @export
#' @examples
#' \donttest{
#' import_owl("http://purl.obolibrary.org/obo/po.owl") 
#' }
import_owl = function(file, relation_type = character(0), inherit_relations = TRUE, verbose = simona_opt$verbose, ...) {
	
	owl = read_xml(file, options = "HUGE")

	####### relation / ObjectProperty ########
	ObjectProperty = xml_find_all(owl, ".//owl:ObjectProperty")
	
	if(verbose) message("Parsing ", length(ObjectProperty), " <owl:ObjectProperty> ...")
	id = xml_attr(ObjectProperty, "about")
	short_id = .owl_get_text(ObjectProperty, ".//*[local-name()='id']", NA_character_)
	short_id = ifelse(is.na(short_id), gsub("^.*#", "", basename(id)), short_id)
	name = .owl_get_text(ObjectProperty, ".//rdfs:label[@xml:lang='en'] | .//rdfs:label[not(@xml:lang)]", NA_character_); name = gsub(" ", "_", name);
	def = .owl_get_text(ObjectProperty, ".//*[local-name()='IAO_0000115']", NA_character_)
	namespace = .owl_get_text(ObjectProperty, ".//*[local-name()='hasOBONamespace']", NA_character_)
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

	relations_id_to_name = structure(ifelse(is.na(relation_meta$name), relation_meta$short_id, relation_meta$name), names = relation_meta$id)
	relations_id_to_name[is.na(relations_id_to_name)] = relation_meta$id[is.na(relations_id_to_name)]
	relations_id_to_name["is_a"] = "is_a"

	if(nrow(relation_relations)) {
		
		relation_relations$parent = unname(relations_id_to_name[relation_relations$parent])
		relation_relations$child = unname(relations_id_to_name[relation_relations$child])

		suppressWarnings(suppressMessages(relations_DAG <- create_ontology_DAG(relation_relations$parent, relation_relations$child)))

		if(inherit_relations) {
			relation_type = merge_offspring_relation_types(relations_DAG, relation_type)
		}
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
		stop("Cannot find any owl:Class. Consider to use `import_ontology()` directly.")
	}

	if(verbose) message("Parsing ", length(Class), " <owl:Class> ...")
	id = xml_attr(Class, "about")
	short_id = .owl_get_text(Class, ".//*[local-name()='id']", NA_character_)
	short_id = ifelse(is.na(short_id), gsub("^.*#", "", basename(id)), short_id)
	name = .owl_get_text(Class, ".//rdfs:label[@xml:lang='en'] | .//rdfs:label[not(@xml:lang)]", NA_character_)
	def = .owl_get_text(Class, ".//*[local-name()='IAO_0000115']", NA_character_)
	namespace = .owl_get_text(Class, ".//*[local-name()='hasOBONamespace']", NA_character_)
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
		if(verbose) message("Parsing ", length(Description), " <rdf:Description> ...")
		id = xml_attr(Description, "about")
		name = .owl_get_text(Description, ".//*[local-name()='prefLabel']", NA_character_)
		def = .owl_get_text(Description, ".//*[local-name()='definition']", NA_character_)

		df_description = data.frame(id = id, name = name, def = def)
		df_description = df_description[!is.na(df_description$id), , drop = FALSE]
		df_description = df_description[!duplicated(df_description$id), , drop = FALSE]
		rownames(df_description) = df_description$id
	} else {
		df_description = data.frame()
	}


	## terms
	lt = .wrap_relations(lt_terms, "term")
	term_meta = lt$meta
	term_relations = lt$relations

	if(nrow(df_description) > 0) {
		ind = which(is.na(term_meta$name) & term_meta$id %in% df_description$id)
		if(length(ind)) {
			term_meta[ind, "name"] = df_description[ term_meta$id[ind], "name"]
			term_meta[ind, "definition"] = df_description[ term_meta$id[ind], "def"]
		}

		ind = which(is.na(term_relations$name) & term_meta$id %in% df_description$id)
		if(length(ind)) {
			term_relations[ind, "name"] = df_description[ term_relations$id[ind], "name"]
			term_relations[ind, "definition"] = df_description[ term_relations$id[ind], "def"]
		}
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

	term_meta = term_meta[term_meta$id %in% c(term_relations$parent, term_relations$child), , drop = FALSE]
	internal_id = term_meta$short_id
	dd = internal_id[duplicated(internal_id)]
	ldd = internal_id %in% dd
	internal_id[ldd] = term_meta$id[ldd]
	idmap = structure(internal_id, names = term_meta$id)

	dag = create_ontology_DAG(parents = idmap[term_relations$parent], children = idmap[term_relations$child], relations = term_relations$relation,
		source = paste0(ontology, ", ", version), relations_DAG = relations_DAG, verbose = verbose, ...)
	rownames(term_meta) = internal_id

	term_meta = term_meta[dag@terms, , drop = FALSE]
	
	if(dag_root(dag) == SUPER_ROOT) {
		term_meta$id[dag@root] = SUPER_ROOT
		term_meta$short_id[dag@root] = SUPER_ROOT
	}
	
	mcols(dag) = term_meta

	dag
	
}

#' @param robot_jar The path of the `robot.jar` file. It can be downloaded from https://github.com/ontodev/robot/releases.
#'         Internally, the file is converted to the obo format and parsed by `import_obo()`. The value of `robot_jar` can be
#'         set as a global option `simona_opt$robot_jar = ...`.
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
#' The description of the ROBOT tool is at \url{http://robot.obolibrary.org/convert}.
#' 
#' @rdname import_obo
#' @export
#' @importFrom utils download.file
#' @examples
#' \dontrun{
#' # The plant ontology: http://obofoundry.org/ontology/po.html 
#' dag = import_ontology("http://purl.obolibrary.org/obo/po.owl", robot_jar = ...)
#' }
import_ontology = function(file, robot_jar = simona_opt$robot_jar, JAVA_ARGS = "", verbose = simona_opt$verbose, ...) {

	if(grepl("\\.(obo|obo.gz)$", file, ignore.case = TRUE)) {
		return(import_obo(file, ...))
	}

	if(Sys.which("java") == "") {
		stop("Java is not available.")
	}

	if(is.null(robot_jar)) {
		stop("'robot.jar' has not been set. It can be downloaded from https://github.com/ontodev/robot/releases.")
	}

	if(!file.exists(robot_jar)) {
		stop("Cannot find 'robot.jar'. It can be downloaded from https://github.com/ontodev/robot/releases.")
	}
	robot_jar = normalizePath(robot_jar)

	if(grepl("^(http|ftp)", file)) {
		if(verbose) message(qq("Downloading @{file}..."))
		file2 = tempfile(fileext = paste0("_", basename(file)))
		download.file(file, destfile = file2, quiet = TRUE)
		on.exit(file.remove(file2))

		file = file2
	}

	if(verbose) message("Converting ", basename(file), " to the obo format.")
	output = tempfile(fileext = ".obo.gz")

	file = normalizePath(file)
	java_path = Sys.which("java")
	cmd = qq("'@{java_path}' @{JAVA_ARGS} -jar '@{robot_jar}' convert --input '@{file}' --format obo --output '@{output}' --check false")
	if(verbose) message("  ", cmd)

	code = system2(java_path, c(JAVA_ARGS, "-jar", robot_jar, "convert", "--input", file, "--format", "obo", "--output", output, "--check", "false"))
	if(code != 0) {
		if(grepl("\\.owl$", file, ignore.case = TRUE)) {
			message("Consider to use `import_owl()`")
		} else if(grepl("\\.ttl$", file, ignore.case = TRUE)) {
			message("Consider to use `import_ttl()`")
		}
		stop("Executing 'robot.jar' failed.")
	}

	lt = import_obo(output, verbose = verbose, ...)

	if(file.exists(output)) {
		file.remove(output)
	}

	lt
}


#' @details
#' `import_ttl()` is a simple parser for the `.ttl` format files. It only recognizes
#' terms that have the `owl:Class` object. The "is_a" relation is recognized by the predicate `rdfs:subClassOf`
#' or an ontology-specific predicate that contains `.*/isa`. Other relation types are defined with
#' the predicate `owl:ObjectProperty`. The format is parsed by a Perl script `system.file("scripts", "parse_ttl.pl", package = "simona")`.
#' @rdname import_obo
#' @export
#' @examples
#' \donttest{
#' # file is from https://bioportal.bioontology.org/ontologies/MSTDE
#' import_ttl("https://jokergoo.github.io/simona/MSTDE.ttl")
#' }
import_ttl = function(file, relation_type = "part_of", verbose = simona_opt$verbose, ...) {

	if(Sys.which("perl") == "") {
		stop("Perl is not available.")
	}

	if(grepl("^(http|ftp)", file)) {

		message(qq("Downloading @{file}..."))
		file2 = tempfile(fileext = paste0("_", basename(file)))
		download.file(file, destfile = file2, quiet = TRUE)
		on.exit(file.remove(file2))

		file = file2
		source = basename(file2)
	} else {
		source = basename(file)
	}

	file = normalizePath(file)

	if(verbose) message("Parsing .ttl file...")
	perl_script = system.file("scripts", "parse_ttl.pl", package = "simona")
	cmd = qq("perl '@{perl_script}' '@{file}'")
	if(length(relation_type)) {
		cmd = paste0(cmd, " ", paste("'", relation_type, "'", sep = "", collapse = " "))
	}
	df = read.csv(pipe(cmd))
		
	if(verbose) message("Constructing the DAG_ontology object...")
	lt_parents = strsplit(df$parent, ",")
	children = rep(df$id, times = vapply(lt_parents, length, FUN.VALUE = integer(1)))
	parents = unlist(lt_parents)

	lt_relations = strsplit(df$relation_type, ",")
	relations = unlist(lt_relations)

	dag = create_ontology_DAG(parents = parents, children = children, relations = relations,
		source = source, verbose = verbose, ...)

	term_meta = df[, seq_len(4), drop = FALSE]
	colnames(term_meta) = c("id", "name", "short_id", "definition")

	rownames(term_meta) = term_meta$id
	term_meta = term_meta[dag@terms, , drop = FALSE]

	if(dag_root(dag) == SUPER_ROOT) {
		term_meta$id[dag@root] = SUPER_ROOT
		term_meta$short_id[dag@root] = SUPER_ROOT
	}
	
	mcols(dag) = term_meta

	dag
}