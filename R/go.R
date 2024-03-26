
#' Create the ontology_DAG object from the GO.db package
#' 
#' @param namespace One of "BP", "CC" and "MF".
#' @param relations Types of the GO term relations. In the **GO.db** package, the GO term relations can be "is_a", "part_of",
#'               "regulates", "negatively regulates", "positively regulates". Note since "regulates" is a parent relation
#'               of "negatively regulates", "positively regulates", if "regulates" is selected, "negatively regulates" and "positively regulates"
#'               are also selected. Note "is_a" is always included.
#' @param org_db The name of the organism package or the corresponding database object, e.g. `"org.Hs.eg.db"` or 
#'            directly the [`org.Hs.eg.db::org.Hs.eg.db`] object for human, then the gene annotation to GO terms will be added
#'            to the object. For other non-model organisms, consider to use the **AnnotationHub** package to find one.
#' @param evidence_code A vector of evidence codes for gene annotation to GO terms. See \url{https://geneontology.org/docs/guide-go-evidence-codes/}.
#' @param retrieve_alternative Whether to retrieve alternative/obsolete GO terms from geneontology.org?
#' @param verbose Whether to print messages.
#' 
#' @return An `ontology_DAG` object.
#' @export
#' @importFrom utils getFromNamespace
#' @examples
#' dag = create_ontology_DAG_from_GO_db()
#' dag
create_ontology_DAG_from_GO_db = function(namespace = "BP", relations = "part of", org_db = NULL, 
	evidence_code = NULL, retrieve_alternative = FALSE, verbose = simona_opt$verbose) {

	check_pkg("GO.db", bioc = TRUE)

	if(namespace == "BP") {
		df = AnnotationDbi::toTable(GO.db::GOBPCHILDREN)
	} else if(namespace == "CC") {
		df = AnnotationDbi::toTable(GO.db::GOCCCHILDREN)
	} else if(namespace == "MF") {
		df = AnnotationDbi::toTable(GO.db::GOMFCHILDREN)
	} else {
		stop("Value of `namespace` can only be one of 'BP', 'MF' and 'CC'.")
	}

	l = df[, 3] == "isa"
	df[l, 3] = "is_a"
	df[, 3] = gsub(" ", "_", df[, 3])

	if(length(relations) == 0) {
		relations = character(0)
	}

	if(identical(relations, NA)) {
		relations = character(0)
	}

	relations = c("is_a", relations)
	relations = normalize_relation_type(relations)
	if("regulates" %in% relations) {
		relations = c(relations, "negatively_regulates", "positively_regulates")
	}
	df = df[df[, 3] %in% relations, , drop = FALSE]

	if(verbose) {
		message("relations: ", paste(relations, collapse = ", "))
	}

	if(!is.null(org_db)) {
		if(is.character(org_db)) {
			check_pkg(org_db, bioc = TRUE)
			org_db = getFromNamespace(org_db, ns = org_db)
		}
	
		suppressMessages(tb <- AnnotationDbi::select(org_db, keys = AnnotationDbi::keys(org_db), columns = c("GO", "EVIDENCE", "ONTOLOGY")))
		tb = tb[tb$ONTOLOGY == namespace, , drop = FALSE]
		if(!is.null(evidence_code)) {
			tb = tb[tb$EVIDENCE %in% evidence_code, , drop = FALSE]
			if(row(tb) == 0) {
				stop("No annotation left after filtering by the evidence codes.")
			}
		}
		tb = tb[, c(1, which(colnames(tb) == "GO")), drop = FALSE]
		annotation = split(tb[, 1], tb[, 2])
		annotation = lapply(annotation, unique)
		annotation = lapply(annotation, as.character)
	} else {
		annotation = NULL
	}

	relations_DAG = create_ontology_DAG(c("regulates", "regulates"), c("negatively regulates", "positively regulates"))

	alternative_terms = list()
	if(retrieve_alternative) {
		alternative_terms = alternative_GO_terms(verbose = verbose)
	}

	go_db_version = read.dcf(system.file("DESCRIPTION", package = "GO.db"))[1, "Version"]
	dag = create_ontology_DAG(parents = df[, 2], children = df[, 1], relations = df[, 3], relations_DAG = relations_DAG,
		annotation = annotation, source = paste0("GO ", namespace, " / GO.db package ", go_db_version), alternative_terms = alternative_terms)

	go = GO.db::GOTERM[dag@terms]
	meta = data.frame(id = AnnotationDbi::GOID(go),
		        name = AnnotationDbi::Term(go),
		        definition = AnnotationDbi::Definition(go))
	rownames(meta) = dag@terms

	mcols(dag) = meta

	dag
}


#' Mappings between alternative GO terms to official GO terms
#' 
#' @param tag In the `go-basic.obo` file, there are three tags which define alternative GO terms: `replaced_by`, `alt_id` and `consider`.
#'        See https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html#S.2.2.1
#' @param version Version of the `go-basic.obo` file. By default it is the version for building **GO.db** package. The value is a string in the format of "2024-01-17".
#' 
#' @return A list of named vectors where names are alternative GO IDs and value vectors are current GO IDs in use.
#' @export
#' @rdname create_ontology_DAG_from_GO_db
alternative_GO_terms = function(tag = c("replaced_by", "alt_id", "consider"), version = NULL, verbose = TRUE) {

	db_info = GO.db::GO_dbInfo()
	if(is.null(version)) {
		obo_version = db_info$value[db_info$name == "GOSOURCEDATE"]
	} else {
		obo_version = version
	}

	obo_url = paste0("https://release.geneontology.org/", obo_version, "/ontology/go-basic.obo")

	if(verbose) {
		message("Use GO obo version: ", obo_version)
		message("Downloading ", obo_url, " to get the alternative terms.")
	}
	op = getOption("timeout")
	options(timeout = 999999999)
	on.exit(options(timeout = op))
	con = url(obo_url)
	ln = readLines(con)
	close(con)

	if(verbose) {
		message("Retrieving mappings between alternative terms to DAG terms.")
	}
	ind_stanza = c(grep("^\\[.*\\]$", ln), length(ln))
	ind1 = grep("^\\[Term\\]$", ln)

	ind2 = cpp_match_index(ind1, ind_stanza[-1]-1)
	n = length(ind1)

	tag = tag
	process_obo_alternative = function(ln) {

		lt = list()
		i = grep("^id:", ln)
		lt$id = gsub("^id: (\\S+)(\\s*.*)$", "\\1", ln[i])[1]

		lt$alt_id = character(0)
		lt$is_obsolete = "false"
		lt$replaced_by = character(0)
		lt$consider = character(0)
		
		if("alt_id" %in% tag) {
			i = grep("^alt_id:", ln)
			if(length(i)) {
				lt$alt_id = gsub("^alt_id: (.*)$", "\\1", ln[i])
			}
		}
		
		i = grep("^is_obsolete:", ln)
		if(length(i)) {
			lt$is_obsolete = gsub("^is_obsolete: (.*)$", "\\1", ln[i])

			if("replaced_by" %in% tag) {
				i = grep("^replaced_by:", ln)
				if(length(i)) {
					lt$replaced_by = gsub("^replaced_by: (.*)$", "\\1", ln[i])
				}
			}

			if("consider" %in% tag) {
				i = grep("^consider:", ln)
				if(length(i)) {
					lt$consider = gsub("^consider: (.*)$", "\\1", ln[i])
				}
			}
		}

		lt
	}

	lt_terms = vector("list", n)
	for(i in seq_len(n)) {
		lt_terms[[i]] = process_obo_alternative(ln[seq(ind1[i], ind2[i])])
	}

	alternative_terms1 = list()
	alternative_terms2 = list()
	alternative_terms3 = list()

	if("replaced_by" %in% tag) {
		l = sapply(lt_terms, function(x) length(x$replaced_by) > 0)
		alternative_terms1 = lt_terms[l]
		lt = lapply(alternative_terms1, function(x) {
			data.frame(id = rep(x$id, length(x$replaced_by)), replaced_by = x$replaced_by)
		})
		alternative_terms1 = do.call(rbind, lt)
		alternative_terms1 = split(alternative_terms1$replaced_by, alternative_terms1$id)
		alternative_terms1 = lapply(alternative_terms1, function(x) structure(x, names = rep("replaced_by", length(x))))
	}

	if("consider" %in% tag) {
		l = sapply(lt_terms, function(x) length(x$consider) > 0)
		alternative_terms2 = lt_terms[l]
		lt = lapply(alternative_terms2, function(x) {
			data.frame(id = rep(x$id, length(x$consider)), consider = x$consider)
		})
		alternative_terms2 = do.call(rbind, lt)
		alternative_terms2 = split(alternative_terms2$consider, alternative_terms2$id)
		alternative_terms2 = lapply(alternative_terms2, function(x) structure(x, names = rep("consider", length(x))))
	}

	if("alt_id" %in% tag) {
		l = sapply(lt_terms, function(x) length(x$alt_id) > 0)
		alternative_terms3 = lt_terms[l]
		alt_lt = lapply(alternative_terms3, function(x) {
			data.frame(id = rep(x$id, length(x$alt_id)), alt_id = x$alt_id)
		})
		alternative_terms3 = do.call(rbind, alt_lt)
		alternative_terms3 = split(alternative_terms3$id, alternative_terms3$alt_id)
		alternative_terms3 = lapply(alternative_terms3, function(x) structure(x, names = rep("alt_id", length(x))))
	}

	cn = unique(c(names(alternative_terms1), names(alternative_terms2), names(alternative_terms3)))
	lt = vector("list", length(cn))
	names(lt) = cn
	for(x in names(alternative_terms1)) {
		lt[[x]] = c(lt[[x]], alternative_terms1[[x]])
	}
	for(x in names(alternative_terms2)) {
		lt[[x]] = c(lt[[x]], alternative_terms2[[x]])
	}
	for(x in names(alternative_terms3)) {
		lt[[x]] = c(lt[[x]], alternative_terms3[[x]])
	}

	lt
}

