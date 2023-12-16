
#' Import ontologies already having gene annotations
#' 
#' @param organism Organism.
#' @param verbose Whether to print messages?
#' @param ... Pass to [`create_ontology_DAG()`].
#' 
#' @details
#' There are the following ontologies:
#' 
#' - `ontology_kw()`: UniProt Keywords. The list of supported organisms can be found in [`UniProtKeywords::load_keyword_genesets()`].
#' - `ontology_chebi()`: Chemical Entities of Biological Interest.
#' - `ontology_hp()`: The Human Phenotype Ontology.
#' - `ontology_pw()`: Pathway Ontology.
#' - `ontology_rdo()`: RGD Disease Ontology.
#' - `ontology_vt()`: Vertebrate Trait Ontology.
#' 
#' The source of the original files can be found with `simona:::RGD_TB`.
#' 
#' @export
#' @importFrom utils data read.table
#' @rdname ontology
ontology_kw = function(organism = "human", verbose = simona_opt$verbose, ...) {

	check_pkg("UniProtKeywords", bioc = TRUE)

	if(verbose) {
		message("Obtaining ontology relations from the UniProtKeywords package...")
	}
	kw_parents = NULL
	load(system.file("data", "kw_parents.rda", package = "UniProtKeywords"))

	parents = unlist(kw_parents)
	children = rep(names(kw_parents), times = sapply(kw_parents, length))

	kw_terms = NULL
	load(system.file("data", "kw_terms.rda", package = "UniProtKeywords"))

	id = sapply(kw_terms, function(x) x$Identifier)
	accession = sapply(kw_terms, function(x) x$Accession)
	map = structure(accession, names = id)
	parents = map[parents]
	children = map[children]

	dag = create_ontology_DAG(parents, children, source = "UniProt Keywords", verbose = verbose, ...)

	meta = data.frame(
		id = sapply(kw_terms, function(x) x$Identifier),
		accession = sapply(kw_terms, function(x) x$Accession),
		name = sapply(kw_terms, function(x) x$Identifier),
		description = sapply(kw_terms, function(x) x$Description),
		category = sapply(kw_terms, function(x) paste(x$Category, collapse = "; "))
	)
	rownames(meta) = meta$accession
	meta = meta[dag@terms, ]
	rownames(meta)[nrow(meta)] = SUPER_ROOT

	mcols(dag) = meta

	annotation = UniProtKeywords::load_keyword_genesets(organism)
	names(annotation) = map[names(annotation)]
	add_annotation(dag, annotation)

}


RGD_TB = read.table(textConnection(
"org	onto	anno_url	onto_url
canis	chebi	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/canis_genes_chebi	https://purl.obolibrary.org/obo/chebi.obo
homo	chebi	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/homo_genes_chebi	https://purl.obolibrary.org/obo/chebi.obo
mus	chebi	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/mus_genes_chebi	https://purl.obolibrary.org/obo/chebi.obo
rattus	chebi	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/rattus_genes_chebi	https://purl.obolibrary.org/obo/chebi.obo
sus	chebi	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/sus_genes_chebi	https://purl.obolibrary.org/obo/chebi.obo
homo	hp	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/homo_genes_hp	https://purl.obolibrary.org/obo/hp.obo
mus	hp	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/mus_genes_hp	https://purl.obolibrary.org/obo/hp.obo
canis	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/canis_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
chinchilla	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chinchilla_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
chlorocebus	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chlorocebus_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
heterocephalus	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/heterocephalus_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
homo	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/homo_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
ictidomys	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/ictidomys_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
mus	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/mus_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
pan	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/pan_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
rattus	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/rattus_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
sus	pw	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/sus_genes_pw	https://download.rgd.mcw.edu/ontology/pathway/pathway.obo
canis	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/canis_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
chinchilla	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chinchilla_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
chlorocebus	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chlorocebus_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
heterocephalus	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/heterocephalus_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
homo	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/homo_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
ictidomys	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/ictidomys_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
mus	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/mus_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
pan	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/pan_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
rattus	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/rattus_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
sus	rdo	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/sus_genes_rdo	https://download.rgd.mcw.edu/ontology/disease/RDO.obo
canis	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/canis_genes_vt	https://purl.obolibrary.org/obo/vt.owl
chinchilla	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chinchilla_genes_vt	https://purl.obolibrary.org/obo/vt.owl
chlorocebus	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/chlorocebus_genes_vt	https://purl.obolibrary.org/obo/vt.owl
heterocephalus	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/heterocephalus_genes_vt	https://purl.obolibrary.org/obo/vt.owl
homo	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/homo_genes_vt	https://purl.obolibrary.org/obo/vt.owl
ictidomys	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/ictidomys_genes_vt	https://purl.obolibrary.org/obo/vt.owl
mus	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/mus_genes_vt	https://purl.obolibrary.org/obo/vt.owl
pan	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/pan_genes_vt	https://purl.obolibrary.org/obo/vt.owl
rattus	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/rattus_genes_vt	https://purl.obolibrary.org/obo/vt.owl
sus	vt	https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/sus_genes_vt	https://purl.obolibrary.org/obo/vt.owl
"), header = TRUE)


organism_map = c("human" = "homo",
	             "mouse" = "mus",
	             "rat" = "rattus",
	             "pig" = "sus", 
	             "dog" = "canis", 
	             "chimpanzee" = "pan")

#' @rdname ontology
#' @export
ontology_chebi = function(
	organism = c("human", "mouse", "rat", "pig", "dog"), 
	verbose = simona_opt$verbose, ...) {
	
	organism = match.arg(organism)[1]
	organism = organism_map[organism]

	i = which(RGD_TB$org == organism & RGD_TB$onto == "chebi")
	onto_url = RGD_TB$onto_url[i]
	anno_url = RGD_TB$anno_url[i]

	if(verbose) {
		message(qq("Downloading ontology file from @{onto_url}..."))
	}
	dag = import_ontology(onto_url, verbose = verbose, ...)

	if(verbose) {
		message(qq("Downloading annotation file from @{anno_url}..."))
	}
	anno_tb = read.table(url(anno_url), comment.char = "!", sep = "\t", quote = "")
	anno = split(anno_tb[[3]], anno_tb[[5]])

	add_annotation(dag, anno)
}


#' @rdname ontology
#' @export
ontology_hp = function(
	organism = c("human", "mouse"), verbose = simona_opt$verbose, ...) {
	
	organism = match.arg(organism)[1]
	organism = organism_map[organism]

	i = which(RGD_TB$org == organism & RGD_TB$onto == "hp")
	onto_url = RGD_TB$onto_url[i]
	anno_url = RGD_TB$anno_url[i]

	if(verbose) {
		message(qq("Downloading ontology file from @{onto_url}..."))
	}
	dag = import_ontology(onto_url, verbose = verbose, ...)

	if(verbose) {
		message(qq("Downloading annotation file from @{anno_url}..."))
	}
	anno_tb = read.table(url(anno_url), comment.char = "!", sep = "\t", quote = "")
	anno = split(anno_tb[[3]], anno_tb[[5]])

	add_annotation(dag, anno)
}

#' @rdname ontology
#' @export
ontology_pw = function(
	organism = c("human", "mouse", "rat", "pig", "dog", "chimpanzee"),
	verbose = simona_opt$verbose, ...) {
	
	organism = match.arg(organism)[1]
	organism = organism_map[organism]

	i = which(RGD_TB$org == organism & RGD_TB$onto == "pw")
	onto_url = RGD_TB$onto_url[i]
	anno_url = RGD_TB$anno_url[i]

	if(verbose) {
		message(qq("Downloading ontology file from @{onto_url}..."))
	}
	dag = import_ontology(onto_url, verbose = verbose, ...)

	if(verbose) {
		message(qq("Downloading annotation file from @{anno_url}..."))
	}
	anno_tb = read.table(url(anno_url), comment.char = "!", sep = "\t", quote = "")
	anno = split(anno_tb[[3]], anno_tb[[5]])

	add_annotation(dag, anno)
}

#' @rdname ontology
#' @export
ontology_rdo = function(
	organism = c("human", "mouse", "rat", "pig", "dog", "chimpanzee"),
	verbose = simona_opt$verbose, ...) {
	
	organism = match.arg(organism)[1]
	organism = organism_map[organism]

	i = which(RGD_TB$org == organism & RGD_TB$onto == "rdo")
	onto_url = RGD_TB$onto_url[i]
	anno_url = RGD_TB$anno_url[i]

	if(verbose) {
		message(qq("Downloading ontology file from @{onto_url}..."))
	}
	dag = import_ontology(onto_url, verbose = verbose, ...)

	if(verbose) {
		message(qq("Downloading annotation file from @{anno_url}..."))
	}
	anno_tb = read.table(url(anno_url), comment.char = "!", sep = "\t", quote = "")
	anno = split(anno_tb[[3]], anno_tb[[5]])

	add_annotation(dag, anno)
}

#' @rdname ontology
#' @export
ontology_vt = function(
	organism = c("human", "mouse", "rat", "pig", "dog", "chimpanzee"),
	verbose = simona_opt$verbose, ...) {
	
	organism = match.arg(organism)[1]
	organism = organism_map[organism]

	i = which(RGD_TB$org == organism & RGD_TB$onto == "vt")
	onto_url = RGD_TB$onto_url[i]
	anno_url = RGD_TB$anno_url[i]

	if(verbose) {
		message(qq("Downloading ontology file from @{onto_url}..."))
	}
	dag = import_ontology(onto_url, verbose = verbose, ...)

	if(verbose) {
		message(qq("Downloading annotation file from @{anno_url}..."))
	}
	anno_tb = read.table(url(anno_url), comment.char = "!", sep = "\t", quote = "")
	anno = split(anno_tb[[3]], anno_tb[[5]])

	add_annotation(dag, anno)
}

#' @details
#' `ontology_go()` is an alias of [`create_ontology_DAG_from_GO_db()`]. All arguments go there.
#' @rdname ontology
#' @export
ontology_go = function(...) {
	create_ontology_DAG_from_GO_db(...)
}
