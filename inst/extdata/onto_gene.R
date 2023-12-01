## uniprot keywords


library(UniProtKeywords)

data(kw_parents)

parents = unlist(kw_parents)
children = rep(names(kw_parents), times = sapply(kw_parents, length))

dag = create_ontology_DAG(parents, children)

data(kw_terms)

meta = data.frame(
	id = sapply(kw_terms, function(x) x$Identifier),
	accession = sapply(kw_terms, function(x) x$Accession),
	name = sapply(kw_terms, function(x) x$Identifier),
	description = sapply(kw_terms, function(x) x$Description),
	category = sapply(kw_terms, function(x) paste(x$Category, collapse = "; "))
)
rownames(meta) = meta$id
meta = meta[dag@terms, ]
rownames(meta)[nrow(meta)] = simona:::SUPER_ROOT

mcols(dag) = meta

annotation = load_keyword_genesets("9606")
dag = add_annotation(dag, annotation)


## The Human Phenotype Ontology

# https://hpo.jax.org/app/data/annotations

dag = import_obo("https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-10-09/hp-base.obo")

tb = read.table(url("https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-10-09/genes_to_phenotype.txt"), sep = "\t", header = TRUE)
annotation = split(tb$ncbi_gene_id, tb$hpo_id)

dag = add_annotation(dag, annotation)


## Pathway ontology and many

# https://download.rgd.mcw.edu/ontology/

library(rvest)

html = read_html("https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/")
tb = html %>% html_element("table") %>% html_table()
fn = tb[grep("_genes_", tb[[2]]), ][[2]]
tb = data.frame(org = gsub("_.*$", "", fn),
	            onto = gsub("^.*_", "", fn))
tb = tb[!tb$onto %in% c("go", "nbo", "mp", "cmo"), ]
tb = tb[order(tb$onto), ]

tb$anno_url = paste0("https://download.rgd.mcw.edu/ontology/annotated_rgd_objects_by_ontology/", tb$org, "_genes_", tb$onto)

onto = c("chebi" = "https://purl.obolibrary.org/obo/chebi.obo",
	     "pw" = "https://download.rgd.mcw.edu/ontology/pathway/pathway.obo",
	     "rdo" = "https://download.rgd.mcw.edu/ontology/disease/RDO.obo",
	     "vt" = "https://purl.obolibrary.org/obo/vt.owl",
	     "hp" = "https://purl.obolibrary.org/obo/hp.obo")

tb$onto_url = onto[tb$onto]
## MeSH

dag = import_ttl("https://data.bioontology.org/ontologies/MESH/submissions/26/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb")


