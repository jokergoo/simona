
setwd("~/workspace/ontology")


# http://obofoundry.org/


dir.create("OBOFoundry", showWarnings = FALSE)
setwd("OBOFoundry")

library(jsonlite)

lt = fromJSON("http://obofoundry.org/registry/ontologies.jsonld")
tb = lt$ontologies
saveRDS(tb, file = "OBOFoundry_meta_table.rds")

library(rvest)
options(timeout = 9999999)

for(i in seq_len(nrow(tb))) {
	qqcat("============ @{tb[i, 'id']} [@{i}/@{nrow(tb)}] ============\n")

	html = read_html(qq("http://obofoundry.org/ontology/@{tb[i, 'id']}.html"))
	nodes = html %>% html_elements(xpath = "//div[contains (text(), 'Products')]/following-sibling::table/*/tr/td[1]/a[@href]")
	file = nodes %>% html_text()
	url = nodes %>% html_attr("href")

	dir.create(tb$id[i], showWarnings = FALSE)

	for(j in seq_along(file)) {
		dest = qq("@{tb$id[i]}/@{basename(file[j])}")

		oe = try(header <- curlGetHeaders(url[j]))
		if(!inherits(oe, "try-error")) {
			ln2 = header[grepl("Content-Length", header, ignore.case = TRUE)]

			filesize = max(as.numeric(gsub('^Content-Length: (\\d+)\\s*$', "\\1", ln2, ignore.case = TRUE)))
			
			if(file.exists(dest)) {
				if(file.info(dest)[1, "size"] == filesize) {
					qqcat("already downloaded, skip.\n")
					next
				}
			}
		}
		oe = try(download.file(url[j], dest = dest))
		if(inherits(oe, "try-error")) {
			file.remove(dest)
		}
	}
}



##########################

setwd("~/workspace/ontology")

dir.create("BioPortal", showWarnings = FALSE)
setwd("BioPortal")

options(timeout = 9999999)

## https://bioportal.bioontology.org/
apikey = 
js = fromJSON(qq("https://data.bioontology.org/ontologies?apikey=@{apikey}"))

saveRDS(js, file = "BioPortal_meta_table.rds")


for(i in seq_len(nrow(js))) {
	acronym = js$acronym[i]
	submissions = js$links$submissions[i]

	qqcat("============ @{acronym} [@{i}/@{nrow(js)}] ============\n")

	oe = try(sub <- fromJSON(qq("@{submissions}?apikey=@{apikey}")))

	if(inherits(oe, "try-error")) {
		next
	}

	if(length(sub) == 0) {
		next
	}

	hasOntologyLanguage = sub$hasOntologyLanguage[1]
	submission_id = sub[1, "@id"]

	url = qq("@{submission_id}/download?apikey=@{apikey}")
	header = curlGetHeaders(url)
	ln = header[grepl("Content-Disposition: attachment; filename", header, ignore.case = TRUE)]
	ln2 = header[grepl("Content-Length", header, ignore.case = TRUE)]

	if(length(ln)) {
		dest = gsub('^Content-Disposition: attachment; filename="(.*)".*$', "\\1", ln, ignore.case = TRUE)
		filesize = as.numeric(gsub('^Content-Length: (\\d+)\\s*$', "\\1", ln2, ignore.case = TRUE))
	} else {
		next
	}

	dest = qq("@{acronym}/@{dest}")

	dir.create(acronym, showWarnings = FALSE)

	if(file.exists(dest)) {
		qqcat("already downloaded, skip.\n")
		if(file.info(dest)[1, "size"] == filesize) {
			next
		}
	}

	qqcat("  hasOntologyLanguage: @{hasOntologyLanguage}\n")
	qqcat("  download: @{submission_id}/download\n")
	qqcat("  local: @{dest}\n")
	cat("\n")

	download.file(qq("@{submission_id}/download?apikey=@{apikey}"), dest = dest)
}

