

# first fork https://github.com/OBOFoundry/purl.obolibrary.org
library(yaml)
options(timeout = 9999999)
for(f in list.files(path = "~/project/development/purl.obolibrary.org/config", full.names = TRUE)) {
	ln = read_yaml(f)
	url = unlist(ln$products)
	for(x in url) {
		dest = paste0("~/project/development/purl.obolibrary.org/files/", basename(x))
		if(file.exists(dest)) {
			next
		}
		oe = try(download.file(x, dest = dest))
		if(inherits(oe, "try-error")) {
			file.remove(dest)
		}
	}

}


options(timeout = 9999999)

## https://bioportal.bioontology.org/
apikey = "5cde006c-26ea-44a2-991b-e60b997d27dc"
js = fromJSON(qq("https://data.bioontology.org/ontologies?apikey=@{apikey}"))

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
	ln = header[grepl("Content-Disposition: attachment; filename", header)]

	if(length(ln)) {
		dest = gsub('^Content-Disposition: attachment; filename="(.*)".*$', "\\1", ln)
	} else {
		next
	}

	if(file.exists(qq("~/temp/@{dest}"))) {
		next
	}

	qqcat("  hasOntologyLanguage: @{hasOntologyLanguage}\n")
	qqcat("  download: @{submission_id}/download\n")
	qqcat("  local: @{dest}\n")
	cat("\n")

	download.file(qq("@{submission_id}/download?apikey=@{apikey}"), dest = qq("~/temp/@{dest}"))
}
