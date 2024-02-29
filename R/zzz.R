
#' @importFrom utils packageDescription
.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

    msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/simona/
Github page: https://github.com/jokergoo/simona
Documentation: https://jokergoo.github.io/simona/

If you use it in published research, please cite:
Gu, Z. simona: a Comprehensive R package for Semantic Similarity 
  Analysis on Bio-Ontologies. bioRxiv 2023.

This message can be suppressed by:
  suppressPackageStartupMessages(library(simona))
========================================
")  

    packageStartupMessage(msg)
}


finalize = function(env) {
   if(!is.null(simona_opt$robot_jar)) {
		if(grepl("robot_temp_", simona_opt$robot_jar)) {
			if(file.exists(simona_opt$robot_jar)) {
				file.remove(simona_opt$robot_jar)
				simona_opt$robot_jar = NULL
			}
		}
	}
}

.onLoad = function(libname, pkgname) {
   parent = parent.env(environment())
   reg.finalizer(parent, finalize, onexit = TRUE)
}

.onUnload = function(libpath) {
	if(!is.null(simona_opt$robot_jar)) {
		if(grepl("robot_temp_", simona_opt$robot_jar)) {
			if(file.exists(simona_opt$robot_jar)) {
				file.remove(simona_opt$robot_jar)
				simona_opt$robot_jar = NULL
			}
		}
	}
}
