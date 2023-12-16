
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
