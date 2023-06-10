
term_to_node_id = function(dag, term, strict = TRUE) {
	if(is.numeric(term)) {
		term
	} else if(length(term) == 1) {
		i = which(dag@terms == term)

		if(length(i) == 0) {
			stop("Cannot find term: ", term)
		}

		i
	} else if(length(term) > 1) {
		i = which(dag@terms %in% term)

		if(length(i) == 0) {
			stop("Cannot find all these terms.")
		}
		if(length(i) != length(term)) {
			if(strict) {
				stop("Cannot find some of the terms in the DAG.")
			} else {
				message("removed", length(term) - length(i), "terms that cannot be found in the DAG.")
			}
		}

		i
	}
}


check_pkg = function(pkg, bioc = FALSE, github = NULL) {
	if(requireNamespace(pkg, quietly = TRUE)) {
		return(NULL)
	} else {

		if(!interactive()) {
			if(bioc) {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}

		if(bioc) {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from Bioconductor? [y|n] "))
		} else {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from CRAN? [y|n] "))
		}

		if(bioc) {
			if(tolower(answer) %in% c("y", "yes")) {
				if(!requireNamespace("BiocManager", quietly = TRUE)) {
					install.packages("BiocManager")
				}
				suppressWarnings(BiocManager::install(pkg, update = FALSE))

				if(!requireNamespace(pkg, quietly = TRUE)) {
					if(is.null(github)) {
						stop_wrap(qq("Cannot find '@{pkg}' from Bioconductor."))
					} else {
						answer = readline(qq("Not on Bioconductor. Install '@{pkg}' from GitHub: '@{github}/@{pkg}'? [y|n] "))
						if(tolower(answer) %in% c("y", "yes")) {
							BiocManager::install(paste0(github, "/", pkg), update = FALSE)
						} else {
							stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
						}
					}
				}
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			}
		} else {
			if(tolower(answer) %in% c("y", "yes")) {
				suppressWarnings(install.packages(pkg))
				if(!requireNamespace(pkg, quietly = TRUE)) {
					if(is.null(github)) {
						stop_wrap(qq("Cannot find '@{pkg}' from CRAN"))
					} else {
						answer = readline(qq("Not on CRAN. Install '@{pkg}' from GitHub: '@{github}/@{pkg}'? [y|n] "))
						if(tolower(answer) %in% c("y", "yes")) {
							BiocManager::install(paste0(github, "/", pkg), update = FALSE)
						} else {
							stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
						}
					}
				}
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}
	}
}


stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

warning_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    warning(x, call. = FALSE)
}

message_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}


add_transparency = function (col, transparency = 0) {
    rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
}
