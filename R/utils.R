
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
		unique_term = unique(term)
		i = which(dag@terms %in% unique_term)

		if(length(i) == 0) {
			stop("Cannot find all these terms.")
		}

		if(length(i) != length(unique_term)) {
			if(strict) {
				stop("Cannot find some of the terms in the DAG.")
			} else {
				message("removed ", length(unique_term) - length(i), " terms that cannot be found in the DAG.")
			}
		}

		unname(structure(i, names = dag@terms[i])[term])
	}
}


#' @importFrom GetoptLong qq
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

#' @importFrom grDevices col2rgb rgb
add_transparency = function (col, transparency = 0, alpha = TRUE) {
	if(alpha) {
    	rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
    } else {

		m = col2rgb(col)
		m = 255 - (255-m)*(1-transparency)
		rgb(t(m), maxColorValue = 255)

    }
}


lt_children_to_lt_parents = function(lt_children) {
	n = length(lt_children)
	parents = rep(seq_len(n), times = sapply(lt_children, length))
	children = unlist(lt_children)

	lt = split(parents, children)

	lt_parents = rep(list(integer(0)))
	lt_parents[ as.integer(names(lt)) ] = lt

	lt_parents
}

lt_parents_to_lt_children = function(lt_parents) {
	n = length(lt_parents)
	children = rep(seq_len(n), times = sapply(lt_parents, length))
	parents = unlist(lt_parents)

	lt = split(children, parents)

	lt_children = rep(list(integer(0)))
	lt_children[ as.integer(names(lt)) ] = lt

	lt_children
}

dag_is_tree = function(dag) {
	n_terms = dag@n_terms
	n_relations = sum(sapply(dag@lt_children, length))

	n_terms == n_relations + 1
}


merge_offspring_relation_types = function(relations_DAG, relations) {
	r1 = relations_DAG@terms
	rc = intersect(r1, relations)

	if(length(rc)) {
		unique(c(setdiff(relations, rc), dag_offspring(relations_DAG, rc, include_self = TRUE)))
	} else {
		relations
	}
}


# all offspring types are assigned to the same value
extend_contribution_factor = function(relations_DAG, contribution_factor) {

	cf = contribution_factor
	
	if(is.null(relations_DAG)) {
		return(cf)
	}

	for(nm in names(contribution_factor)) {
		if(nm %in% relations_DAG@terms) {
			offspring = dag_offspring(relations_DAG, nm)
			if(length(offspring)) {
				cf[offspring] = contribution_factor[nm]
			}
		}
	}

	cf
}


normalize_relation_type = function(x) {

	x = tolower(x)
	x = gsub("[- ~]", "_", x)
	x[x == "isa"] = "is_a"
	x[x == "part_a"] = "part_of"

	x
}
