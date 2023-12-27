

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
