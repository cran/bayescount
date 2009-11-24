.onLoad <- function(libname, pkgname) library.dynam("bayescount", pkgname, libname)
.onUnload <- function(libpath) library.dynam.unload("bayescount", libpath)

# from ?.onLoad
# Anything needed for the functioning of the name space should be handled at load/unload times by the .onLoad and .onUnload hooks. For example, DLLs can be loaded (unless done by a useDynLib directive in the ‘NAMESPACE’ file) and initialized in .onLoad and unloaded in .onUnload. Use .onAttach only for actions that are needed only when the package becomes visible to the user, for example a start-up message.

#.First.lib <- function(lib,pkg)
#{
#   library.dynam("bayescount",pkg,lib)
#   cat("power loaded\n")
#}

#.Last.lib <- function(){
#	library.dynam.unload("bayescount")	
#}

##  If I weren't using a NAMESPACE then I'd use .First.lib and .last.lib instead of .onLoad and .onUnload