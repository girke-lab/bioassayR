setMethod("show", signature=signature(
    object="BioassayDB"),
    function(object) {
        cat("class:\t\t\t", class(object), "\n")
        cat("assays:\t\t\t", queryBioassayDB(object, "SELECT COUNT(*) FROM assays")[[1]], "\n")
    	cat("sources:\t\t", paste(queryBioassayDB(object, "SELECT description FROM sources")[[1]], sep=", "), "\n")
    	cat("source versions:\t", paste(queryBioassayDB(object, "SELECT version FROM sources")[[1]], sep=", "), "\n")
        cat("writeable:\t\t")
        if(.writeable(object)){
            cat("yes\n")
        } else {
            cat("no\n")
        }
    }
)

setMethod("queryBioassayDB", signature(object="BioassayDB"),
    function(object, query) {
        # input tests:
        if(! is.character(query))
            stop("sql query not class character")
        if(length(query) > 1){
            warning("too many sql queries, only the first was processed")
            query <- query[1]
        }

        dbGetQuery(slot(object, "database"), query)
    }
)

###############################
# internal use only functions #
###############################

.writeable <- function(object){
    if(slot(slot(object, "database"), "flags") == SQLITE_RO){
        return(FALSE)
    } else {
        return(TRUE)
    }
}
