setMethod("show", signature=signature(
    object="bioassaySet"),
    function(object) {
        cat("class:\t\t", class(object), "\n")
        cat("assays:\t\t", dim(slot(object, "activity"))[1], "\n")
        cat("compounds:\t", dim(slot(object, "activity"))[2], "\n")
        cat("targets:\t", length(colnames(slot(object, "targets"))), "\n")
        cat("sources:\t", paste(unique(slot(object, "sources")$description), sep=", "), "\n")
    }
)

# activity
setMethod(f="activity", signature="bioassaySet", definition=function(x) {return(x@activity)})
setReplaceMethod(f="activity", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "dgCMatrix")
        stop("invalid input: must be of class dgCMatrix") 
    x@activity <- value
    return(x)
})

# scores
setMethod(f="scores", signature="bioassaySet", definition=function(x) {return(x@scores)})
setReplaceMethod(f="scores", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "dgCMatrix")
        stop("invalid input: must be of class dgCMatrix") 
    x@scores <- value
    return(x)
})

# targets
setMethod(f="targets", signature="bioassaySet", definition=function(x) {return(x@targets)})
setReplaceMethod(f="targets", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "dgCMatrix")
        stop("invalid input: must be of class dgCMatrix") 
    x@targets <- value
    return(x)
})

# sources
setMethod(f="sources", signature="bioassaySet", definition=function(x) {return(x@sources)})
setReplaceMethod(f="sources", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "data.frame")
        stop("invalid input: must be of class data.frame") 
    x@sources <- value
    return(x)
})

# source_id
setMethod(f="source_id", signature="bioassaySet", definition=function(x) {return(x@source_id)})
setReplaceMethod(f="source_id", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "integer")
        stop("invalid input: must be of class integer") 
    x@source_id <- value
    return(x)
})

# assay_type
setMethod(f="assay_type", signature="bioassaySet", definition=function(x) {return(x@assay_type)})
setReplaceMethod(f="assay_type", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "character")
        stop("invalid input: must be of class character") 
    x@assay_type <- value
    return(x)
})

# organism
setMethod(f="organism", signature="bioassaySet", definition=function(x) {return(x@organism)})
setReplaceMethod(f="organism", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "character")
        stop("invalid input: must be of class character") 
    x@organism <- value
    return(x)
})

# scoring
setMethod(f="scoring", signature="bioassaySet", definition=function(x) {return(x@scoring)})
setReplaceMethod(f="scoring", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "character")
        stop("invalid input: must be of class character") 
    x@scoring <- value
    return(x)
})

# target_types
setMethod(f="target_types", signature="bioassaySet", definition=function(x) {return(x@target_types)})
setReplaceMethod(f="target_types", signature="bioassaySet", definition=function(x, value) {
    if(class(value) != "character")
        stop("invalid input: must be of class character") 
    x@target_types <- value
    return(x)
})