## Define print behavior for bioassay
setMethod("show", signature=signature(
    object="bioassay"),
    function(object) {
        cat("class:\t\t", class(object), "\n")
        cat("aid:\t\t", slot(object, "aid"), "\n")
        cat("source_id:\t", slot(object, "source_id"), "\n")
        cat("assay_type:\t", slot(object, "assay_type"), "\n")
        cat("organism:\t", slot(object, "organism"), "\n")
        cat("scoring:\t", slot(object, "scoring"), "\n")
        cat("targets:\t", slot(object, "targets"), "\n")
        cat("target_types:\t", slot(object, "target_types"), "\n")
        cat("total scores:\t", dim(slot(object, "scores"))[1], "\n")
    }
)

## methods for bioassay
setMethod(f="aid", signature="bioassay", definition=function(x) {return(x@aid)})
setReplaceMethod(f="aid", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    if(length(value) > 1){
        warning("too many inputs, only the first was kept")
        value <- value[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", value, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace") 

    x@aid <- value
    return(x)
})

setMethod(f="source_id", signature="bioassay", definition=function(x) {return(x@source_id)})
setReplaceMethod(f="source_id", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    if(length(value) > 1){
        warning("too many inputs, only the first was kept")
        value <- value[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", value, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace") 

    x@source_id <- value
    return(x)
})

setMethod(f="assay_type", signature="bioassay", definition=function(x) {return(x@assay_type)})
setReplaceMethod(f="assay_type", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    if(length(value) > 1){
        warning("too many inputs, only the first was kept")
        value <- value[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", value, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace") 

    x@assay_type <- value
    return(x)
})

setMethod(f="organism", signature="bioassay", definition=function(x) {return(x@organism)})
setReplaceMethod(f="organism", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    if(length(value) > 1){
        warning("too many inputs, only the first was kept")
        value <- value[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", value, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace") 

    x@organism <- value
    return(x)
})

setMethod(f="scoring", signature="bioassay", definition=function(x) {return(x@scoring)})
setReplaceMethod(f="scoring", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    if(length(value) > 1){
        warning("too many inputs, only the first was kept")
        value <- value[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", value, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace") 

    x@scoring <- value
    return(x)
})

setMethod(f="targets", signature="bioassay", definition=function(x) {return(x@targets)})
setReplaceMethod(f="targets", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    for(i in value){
        if(! grepl("^[a-zA-Z_0-9\\s]+$", i, perl=T))
            stop("invalid input: must contain only alphanumerics and/or whitespace") 
    }

    x@targets <- value
    return(x)
})

setMethod(f="target_types", signature="bioassay", definition=function(x) {return(x@target_types)})
setReplaceMethod(f="target_types", signature="bioassay", definition=function(x, value) {
    if(is.numeric(value)){
        value <- as.character(value)
    } else if(! is.character(value)){
        stop("input not class character")
    }
    for(i in value){
        if(! grepl("^[a-zA-Z_0-9\\s]+$", i, perl=T))
            stop("invalid input: must contain only alphanumerics and/or whitespace") 
    }

    x@target_types <- value
    return(x)
})

setMethod(f="scores", signature="bioassay", definition=function(x) {return(x@scores)})
setReplaceMethod(f="scores", signature="bioassay", definition=function(x, value) {
    if(! is.data.frame(value)){
        stop("input not data frame")
    }
    if(ncol(value) != 4){
        stop("scores must have 4 columns")
    }

    x@scores <- value
    return(x)
})
