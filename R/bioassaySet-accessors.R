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
