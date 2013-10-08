# takes a cid and returns a table of the proteins it shows activity against
activeTargets <- function(database, cid) {
    if(class(database) != "BioAssaySet")
        stop("database not of class BioAssaySet")
    if(is.numeric(cid)){
        cid <- as.character(cid)
    } else if(! is.character(cid)){
        stop("input not class character")
    }
    if(length(cid) > 1){
        warning("too many inputs, only the first was kept")
        cid <- cid[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", cid, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")

    queryResult <- .targetsByCid(database, cid) 
    if(nrow(queryResult) == 0){ return (NA) }
    queryByTarget <- split(queryResult$'activity', queryResult$'target')
    resultByTarget <- t(sapply(queryByTarget, function(x){
        x <- x[! is.na(x)]
        c(sum(x)/length(x), length(x))
    }))
    resultByTarget <- resultByTarget[resultByTarget[,1] > 0,,drop=F]
    if(nrow(resultByTarget) == 0){ return (NA) }
    resultByTarget <- resultByTarget[!is.na(rownames(resultByTarget)),,drop=F]
    if(nrow(resultByTarget) == 0){ return (NA) }
    colnames(resultByTarget) <- c("fraction_active", "total_screens")
    return(resultByTarget)
}

# takes a target id and returns a table of compounds that show activity against it
activeAgainst <- function(database, target){
    if(class(database) != "BioAssaySet")
        stop("database not of class BioAssaySet")
    if(is.numeric(target)){
        target <- as.character(target)
    } else if(! is.character(target)){
        stop("input not class character")
    }
    if(length(target) > 1){
        warning("too many inputs, only the first was kept")
        target <- target[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", target, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")

    queryResult <- .cidsByTarget(database, target)
    if(nrow(queryResult) == 0){ return(NA) }
    queryByCid <- split(queryResult$'activity', queryResult$'cid')
    resultByCid <- t(sapply(queryByCid, function(x){
        x <- x[! is.na(x)]
        c(sum(x)/length(x), length(x))
    }))
    resultByCid <- resultByCid[resultByCid[,1] > 0,,drop=F]
    if(nrow(resultByCid) == 0){ return (NA) }
    resultByCid <- resultByCid[!is.na(rownames(resultByCid)),,drop=F]
    if(nrow(resultByCid) == 0){ return (NA) }
    colnames(resultByCid) <- c("fraction_active", "total_assays")
    return(resultByCid)
}

# takes a target id and returns compounds ordered by decreasing selectivity for
# this target
selectiveAgainst <- function(database, target, maxCompounds = 10, minimumTargets = 10){
    if(class(database) != "BioAssaySet")
        stop("database not of class BioAssaySet")
    if(is.numeric(target)){
        target <- as.character(target)
    } else if(! is.character(target)){
        stop("input not class character")
    }
    if(length(target) > 1){
        warning("too many inputs, only the first was kept")
        target <- target[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", target, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")
    if(! is.numeric(maxCompounds))
        stop("non numeric maxCompounds input")
    if(! is.numeric(minimumTargets))
        stop("non numeric minimumTargets input")

    activeCids <- unique(.cidsByTarget(database, target, activeOnly = TRUE)$cid)
    if(length(activeCids) == 0){ return(NA) }
    targetSums <- t(sapply(activeCids, function(cid){
        queryResult <- .targetsByCid(database, cid)
        queryByCid <- split(queryResult$'activity', queryResult$'target')
        resultByCid <- t(sapply(queryByCid, function(x){
            x <- x[! is.na(x)]
            c(sum(x), length(x))
        }))
        c(sum(resultByCid[,1] > 0), length(rownames(resultByCid)))
    }))
    if(dim(targetSums)[2] == 0){ return(NA) }
    row.names(targetSums) <- activeCids
    colnames(targetSums) <- c("active_targets", "tested_targets")
    targetSums <- targetSums[targetSums[,2] >= minimumTargets,,drop=F]
    if(dim(targetSums)[1] == 0){ return(NA) }
    targetSums <- targetSums[order(targetSums[,1], targetSums[,2], na.last = NA),,drop=F]
    if(dim(targetSums)[1] > maxCompounds){
        targetSums <- targetSums[1:maxCompounds,]
    }
    return(targetSums)
}

# this takes an assay id (aid) and returns a bioassay object for the corresponding assay
getAssay <- function(database, aid){
    if(class(database) != "BioAssaySet")
        stop("database not of class BioAssaySet")
    if(is.numeric(aid)){
        aid <- as.character(aid)
    } else if(! is.character(aid)){
        stop("input not class character")
    }
    if(length(aid) > 1){
        warning("too many inputs, only the first was kept")
        aid <- aid[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", aid, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")

    aid <- as.character(aid)
    assay <- queryBioassayDB(database, paste("SELECT * FROM assays WHERE aid =", aid))
    if(nrow(assay) < 1){
        stop("assay not found")
    }
    source_id <- queryBioassayDB(database, paste("SELECT description FROM assays NATURAL JOIN sources WHERE aid =", aid))$description
    target_list <- queryBioassayDB(database, paste("SELECT * FROM targets WHERE aid =", aid))
    targets <- target_list$target
    target_types <- target_list$target_type
    scores <- queryBioassayDB(database, paste("SELECT cid, sid, activity, score FROM activity WHERE aid =", aid))

    new("bioassay",
        aid = aid,
        source_id = source_id,
        assay_type = assay$assay_type,
        organism = assay$organism,
        scoring = assay$scoring,
        targets = targets,
        target_types = target_types,
        scores = scores
    )
}

# this returns a matrix with activity values for each compound in the database
activityMatrix <- function(database, maxAssayLimit=100){
    if(class(database) != "BioAssaySet")
        stop("database not of class BioAssaySet")
    if(! is.numeric(maxAssayLimit))
        stop("maxAssayLimit not numeric")

    totalAssays <- queryBioassayDB(database, "SELECT COUNT(DISTINCT aid) FROM assays")[[1]]
    if(totalAssays > maxAssayLimit){
        stop("Database has too many assays. Try a larger maxAssayLimit or smaller database.")
    }
    result <- queryBioassayDB(database, "SELECT cid, aid, activity FROM activity")
    cids <- unique(result$cid)
    aids <- unique(result$aid)
    fullMatrix <- matrix(nrow=length(cids), ncol=length(aids), dimnames=list(cids, aids))
    splitResult <- split(result, result$cid)
    for(x in splitResult){
        fullMatrix[as.character(x$cid[1]),as.character(x$aid)] <- x$activity
    }
    return(fullMatrix)
}

###########################
# index optimized queries #
###########################

# Performance: indexes activity_cid, and targets_aid provide ~12x speedup
.targetsByCid <- function(database, cid){
    queryBioassayDB(database, paste(
        "SELECT activity.aid, activity.activity, targets.target",
        " FROM activity",
        " NATURAL JOIN targets",
        " WHERE activity.cid = '", cid, "'",
        " AND activity.activity NOT NULL",
        " AND targets.target_type = 'protein'",
        sep = ""
    ))
}

# Performance: indexes targets_target, activity_aid provide MASSIVE speedup
.cidsByTarget <- function(database, target, activeOnly = FALSE){
    if(activeOnly){
        queryString <- " AND activity.activity = 1"
    } else {
        queryString <- ""
    }
    queryBioassayDB(database, paste(
        "SELECT activity.aid, activity.activity, activity.cid",
        " FROM targets",
        " NATURAL JOIN activity",
        " WHERE targets.target = '", target, "'",
        queryString,
        sep = ""
    ))
}

# # Performance: uses indexes activity_cid, targets_aid
# .countActiveTargets <- function(database, cid) {
#     queryBioassayDB(database, paste(
#         "SELECT COUNT(DISTINCT targets.target)",
#         " FROM activity",
#         " NATURAL JOIN targets",
#         " WHERE activity.cid = '", cid, "'",
#         " AND activity.activity = '1'",
#         " AND targets.target_type = 'protein'",
#         sep = ""
#     ))
# }
