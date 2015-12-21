# returns a ChemmineR style FPset object for a given set of cids and targets
bioactivityFingerprint <- function(bioassaySet, targets = FALSE, summarizeReplicates = "activesFirst"){
  if(class(bioassaySet) != "bioassaySet")
    stop("input not of class bioassaySet")
  if(class(summarizeReplicates) != "function" && summarizeReplicates != "mode" && summarizeReplicates != "activesFirst")
    stop("invalid conflict resolver option")
  if(targets){
      if(length(targets) < 1){
        stop("target list must have at least one valid entry")  
      }
      if(is.numeric(targets)){
        targets <- as.character(targets)
      } else if(! is.character(targets)){
        stop("targets not class numeric or character")
      } 
  } else {
      targets <- allTargets(bioassaySet)
  }
  activityMatrix <- perTargetMatrix(bioassaySet, inactives = TRUE, targetOrder = targets, summarizeReplicates = summarizeReplicates)
  binaryMatrix <- t(as.matrix(1*(activityMatrix > 1)))
  fingerPrint <- as(binaryMatrix, "FPset")
  slot(fingerPrint, "type") <- "bioactivity"
  return(fingerPrint)
}

# lists all CIDS present in a BioassayDB, bioassay, bioassaySet, or target matrix (dgCMatrix) object
allCids <- function(inputObject, activesOnly = FALSE){
    if(class(activesOnly) != "logical")
        stop("activesOnly not of class logical (TRUE / FALSE)")
    if(class(inputObject) == "BioassayDB"){
        return(.allCids(inputObject, activesOnly)[[1]])
    } else if(class(inputObject) == "bioassay"){
        if(activesOnly)
            return(unique(slot(inputObject, "scores")$cid[slot(inputObject, "scores")$activity == 1]))
        else
            return(unique(slot(inputObject, "scores")$cid))
    } else if(class(inputObject) == "bioassaySet"){
        if(activesOnly)
            return(unique(colnames(slot(inputObject, "activity"))[which(slot(inputObject, "activity") == 2, arr.ind=T)[,2]]))
        else
            return(unique(colnames(slot(inputObject, "activity"))))
    } else if(class(inputObject) == "dgCMatrix"){
        if(activesOnly)
            return(unique(colnames(inputObject)[which(inputObject == 2, arr.ind=T)[,2]]))
        else
            return(unique(colnames(inputObject)))
    } else
        stop("inputObject not of BioassayDB, bioassay, bioassaySet, or target matrix (dgCMatrix) class")
}

# lists all Target IDs present in a BioassayDB, bioassay, bioassaySet, or target matrix (dgCMatrix) object
allTargets <- function(inputObject){
    if(class(inputObject) == "BioassayDB"){
        return(.allTargets(inputObject)[[1]])
    } else if(class(inputObject) == "bioassay"){
        return(unique(slot(inputObject, "targets")))
    } else if(class(inputObject) == "bioassaySet"){
        return(unique(colnames(slot(inputObject, "targets"))))
    } else if(class(inputObject) == "dgCMatrix"){
        return(unique(rownames(inputObject)))
    } else
        stop("inputObject not of BioassayDB, bioassay, bioassaySet, or target matrix (dgCMatrix) class")
}

# takes a cid and returns a table of the targets it shows activity against
activeTargets <- function(database, cid){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
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

# takes a cid and returns a table of the targets it has been found inactive against 
inactiveTargets <- function(database, cid){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
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

    queryResult <- .targetsByCid(database, cid, activity = 0) 
    if(nrow(queryResult) == 0){ return (NA) }
    inactiveTargets <- table(queryResult[,3])
    inactiveTargets <- as.data.frame(inactiveTargets, row.names=(names(inactiveTargets)))[,2, drop=FALSE]
    colnames(inactiveTargets) <- "inactive screens"
    return(inactiveTargets)
}

# takes a target id and returns a table of compounds that show activity against it
activeAgainst <- function(database, target){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
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
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
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

# This takes a bioassaySet of multiple assays and returns
# a vector of the targets of each, with the assay identifiers themselves (aids)
# as names. If a single assay contains multiple targets, these will all be listed.
assaySetTargets <- function(assays){
    # check input sanity
    if(class(assays) != "bioassaySet")
        stop("database not of class bioassaySet")
    targetMatrix <- slot(assays, "targets")
    targetCoords <- which(targetMatrix == 1, arr.ind=TRUE)
    assayTargets <- colnames(targetMatrix)[targetCoords[,2]]
    names(assayTargets) <- rownames(targetMatrix)[targetCoords[,1]]
    return(assayTargets)
}

# This takes a bioassaySet of multiple assays, and returns one with a single score
# per target- collapsing multiple assays against distinct targets into a single assay
# values:
#   0 = untested or inconclusive
#   1 = inactive
#   2 = active in 1 or more assays (can still be inactive in some)
perTargetMatrix <- function(assays, inactives = TRUE, assayTargets = FALSE, targetOrder = FALSE, summarizeReplicates = "activesFirst"){
    # check input sanity
    if(class(assays) != "bioassaySet")
        stop("database not of class bioassaySet")
    if(class(inactives) != "logical")
        stop("inactives option not of class logical (TRUE or FALSE)")
    if(! is.logical(assayTargets) && (class(assayTargets) != "character"))
        stop("assayTargets not of class character")
    if(is.logical(targetOrder) && isTRUE(targetOrder))
        stop("TRUE is an invalid option for targetOrder- it must be FALSE or character")
    if(! is.logical(targetOrder) && (class(targetOrder) != "character"))
        stop("targetOrder not of class character")
    if(is.logical(targetOrder) && isTRUE(targetOrder))
        stop("TRUE is an invalid option for targetOrder- it must be FALSE or character")
    if(class(summarizeReplicates) != "function" && summarizeReplicates != "mode" && summarizeReplicates != "activesFirst")
        stop("invalid conflict resolver option")

    message("Note: in this version active scores now use a 2 instead of a 1")

    # Get coordinates of scores
    # NOTE: inactives must come after active values
    #   for later duplicate removal to give preference to active scores
    activityMatrix <- slot(assays, "activity")
    coords <- which(activityMatrix == 2, arr.ind=TRUE)
    suppressWarnings(coords <- cbind(coords,2))
    if(inactives){
        inactiveCoords <- which(activityMatrix == 1, arr.ind=TRUE)
        suppressWarnings(inactiveCoords <- cbind(inactiveCoords,1))
        coords <- rbind(coords, inactiveCoords)
    }
    
    # get assay to target mappings 
    if(is.logical(assayTargets)){
        assayTargets <- assaySetTargets(assays)
    }

    # get target/cid pairs per assay and drop those without targets
    targetsPerAssay <- assayTargets[as.character(rownames(activityMatrix)[coords[,1]])]
    cidsPerAssay <- colnames(activityMatrix)[coords[,2]]
    cidsPerAssay <- cidsPerAssay[! is.na(targetsPerAssay)]
    activityScores <- coords[,3][! is.na(targetsPerAssay)]
    targetsPerAssay <- targetsPerAssay[! is.na(targetsPerAssay)]    

    if(summarizeReplicates != "activesFirst"){
        if(summarizeReplicates == "mode")
            summarizeReplicates <- function(x) { as.numeric(names(which.max(table(x)))) }

        # identify cid/target pairs that are duplicated and resolve according to specified rule
        mergedPairs <- paste(cidsPerAssay, targetsPerAssay, sep="______")
        duplicatedPairs <- duplicated(mergedPairs) | duplicated(mergedPairs, fromLast=TRUE)
        if(sum(duplicatedPairs) > 0){
            pairFactor <- factor(mergedPairs[duplicatedPairs], levels=unique(mergedPairs[duplicatedPairs]))
            splitScores <- split(activityScores[duplicatedPairs], pairFactor)
            resolvedScores <- sapply(splitScores, summarizeReplicates)
            firstDuplicates <- duplicated(mergedPairs, fromLast=TRUE) & !duplicated(mergedPairs)
            activityScores[firstDuplicates] <- resolvedScores
        }
    }

    # make cid/target pairs unique
    nonDuplicates <- ! duplicated(paste(cidsPerAssay, targetsPerAssay, sep="______"))
    cidsPerAssay <- cidsPerAssay[nonDuplicates]
    targetsPerAssay <- targetsPerAssay[nonDuplicates]
    activityScores <- activityScores[nonDuplicates]

    # get unique values for cols and rows
    if(is.logical(targetOrder)){
        uniqueTargets <- unique(targetsPerAssay)
    } else {
        uniqueTargets <- targetOrder
        inTargetList <- targetsPerAssay %in% targetOrder
        targetsPerAssay <- targetsPerAssay[inTargetList]
        cidsPerAssay <- cidsPerAssay[inTargetList]
        activityScores <- activityScores[inTargetList]
    }
    uniqueCids <- unique(cidsPerAssay)
    
    # return matrix
    sparseMatrix(
        i = match(targetsPerAssay, uniqueTargets),
        j = match(cidsPerAssay, uniqueCids),
        x = activityScores,
        dims = c(length(uniqueTargets), length(uniqueCids)),
        dimnames = list(uniqueTargets, uniqueCids),
        symmetric = FALSE,
        index1 = TRUE)
}

# this takes a list of compound ids (cids) and returns a bioassaySet object
# with only the activity of those compounds listed
getBioassaySetByCids <- function(database, cids){
    # check input sanity
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(is.numeric(cids)){
        cids <- as.character(cids)
    } else if(! is.character(cids)){
        stop("cids not class numeric or character")
    }
    if(length(unique(cids)) != length(cids)){
        stop("cid list contains duplicates")
    }

    # create activity and score matrix
    result <- .activityByCids(database, cids)
    uniqueAssays <- unique(result$aid)
    uniqueCids <- unique(result$cid)
    activity <- sparseMatrix(
        i = match(result$aid, uniqueAssays),
        j = match(result$cid, uniqueCids),
        x = result$activity + 1,
        dims = c(length(uniqueAssays), length(uniqueCids)),
        dimnames = list(uniqueAssays, uniqueCids),
        symmetric = FALSE,
        index1 = TRUE)
    scores <- sparseMatrix(
        i = match(result$aid, uniqueAssays),
        j = match(result$cid, uniqueCids),
        x = result$score,
        dims = c(length(uniqueAssays), length(uniqueCids)),
        dimnames = list(uniqueAssays, uniqueCids),
        symmetric = FALSE,
        index1 = TRUE)
    aids <- uniqueAssays

    # create target matrix
    resultTargets <- .targetsByAids(database, aids)
    uniqueTargets <- unique(resultTargets$target)
    targets <- sparseMatrix(
        i = match(resultTargets$aid, uniqueAssays),
        j = match(resultTargets$target, uniqueTargets),
        x = 1,
        dims = c(length(uniqueAssays), length(uniqueTargets)),
        dimnames = list(uniqueAssays, uniqueTargets),
        symmetric = FALSE,
        index1 = TRUE)
    target_types <- resultTargets$target_type
    names(target_types) <- resultTargets$aid

    # get data sources
    resultSource <- .sourcesByAids(database, aids)

    # get assay annotation
    resultAnnotation <- .assaysByAids(database, aids)
    source_id <- resultAnnotation$source_id
    names(source_id) <- resultAnnotation$aid
    assay_type <- resultAnnotation$assay_type
    names(assay_type) <- resultAnnotation$aid
    organism <- resultAnnotation$organism
    names(organism) <- resultAnnotation$aid
    scoring <- resultAnnotation$scoring
    names(scoring) <- resultAnnotation$aid 

    # create bioassaySet object
    new("bioassaySet",
        activity = activity,
        scores = scores,
        targets = targets,
        sources = resultSource,
        source_id = source_id,
        assay_type = assay_type,
        organism = organism,
        scoring = scoring,
        target_types = target_types)
}

# this takes a list of assay ids (aids) and returns a bioassaySet object
getAssays <- function(database, aids){
    # check input sanity
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(is.numeric(aids)){
        aids <- as.character(aids)
    } else if(! is.character(aids)){
        stop("aids not class numeric or character")
    }

    # create activity and score matrix
    result <- .activityByAids(database, aids)
    uniqueAssays <- unique(result$aid)
    uniqueCids <- unique(result$cid)
    activity <- sparseMatrix(
        i = match(result$aid, uniqueAssays),
        j = match(result$cid, uniqueCids),
        x = result$activity + 1,
        dims = c(length(uniqueAssays), length(uniqueCids)),
        dimnames = list(uniqueAssays, uniqueCids),
        symmetric = FALSE,
        index1 = TRUE)
    scores <- sparseMatrix(
        i = match(result$aid, uniqueAssays),
        j = match(result$cid, uniqueCids),
        x = result$score,
        dims = c(length(uniqueAssays), length(uniqueCids)),
        dimnames = list(uniqueAssays, uniqueCids),
        symmetric = FALSE,
        index1 = TRUE)

    # create target matrix
    resultTargets <- .targetsByAids(database, aids)
    uniqueTargets <- unique(resultTargets$target)
    targets <- sparseMatrix(
        i = match(resultTargets$aid, uniqueAssays),
        j = match(resultTargets$target, uniqueTargets),
        x = 1,
        dims = c(length(uniqueAssays), length(uniqueTargets)),
        dimnames = list(uniqueAssays, uniqueTargets),
        symmetric = FALSE,
        index1 = TRUE)
    target_types <- resultTargets$target_type
    names(target_types) <- resultTargets$aid

    # get data sources
    resultSource <- .sourcesByAids(database, aids)

    # get assay annotation
    resultAnnotation <- .assaysByAids(database, aids)
    source_id <- resultAnnotation$source_id
    names(source_id) <- resultAnnotation$aid
    assay_type <- resultAnnotation$assay_type
    names(assay_type) <- resultAnnotation$aid
    organism <- resultAnnotation$organism
    names(organism) <- resultAnnotation$aid
    scoring <- resultAnnotation$scoring
    names(scoring) <- resultAnnotation$aid 

    # create bioassaySet object
    new("bioassaySet",
        activity = activity,
        scores = scores,
        targets = targets,
        sources = resultSource,
        source_id = source_id,
        assay_type = assay_type,
        organism = organism,
        scoring = scoring,
        target_types = target_types)
}

# this takes an assay id (aid) and returns a bioassay object for the corresponding assay
getAssay <- function(database, aid){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
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
    scores <- queryBioassayDB(database, paste("SELECT cid, activity, score FROM activity WHERE aid =", aid))

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

# takes cids and a database and returns the target selectivity of the query cids
# scoring method is either as total active targets or fraction of active targets
targetSelectivity <- function(database, cids, scoring="total"){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    if(is.numeric(cids)){
        cids <- as.character(cids)
    } else if(! is.character(cids)){
        stop("'cids' not class 'character'")
    }
    if(! scoring %in% c("total", "fraction")){
        stop("'scoring' must be of type 'total' or 'fraction'")
    }
    activeTargets <- sapply(cids, function(cid){
        length(unique(.targetsByCid(database, cid, activity = 1)$target))    
    })
    if(scoring == "total"){
        return(activeTargets)
    }
    inactiveTargets <- sapply(cids, function(cid){
        length(unique(.targetsByCid(database, cid, activity = 0)$target))    
    })
    return(activeTargets/(activeTargets+inactiveTargets))    
}

# takes a minimum number of distinct protein targets
# and returns all CIDs screened against at least that many proteins
screenedAtLeast <- function(database, minTargets){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    if(class(minTargets) != "numeric")
        stop("'minTargets' not of class 'numeric'")
    counts <- .screenedAtLeast(database, minTargets)
    return(counts$cid[! is.na(counts$cid)])
}

# takes a target identifier and category and returns its translation
# to another target type
translateTargetId <- function(database, target, category){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    result <- .translateTargetId(database, target, category)[[1]]
    if(length(result) == 0)
        result <- NA
    return(result)
}

###########
# queries #
###########

.allCids <- function(database, activesOnly = FALSE){
    if(activesOnly){
        queryString <- " AND activity == 1"
    } else {
        queryString <- ""
    }
    queryBioassayDB(database, paste(
        "SELECT DISTINCT cid FROM activity WHERE cid NOT NULL",
        queryString,
        sep = ""
    ))
}

.allTargets <- function(database){
    queryBioassayDB(database,
        "SELECT DISTINCT target FROM targets WHERE target NOT NULL")
}

.translateTargetId <- function(database, target, category){
    queryBioassayDB(database, paste(
        "SELECT identifier",
        " FROM targetTranslations",
        " WHERE target = '", target, "'",
        " AND category = '", category, "'",
        sep=""))
}

.screenedAtLeast <- function(database, minTargets){
    queryBioassayDB(database, paste(
        "SELECT cid, COUNT(DISTINCT(target)) AS distinctTargets",
        " FROM activity",
        " NATURAL JOIN targets",
        " GROUP BY cid",
        " HAVING distinctTargets >= ", minTargets, sep=""))
}

.targetsByAids  <- function(database, aids){
    con <- slot(database, "database")
    sql <- paste("SELECT * FROM targets WHERE aid = $AID")
    dbBegin(con)
    result <- dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aids))
    dbCommit(con)
    return(result)
}

.activityByAids <- function(database, aids){
    con <- slot(database, "database")
    sql <- paste("SELECT * FROM activity WHERE aid = $AID")
    dbBegin(con)
    result <- dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aids))
    dbCommit(con)
    return(result)
}

.activityByCids <- function(database, cids){
    con <- slot(database, "database")
    sql <- paste("SELECT * FROM activity WHERE cid = $CID")
    dbBegin(con)
    result <- dbGetPreparedQuery(con, sql, bind.data = data.frame(CID=cids))
    dbCommit(con)
    return(result)
}

.assaysByAids <- function(database, aids){
    con <- slot(database, "database")
    sql <- paste("SELECT * FROM assays WHERE aid = $AID")
    dbBegin(con)
    result <- dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aids))
    dbCommit(con)
    return(result)
}

.sourcesByAids <- function(database, aids){
    con <- slot(database, "database")
    sql <- paste("SELECT DISTINCT source_id, description, version FROM assays NATURAL JOIN sources WHERE aid = $AID")
    dbBegin(con)
    result <- dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aids))
    dbCommit(con)
    result[unique(result$source_id),]
}

# Performance: indexes activity_cid, and targets_aid provide ~12x speedup
.targetsByCid <- function(database, cid, activity = "NOT NULL"){
    if(activity != "NOT NULL"){
        activity <- paste("= '", activity, "'", sep="")
    }
    queryBioassayDB(database, paste(
        "SELECT activity.aid, activity.activity, targets.target",
        " FROM activity",
        " NATURAL JOIN targets",
        " WHERE activity.cid = '", cid, "'",
        " AND activity.activity ", activity,
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
