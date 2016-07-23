# centers and standardizes the numeric activity scores for
# a bioassaySet object
# only consideres '0' values actually encoded in the i,j,x representation
# of the sparse matrix- omits empty '0' entries
# these can be computationally distinguished because the sparse matrix actually returns
# zero values it was loaded with.
scaleBioassaySet <- function(bioassaySet, center=TRUE, scale=TRUE){
    if(class(bioassaySet) != "bioassaySet")
        stop("input not of class bioassaySet")    
    if(class(center) != "logical")
        stop("center option not of class logical (TRUE or FALSE)")
    if(class(scale) != "logical")
        stop("scale option not of class logical (TRUE or FALSE)")
    
    tMatrix <- as(bioassaySet@scores, "TsparseMatrix")
    sortIndex <- sort(tMatrix@i, decreasing=F, index.return=T)$ix
    i <- tMatrix@i[sortIndex]
    j <- tMatrix@j[sortIndex]
    x <- tMatrix@x[sortIndex]
    
    standardScores <- tapply(x, i, scale, center = center, scale = scale, simplify=FALSE)
    standardScores <- do.call(c, standardScores)
    
    bioassaySet@scores <- sparseMatrix(
        i = i,
        j = j,
        x = standardScores,
        dims = tMatrix@Dim,
        dimnames = tMatrix@Dimnames,
        symmetric = FALSE,
        index1 = FALSE)
    
    return(bioassaySet)
}

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
# values if useNumericScores = FALSE:
#   0 = untested or inconclusive
#   1 = inactive
#   2 = active in 1 or more assays (can still be inactive in some)
perTargetMatrix <- function(assays, inactives = TRUE, assayTargets = FALSE, targetOrder = FALSE, summarizeReplicates = "activesFirst", useNumericScores = FALSE){
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
    if(class(summarizeReplicates) != "function" && class(summarizeReplicates) != "standardGeneric")
        if(summarizeReplicates != "mode" && summarizeReplicates != "activesFirst")
            stop("invalid conflict resolver option")
    if(class(useNumericScores) != "logical")
        stop("useNumericScores option not of class logical (TRUE or FALSE)")
    
    # message("Note: in this version active scores now use a 2 instead of a 1")

    if(useNumericScores){
        # assays <- scaleBioassaySet(getAssays(db, c(348, 360, 368))) # test code
        activityMatrix <- assays@scores
        tMatrix <- as(activityMatrix, "TsparseMatrix")
        i <- tMatrix@i
        j <- tMatrix@j
        x <- tMatrix@x
        x[x == "NaN"] <- 0
        sortIndex <- sort.int(abs(x), decreasing=T, index.return=T)$ix
        i <- i[sortIndex]
        j <- j[sortIndex]
        x <- x[sortIndex]
        coords <- cbind(i+1, j+1, x)
    } else {
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
    
    if(class(summarizeReplicates) == "character")
        if(summarizeReplicates == "mode")
            summarizeReplicates <- function(x) { as.numeric(names(which.max(table(x)))) }
    
    if(class(summarizeReplicates) != "character"){

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
targetSelectivity <- function(database, cids, scoring="total", category=FALSE, multiTarget="keepOne"){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    if(is.numeric(cids)){
        cids <- as.character(cids)
    } else if(! is.character(cids)){
        stop("'cids' not class 'character'")
    }
    if(! scoring %in% c("total", "fraction")){
        stop("'scoring' must be one of the following: total or fraction")
    }
    allCategories <- .allCategories(database) 
    if(is.logical(category) && category)
        stop(paste(c("'category' must be either 'FALSE' or an available translation category.
             valid options: FALSE", allCategories), collapse=" "))
    if(! is.logical(category) && (! category %in% allCategories))
        stop(paste(c("'category' must be either 'FALSE' or an available translation category.
             valid options: FALSE", allCategories), collapse=" "))
    if((multiTarget != "drop") && (multiTarget != "keepOne") && (multiTarget != "all"))
        stop("multiTarget must be one of the following: drop, keepOne, all")
    
    return(sapply(cids, function(cid){
        activityData <- .targetsByCid(database, cid)
        if(multiTarget == "drop"){
            duplicatedAids <- activityData$aid[duplicated(activityData$aid)]
            activityData <- activityData[! activityData$aid %in% duplicatedAids,,drop=F]
        }
        if(multiTarget == "keepOne"){
            duplicatedRows <- duplicated(activityData$aid)
            dropAidTarget <-  paste(activityData$aid[duplicatedRows], activityData$target[duplicatedRows], sep="______") 
            activityData <- activityData[! duplicatedRows,,drop=F]    
        }
        if(! is.logical(category)){
            categoryTargets <- .targetsByCid(database, cid, category=category)
            if(multiTarget == "drop"){
                categoryTargets <- categoryTargets[! categoryTargets$aid %in% duplicatedAids,]
            }
            if(multiTarget == "keepOne"){
                aidTargetByRow <-  paste(categoryTargets$aid, categoryTargets$target, sep="______")   
                categoryTargets <- categoryTargets[! aidTargetByRow %in% dropAidTarget,,drop=F]
            }
            categoryTargets <- categoryTargets[sort(categoryTargets$activity, decreasing=TRUE, index.return=T)$ix,,drop=F]
            categoryTargets <- categoryTargets[! duplicated(categoryTargets$target),,drop=F]
            categoryDuplicates <- unique(categoryTargets$target[duplicated(categoryTargets$identifier)])
            activityData <- activityData[! activityData$target %in% categoryDuplicates,,drop=F]
        }
        activeTargetIds <- unique(activityData$target[activityData$activity == 1])
        
        if(scoring == "total"){
            return(length(activeTargetIds))
        } else {
            inactiveTargetIds <- unique(activityData$target[activityData$activity == 0])
            inactiveTargetIds <- inactiveTargetIds[! inactiveTargetIds %in% activeTargetIds]
            return(length(activeTargetIds)/(length(activeTargetIds)+length(inactiveTargetIds)))
        }
    }))  
}

# takes a minimum number of distinct protein targets
# and returns all CIDs screened against at least that many proteins
screenedAtLeast <- function(database, minTargets, inconclusives=TRUE){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    if(class(minTargets) != "numeric")
        stop("'minTargets' not of class 'numeric'")
    if(class(inconclusives) != "logical")
        stop("inconclusives option not of class logical (TRUE or FALSE)")
    counts <- .screenedAtLeast(database, minTargets, inconclusives)
    return(counts$cid[! is.na(counts$cid)])
}

# takes a target identifier and category and returns its translation
# to another target type
translateTargetId <- function(database, target, category, fromCategory="GI"){
    if(class(database) != "BioassayDB")
        stop("'database' not of class 'BioassayDB'")
    allCategories <- c("GI", .allCategories(database))
    if(! category %in% allCategories)
        stop(paste(c("'category' must be either an available translation category.
             valid options: ", allCategories), collapse=" "))
    if(! fromCategory %in% allCategories)
        stop(paste(c("'fromCategory' must be either an available translation category.
             valid options: ", allCategories), collapse=" "))
    if(length(target) > 1)
        stop("'target' has length greater than 1, use lapply to run multiple queries")
    
    if(fromCategory != "GI")
        targetGIs <- unique(.reverseTranslateTargetId(database, target, fromCategory))
    else
        targetGIs <- target

    if(category == "GI")
        result <- targetGIs
    else
        result <- unique(unlist(lapply(targetGIs, .translateTargetId, database=database, category=category)))
    
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
        sep=""))[[1]]
}

.reverseTranslateTargetId <- function(database, identifier, category){
    queryBioassayDB(database, paste(
        "SELECT target",
        " FROM targetTranslations",
        " WHERE identifier = '", identifier, "'",
        " AND category = '", category, "'",
        sep=""))[[1]]
}

.screenedAtLeast <- function(database, minTargets, inconclusives){
    if(inconclusives){
        queryString <- paste(
            "SELECT cid, COUNT(DISTINCT(target)) AS distinctTargets",
            " FROM activity",
            " NATURAL JOIN targets",
            " GROUP BY cid",
            " HAVING distinctTargets >= ", minTargets, sep="")        
    } else {
        queryString <- paste(
            "SELECT cid, COUNT(DISTINCT(target)) AS distinctTargets",
            " FROM activity",
            " NATURAL JOIN targets",
            " WHERE activity NOT NULL",
            " GROUP BY cid",
            " HAVING distinctTargets >= ", minTargets, sep="")
    }
    queryBioassayDB(database, queryString)
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

.allCategories <- function(database){
    queryBioassayDB(database, "SELECT DISTINCT category FROM targetTranslations")[[1]]
}

# Performance: indexes activity_cid, and targets_aid provide ~12x speedup
.targetsByCid <- function(database, cid, activity = "NOT NULL", category=FALSE){
    if(activity != "NOT NULL"){
        activity <- paste("= '", activity, "'", sep="")
    }
    if(! is.logical(category)){
        sql <- paste(
            "SELECT activity.aid, activity.activity, targets.target, targetTranslations.identifier",
            " FROM activity",
            " NATURAL JOIN targets",
            " NATURAL JOIN targetTranslations",
            " WHERE activity.cid = '", cid, "'",
            " AND activity.activity ", activity,
            " AND targetTranslations.category = '", category, "'",
            sep = ""
        )
    } else {
        sql <- paste(
            "SELECT activity.aid, activity.activity, targets.target",
            " FROM activity",
            " NATURAL JOIN targets",
            " WHERE activity.cid = '", cid, "'",
            " AND activity.activity ", activity,
            sep = ""
        )
    }
    queryBioassayDB(database, sql)
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
