trinarySimilarity <- function(queryMatrix, targetMatrix, minSharedScreenedTargets = 12, minSharedActiveTargets = 3){
    if(class(queryMatrix) != "dgCMatrix")
        stop("'queryMatrix' not of class 'dgCMatrix' as created by perTargetMatrix, if subsetting use drop=F")
    if(class(targetMatrix) != "dgCMatrix")
        stop("'targetMatrix' not of class 'dgCMatrix' as created by perTargetMatrix")
    if(ncol(queryMatrix) > 1)
        stop("'queryMatrix' has more than one column (compound)- subset first i.e. queryMatrix[,1,drop=F]")
    if(class(minSharedScreenedTargets) != "numeric")
        stop("'minSharedScreenedTargets' not of class 'numeric'")
    if(class(minSharedActiveTargets) != "numeric")
        stop("'minSharedActiveTargets' not of class 'numeric'")
    
    queryActives <- row.names(queryMatrix)[queryMatrix@i[queryMatrix@x == 2] + 1]
    queryInactives <- row.names(queryMatrix)[queryMatrix@i[queryMatrix@x == 1] + 1]    
    tMatrix <- as(targetMatrix, "TsparseMatrix")
    
    scores <- tapply(1:length(tMatrix@i), tMatrix@j, function(x){
        targetList <- tMatrix@i[x]
        targetScores <- tMatrix@x[x]
        targetActives <- row.names(tMatrix)[targetList[targetScores == 2] + 1]
        targetInactives <- row.names(tMatrix)[targetList[targetScores == 1] + 1]    
        intersectSize <- length(intersect(queryActives, targetActives)) + length(intersect(queryInactives, targetInactives))
        unionSize <- length(intersect(c(queryActives, queryInactives), c(targetActives, targetInactives)))
        if(length(intersect(queryActives, targetActives)) < minSharedActiveTargets && unionSize < minSharedScreenedTargets)
            return(NA)
        return(intersectSize / unionSize)
    }, simplify = T)
    scoreLabels <- as.numeric(names(scores))
    scores <- as.numeric(scores)
    names(scores) <- colnames(targetMatrix)[scoreLabels + 1]
    return(scores)
}