trinarySimilarity <- function(queryMatrix, targetMatrix, minSharedAssays = 12, minSharedActives = 3){
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
        if(length(intersect(queryActives, targetActives)) < minSharedActives && unionSize < minSharedAssays)
            return(NA)
        if(unionSize == 0)
            return(NA)
        return(intersectSize / unionSize)
    }, simplify = T)
    names(scores) <- colnames(targetMatrix)[as.numeric(names(scores)) + 1]
    return(scores)
}