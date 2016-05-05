# computes a bayesian cross-reactivity score for every
# cid in a perTargetMatrix
# based on method published here:
# Dančík, V. et al. Connecting Small Molecules with Similar Assay Performance 
# Profiles Leads to New Biological Hypotheses. J Biomol Screen 19, 771-781 (2014).
# note that the prior hit_ratio_sd mentioned in the paper is actually the variance
crossReactivityProbability <- function(inputMatrix, 
                                       threshold=0.25,
                                       prior=list(hit_ratio_mean=0.126, hit_ratio_sd=0.1072)){
    if(class(inputMatrix) != "dgCMatrix")
        stop("'inputMatrix' not of class 'dgCMatrix'")
    if(class(threshold) != "numeric")
        stop("'threshold' not of class 'numeric'")
    if(threshold > 1 || threshold < 0)
        stop("'threshold' must be between 0 and 1")
    if(class(prior) != "list")
        stop("'prior' not of class 'list'")
    if(class(prior$hit_ratio_sd) != "numeric")
        stop("'prior$hit_ratio_sd' not of class 'numeric'")
    if(prior$hit_ratio_sd > 1 || prior$hit_ratio_sd < 0)
        stop("'prior$hit_ratio_sd' must be between 0 and 1")
    if(class(prior$hit_ratio_mean) != "numeric")
        stop("'prior$hit_ratio_mean' not of class 'numeric'")
    if(prior$hit_ratio_mean > 1 || prior$hit_ratio_mean < 0)
        stop("'threshold' must be between 0 and 1")
    
    alpha <- ((1 - prior$hit_ratio_mean) / prior$hit_ratio_sd ^ 2 - 1 / prior$hit_ratio_mean) * prior$hit_ratio_mean ^ 2
    beta <- alpha * (1 / prior$hit_ratio_mean - 1)
    colnames(inputMatrix) <- as.character(colnames(inputMatrix))
    cids <- colnames(inputMatrix)
    sapply(cids, function(cid){
        n <- sum(inputMatrix[,cid] == 2)
        N <- sum(inputMatrix[,cid] > 0)
        pbeta(threshold, n + alpha, N - n + beta, lower.tail=F)
    })
}

# computes a prior for the crossReactivityProbability function
# by computing the hit ratio distribution for an entire bioassayDB database
crossReactivityPrior <- function(database, minTargets=20, category=FALSE){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(class(minTargets) != "numeric")
        stop("'minTargets' not of class 'numeric'")
    
    highlyScreened <- screenedAtLeast(database, minTargets) 
    scores <- sapply(highlyScreened, targetSelectivity, database=database, category=category, scoring="fraction")
    scores <- scores[! is.na(scores)]
    list(hit_ratio_mean=mean(scores), hit_ratio_sd=sd(scores))
}