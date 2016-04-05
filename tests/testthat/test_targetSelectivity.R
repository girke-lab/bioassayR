context("targetSelectivity")

# load sample database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

# get two most highly screened cids
topCids <- queryBioassayDB(sampleDB, paste("SELECT cid, COUNT(DISTINCT aid) AS assays", 
                                           "FROM activity WHERE cid NOT NULL GROUP BY cid", 
                                           "ORDER BY assays DESC LIMIT 2"))$cid

# test code: try ALL cids
# topCids <- queryBioassayDB(sampleDB, paste("SELECT cid, COUNT(DISTINCT aid) AS assays", 
#                                           "FROM activity WHERE cid NOT NULL GROUP BY cid", 
#                                           "ORDER BY assays DESC"))$cid

test_that("total selectivity results are correct with all targets", {
    target <- targetSelectivity(sampleDB, topCids, scoring="total", category=FALSE, multiTarget="all")
    current <- sapply(topCids, function(cid){
        length(unique(row.names(activeTargets(sampleDB, cid))))    
    })
    names(current) <- topCids
    expect_equal(target,current)
})

test_that("fraction selectivity results are correct with all targets", {
    target <- targetSelectivity(sampleDB, topCids, scoring="fraction", category=FALSE, multiTarget="all")
    current <- sapply(topCids, function(cid){
        activeTargets <- unique(row.names(activeTargets(sampleDB, cid)))
        inactiveTargets <- unique(row.names(inactiveTargets(sampleDB, cid)))
        inactiveTargets <- inactiveTargets[! inactiveTargets %in% activeTargets]
        length(activeTargets) / (length(activeTargets) + length(inactiveTargets))
    })
    names(current) <- topCids
    expect_equal(target,current)
})

test_that("fraction selectivity results are correct with multi-target assays dropped", {
    target <- targetSelectivity(sampleDB, topCids, scoring="fraction", category=FALSE, multiTarget="drop")
    current <- sapply(topCids, function(cid){
        activityData <- queryBioassayDB(sampleDB, paste("SELECT aid, activity, target FROM activity NATURAL JOIN targets WHERE activity NOT NULL AND cid =", cid))
        duplicateAids <- activityData$aid[duplicated(activityData$aid)]
        activityData <- activityData[! activityData$aid %in% duplicateAids,,drop=F]
        activeTargets <- unique(activityData$target[activityData$activity == 1])
        inactiveTargets <- unique(activityData$target[activityData$activity == 0])
        inactiveTargets <- inactiveTargets[! inactiveTargets %in% activeTargets]
        length(activeTargets) / (length(activeTargets) + length(inactiveTargets))
    })
    names(current) <- topCids
    expect_equal(target,current)
})

test_that("fraction selectivity results are correct with multi-target assays keep one", {
    target <- targetSelectivity(sampleDB, topCids, scoring="fraction", category=FALSE, multiTarget="keepOne")
    current <- sapply(topCids, function(cid){
        activityData <- queryBioassayDB(sampleDB, paste("SELECT aid, activity, target FROM activity NATURAL JOIN targets WHERE activity NOT NULL AND cid =", cid))
        activityData <- activityData[! duplicated(activityData$aid),,drop=F]
        activeTargets <- unique(activityData$target[activityData$activity == 1])
        inactiveTargets <- unique(activityData$target[activityData$activity == 0])
        inactiveTargets <- inactiveTargets[! inactiveTargets %in% activeTargets]
        length(activeTargets) / (length(activeTargets) + length(inactiveTargets))
    })
    names(current) <- topCids
    expect_equal(target,current)
})

test_that("total kclust selectivity results are correct with all targets", {
    target <- targetSelectivity(sampleDB, topCids, scoring="total", category="kClust", multiTarget="all")
    current <- sapply(topCids, function(cid){
        allActiveTargets <- unique(row.names(activeTargets(sampleDB, cid)))
        allClusters <- lapply(allActiveTargets, translateTargetId, database=sampleDB, category="kClust")
        names(allClusters) <- allActiveTargets
        allClusters <- unlist(allClusters)
        length(unique(names(allClusters)[! duplicated(allClusters)]))
    })
    names(current) <- topCids
    expect_equal(target,current)
})

test_that("fraction kclust selectivity results are correct with all targets", {
    target <- targetSelectivity(sampleDB, topCids, scoring="fraction", category="kClust", multiTarget="all")
    current <- sapply(topCids, function(cid){
        allActiveTargets <- unique(row.names(activeTargets(sampleDB, cid)))
        allInactiveTargets <- unique(row.names(inactiveTargets(sampleDB, cid)))
        allInactiveTargets <- allInactiveTargets[! allInactiveTargets %in% allActiveTargets]
        allClusters <- lapply(c(allActiveTargets, allInactiveTargets), translateTargetId, database=sampleDB, category="kClust")
        names(allClusters) <- c(allActiveTargets, allInactiveTargets)
        allClusters <- unlist(allClusters)
        naNames <- names(allClusters)[is.na(allClusters)]
        allClusters <- allClusters[! is.na(allClusters)]
        keepNames <- c(names(allClusters)[! duplicated(allClusters)], naNames)
        allActiveTargets <- allActiveTargets[allActiveTargets %in% keepNames]
        allInactiveTargets <- allInactiveTargets[allInactiveTargets %in% keepNames]
        length(allActiveTargets) / (length(allActiveTargets) + length(allInactiveTargets))
    })
    names(current) <- topCids
    expect_equal(target,current)
})

# disconnect from sample database
disconnectBioassayDB(sampleDB)