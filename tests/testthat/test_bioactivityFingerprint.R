context("bioactivityFingerprint")

# check correctness of example matrix
test_that("example matrix is correct", {
    ## connect to a test database
    extdata_dir <- system.file("extdata", package="bioassayR")
    sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
    sampleDB <- connectBioassayDB(sampleDatabasePath)
    
    ## retrieve all targets in database
    targetList <- allTargets(sampleDB)
    
    ## get an activity fingerprint object for selected CIDs
    queryCids <- c("2244", "3715", "2662", "3033", "133021", 
                   "44563999", "44564000", "44564001", "44564002") 
    myAssaySet <- getBioassaySetByCids(sampleDB, queryCids)
    myFp <- bioactivityFingerprint(bioassaySet=myAssaySet)
    
    ## convert to matrix
    fpMatrix <- as.matrix(myFp)
    
    # check that active targets in database match
    # those in activity matrix
    for(cid in queryCids){
        # note for now this uses duplicated() to remove multiple targets from
        # multiple target assays. If perTargetMatrix() is updated, this needs
        # changing
        targetVsAid <- queryBioassayDB(sampleDB, paste("SELECT aid, target FROM activity NATURAL JOIN targets WHERE cid =", cid, " AND activity = 1"))
        targetsOnePerAssay <- unique(targetVsAid$target[! duplicated(targetVsAid$aid)])
        expect_equal(table(fpMatrix[cid,] == 1)[["TRUE"]], length(targetsOnePerAssay))
    }
    
    ## disconnect from sample database
    disconnectBioassayDB(sampleDB)
})