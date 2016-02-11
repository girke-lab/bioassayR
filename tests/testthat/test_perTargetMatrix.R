context("perTargetMatrix")

activityMatrix <- Matrix(c(
    # target1 target2
    2,2,2,2   ,2,2, # all actives
    0,0,0,0   ,0,0, # all untested
    2,2,1,1   ,1,2, # mixed
    2,1,1,1   ,1,2, # mixed
    2,2,2,1   ,1,1, # mixed, inactive
    0,0,1,0   ,0,2, # some untested
    0,0,2,0   ,1,0  # some untested
    ), ncol=7, sparse = TRUE,
    dimnames=list(c(10,20,30,40,50,60),c(1,2,3,4,5,6,7)))

targets <- Matrix(c(
    # target1 target2
    1,1,1,1,  0,0, 
    0,0,0,0,  1,1
    ), ncol=2, sparse = TRUE,
    dimnames=list(c(10,20,30,40,50,60),c(100,200)))

sample_bioassaySet <- new("bioassaySet",
    activity = activityMatrix,
    scores = activityMatrix,
    targets = targets,
    sources = data.frame(),
    source_id = integer(),
    assay_type = character(),
    organism = character(),
    scoring = character(),
    target_types = character())

test_that("discrete mode matrix is correct", {
    expectedTargetMatrix <- Matrix(c(
        2,2,
        1,1,
        1,1,
        2,1,
        1,2,
        2,1
        ), ncol=6, sparse=T,
        dimnames=list(c(100,200), c(1,3,4,5,6,7)))

    targetMatrix <- perTargetMatrix(sample_bioassaySet, inactives = T,
                                    summarizeReplicates = "mode")
    expect_equal(targetMatrix, expectedTargetMatrix)
})

test_that("discrete activesFirst matrix is correct", {
    expectedTargetMatrix <- Matrix(c(
        2,2,
        2,2,
        2,2,
        2,1,
        1,2,
        2,1
        ), ncol=6, sparse=T,
        dimnames=list(c(100,200), c(1,3,4,5,6,7)))

    targetMatrix <- perTargetMatrix(sample_bioassaySet, inactives = T,
                                    summarizeReplicates = "activesFirst")
    expect_equal(targetMatrix, expectedTargetMatrix)
})

scaledSample_bioassaySet <- scaleBioassaySet(sample_bioassaySet)

test_that("numeric Z-score activesFirst matrix is correct", {
    expectedTargetMatrix <- Matrix(c(
        1.78885438199983,1.5,
        -1.78885438199983,0.912870929175277,
        -0.447213595499958,-1.5,
        -0.447213595499958,-0.912870929175277,
        0.447213595499958,-0.912870929175277,
        -0.447213595499958,0.912870929175277
    ), ncol=6, sparse=T,
    dimnames=list(c(200,100), c(1,5,4,3,6,7)))
    
    targetMatrix <- perTargetMatrix(scaledSample_bioassaySet, useNumericScores = T,
                                    summarizeReplicates = "activesFirst")
    expect_equal(as.matrix(targetMatrix), as.matrix(expectedTargetMatrix))
})

test_that("numeric Z-score mean matrix is correct", {
    expectedTargetMatrix <- Matrix(c(
        1.11803398874989,0.728217732293819,
        -1.11803398874989,0.228217732293819,
        0,-0.728217732293819,
        0,-0.228217732293819,
        0.447213595499958,-0.912870929175277,
        -0.447213595499958,0.912870929175277
    ), ncol=6, sparse=T,
    dimnames=list(c(200,100), c(1,5,4,3,6,7)))
    
    targetMatrix <- perTargetMatrix(scaledSample_bioassaySet, useNumericScores = T,
                                    summarizeReplicates = mean)
    expect_equal(as.matrix(targetMatrix), as.matrix(expectedTargetMatrix))
})
