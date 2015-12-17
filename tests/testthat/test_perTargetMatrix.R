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

test_that("mode matrix is correct", {
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
                                    conflictResolver = "mode")
    expect_equal(targetMatrix, expectedTargetMatrix)
})

test_that("activesFirst matrix is correct", {
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
                                    conflictResolver = "activesFirst")
    expect_equal(targetMatrix, expectedTargetMatrix)
})
