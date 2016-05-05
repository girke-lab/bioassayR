context("crossReactivityProbability")

activityMatrix <- Matrix(c(
    # target1 target2
    2,2,2,2   ,2,2, # all actives
    0,0,0,0   ,0,0, # all untested
    2,2,1,1   ,1,2, # mixed
    2,1,1,1   ,1,2, # mixed
    2,2,2,1   ,1,1, # mixed, inactive
    0,0,1,0   ,0,0, # some untested
    0,0,2,0   ,0,0  # some untested
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

targetMatrix <- perTargetMatrix(sample_bioassaySet, inactives = T,
                                summarizeReplicates = "mode")

test_that("cross-reactivites are correct", {
    # check that numbers match those published here:
    # Dančík, V. et al. Connecting Small Molecules with Similar Assay Performance 
    # Profiles Leads to New Biological Hypotheses. J Biomol Screen 19, 771–781 (2014).
    # note that the prior hit_ratio_sd mentioned in the paper is actually the variance
    expected <- c(0.57965494, 0.07445909, 0.07445909, 0.28839271, 0.35165635, 0.09860181)
    computed <- crossReactivityProbability(targetMatrix,
                                           threshold=0.25,
                                           prior=list(hit_ratio_mean=0.126, hit_ratio_sd=0.1072))
    names(expected) <- names(computed)
    expect_equal(computed, expected)
})