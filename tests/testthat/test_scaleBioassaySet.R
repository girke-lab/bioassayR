context("scaleBioassaySet")

# create sample sparse activity matrix
assays <- c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,2,5,2,4,3)
cids <- c(0,0,0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,6,6,5)
scores <- c(2,2,2,2,2,2,2,2,1,1,1,2,2,1,1,1,1,-2.5,2,2,2,1,1,1,1,2,2,1,0)
dimnames <- list(c(10,20,30,40,50,60),c(1,2,3,4,5,6,7))

activityMatrix <- sparseMatrix(
    i = assays,
    j = cids,
    x = scores,
    dims = sapply(dimnames, length),
    dimnames = dimnames,
    symmetric = FALSE,
    index1 = FALSE)

sample_bioassaySet <- new("bioassaySet",
                          activity = activityMatrix,
                          scores = activityMatrix,
                          targets = as(Matrix(sparse=T), "dgCMatrix"),
                          sources = data.frame(),
                          source_id = integer(),
                          assay_type = character(),
                          organism = character(),
                          scoring = character(),
                          target_types = character())

standard_scores <- scaleBioassaySet(sample_bioassaySet)@scores

test_that("z-score matrix is correct", {
    expected_scores <- sparseMatrix(
        i = c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,2,3,5,2,4),
        j = c(0,0,0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,6,6),
        x = c(NaN,0.5,0.912870929175277,1.41421356237309,1.78885438199983,0.564288093646835,
              NaN,0.5,-0.912870929175277,0,-0.447213595499958,0.564288093646835,
              NaN,-1.5,-0.912870929175277,0,-0.447213595499958,-1.74416319854476,
              NaN,0.5,0.912870929175277,0,-0.447213595499958,0.0512989176042577,
              -0.912870929175277,-1.41421356237309,0.564288093646835,
              0.912870929175277,-0.447213595499958), # each row here is a cid
        dims = sapply(dimnames, length),
        dimnames = dimnames,
        symmetric = FALSE,
        index1 = FALSE)
    expect_equal(standard_scores, expected_scores)
})
