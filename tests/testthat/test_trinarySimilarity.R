context("trinarySimilarity")

# create sample sparse matrix
# cids:      1              2              3              4              5              6              7
targets <- c(0,1,2,3,4,5,   0,1,2,3,4,5,   0,1,2,3,4,5,   0,1,2,3,4,5,   0,1,2,3,4,5,   0,1,2,3,4,5,   0,1,2,3,4,5)
cids <-    c(0,0,0,0,0,0,   1,1,1,1,1,1,   2,2,2,2,2,2,   3,3,3,3,3,3,   4,4,4,4,4,4,   5,5,5,5,5,5,   6,6,6,6,6,6)
scores <-  c(2,2,1,2,1,0,   0,0,0,0,0,0,   0,2,0,2,0,2,   0,1,0,0,2,2,   1,1,0,0,2,2,   1,1,0,2,2,1,   2,2,1,1,2,2)  

dimnames <- list(c(10,20,30,40,50,60),c(1,2,3,4,5,6,7))

testMatrix <- sparseMatrix(
    i = targets,
    j = cids,
    x = scores,
    dims = sapply(dimnames, length),
    dimnames = dimnames,
    symmetric = FALSE,
    index1 = FALSE)

test_that("similarity search result is correct", {
    # 1 = query, 2 = too few shared actives, 3 = enough shared actives, 
    # 4 = too few shared targets, 5 = enough shared targets
    # 6 = 0.25, 7 = 0.50
    expected_result <- c(1, NA, 1, NA, 0, 0.25, 0.5)
    
    computed_result <- trinarySimilarity(testMatrix[,1,drop=F], testMatrix, minSharedScreenedTargets = 3, minSharedActiveTargets = 2)
    names(expected_result) <- names(computed_result)
    
    expect_equal(computed_result, expected_result)
})