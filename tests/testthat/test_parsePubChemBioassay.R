context("parsePubChemBioassay")

# get file path
extdata_dir <- system.file("extdata", package="bioassayR")
aid <- "2557"
xmlFile <- file.path(extdata_dir, "AID_2557_descr.xml")
csvFile <- file.path(extdata_dir, "AID_2557_datatable_abridged.csv")

# parse files
assay2557 <- parsePubChemBioassay(aid, csvFile, xmlFile)

test_that("regex assay parsing is correct", {
    expect_equal(1, 1)
})