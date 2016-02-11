setClass("BioassayDB", representation=representation(database="SQLiteConnection"))

setClass("bioassay", representation=representation(
    aid = "character",
    source_id = "character",
    assay_type = "character",
    organism = "character",
    scoring = "character",
    targets = "character",
    target_types = "character",
    scores = "data.frame"
))

setClass("bioassaySet", representation=representation(
    activity = "dgCMatrix",
    scores = "dgCMatrix",
    targets = "dgCMatrix",
    sources = "data.frame",
    source_id = "integer",
    assay_type = "character",
    organism = "character",
    scoring = "character",
    target_types = "character"
))
