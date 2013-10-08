setClass("BioAssaySet", representation=representation(database="SQLiteConnection"))

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
