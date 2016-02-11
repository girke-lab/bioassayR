# connects an existing bioassay database
connectBioassayDB <- function (databasePath, writeable = F) {
    if(! file.exists(databasePath)){
        stop("input database does not exist")
    }

    drv <- dbDriver("SQLite")
    if(writeable){
        con <- dbConnect(drv, dbname=databasePath, flags=SQLITE_RW)
    } else {
        con <- dbConnect(drv, dbname=databasePath, flags=SQLITE_RO)
    }
    new("BioassayDB", database = con)
}

# disconnects a bioassay database
disconnectBioassayDB <- function (database) {
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")

    dbDisconnect(slot(database, "database"))
    invisible()
}

# creates a new bioassay database
newBioassayDB <- function(databasePath, writeable = T, indexed = F){
    if(file.exists(databasePath)){
        stop("database filename already exists")
    }   
    if(! is.logical(writeable))
        stop("input must be of type 'logical'")
    if(! is.logical(indexed))
        stop("input must be of type 'logical'")

    drv <- dbDriver("SQLite")
    con <- dbConnect(drv, dbname=databasePath)
    dbGetQuery(con, paste("CREATE TABLE activity",
        "(aid INTEGER, cid INTEGER,",
        "activity INTEGER, score INTEGER)"))
    dbGetQuery(con, paste("CREATE TABLE assays",
        "(source_id INTEGER, aid INTEGER,",
        "assay_type TEXT, organism TEXT, scoring TEXT)"))
    dbGetQuery(con, paste("CREATE TABLE sources",
        "(source_id INTEGER PRIMARY KEY ASC,",
        "description TEXT, version TEXT)"))
    dbGetQuery(con, paste("CREATE TABLE targets",
        "(aid INTEGER, target TEXT, target_type TEXT)"))
    dbGetQuery(con, paste("CREATE TABLE targetTranslations",
        "(target TEXT, category TEXT, identifier TEXT)"))
    dbDisconnect(con)
    database <- connectBioassayDB(databasePath, writeable = T)
    if(indexed){
        addBioassayIndex(database)
    }
    if(! writeable){
        disconnectBioassayDB(database)
        database <- connectBioassayDB(databasePath, writeable = F)
    }
    return(database)
}

# adds a new data source
addDataSource <- function(database, description, version){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    if(is.numeric(version)){
        version <- as.character(version)
    } else if(! is.character(version)){
        stop("input not class character")
    }
    if(length(version) > 1){
        warning("too many inputs, only the first was kept")
        version <- version[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s\\.]+$", version, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")
    if(is.numeric(description)){
        description <- as.character(description)
    } else if(! is.character(description)){
        stop("input not class character")
    }
    if(length(description) > 1){
        warning("too many inputs, only the first was kept")
        description <- description[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s\\.]+$", description, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")

    con <- slot(database, "database")
    sql <- "INSERT INTO sources VALUES (NULL, $DESCRIPTION, $VERSION)"
    dbBegin(con)
    dbGetPreparedQuery(con, sql, bind.data = data.frame(DESCRIPTION=description, VERSION=version))
    dbCommit(con)
    invisible()
}

# parses input files from PubChem Bioassay
parsePubChemBioassay <- function(aid, csvFile, xmlFile, duplicates = "drop"){
    if(! file.exists(csvFile)){
        stop("csv file doesn't exist")
    }   
    if(! file.exists(xmlFile)){
        stop("xml file doesn't exist")
    }   
    if(is.numeric(aid)){
        aid <- as.character(aid)
    } else if(! is.character(aid)){
        stop("input not class character")
    }
    if(length(aid) > 1){
        warning("too many inputs, only the first was kept")
        aid <- aid[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", aid, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")
    if((duplicates != "drop") && (! is.logical(duplicates)))
        stop("duplicates option is not one of the allowed values")

    aid <- as.character(aid)

    # parse csv
    # note: custom CSV parser was written here to accomodate some files from PubChem that have
    #   hard to parse unescaped commas in CSV comment lines
    csvLines <- readLines(csvFile)
    csvLines <- csvLines[! grepl("^RESULT_", csvLines)]
    csvLines <- csvLines[! grepl("^\\s*$", csvLines)]
    if(length(csvLines) < 2){
        tempAssay <- t(data.frame(row.names=c("cid", "activity", "score")))
        tempAssay <- as.data.frame(tempAssay)
    } else {
        header <- strsplit(csvLines[[1]], "\\s*,\\s*")[[1]]
        if(sum(c("PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME", "PUBCHEM_ACTIVITY_SCORE") %in% header) < 3)
            stop("invalid header line")
        csvLines <- csvLines[2:length(csvLines)]
        tempAssay <- do.call(rbind, 
            lapply(csvLines, function(line) {
                newFields <- rep("", length(header))
                newData <- strsplit(line, "\\s*,\\s*")[1:length(header)][[1]]
                if(length(newData) > length(newFields)){
                    newFields <- newData[1:length(newFields)]
                } else {
                    newFields[1:length(newData)] <- newData
                }
                return(newFields)
            })
        )
        tempAssay <- as.data.frame(tempAssay, stringsAsFactors=FALSE)
        colnames(tempAssay) <- header
        tempAssay <- tempAssay[,c("PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME", "PUBCHEM_ACTIVITY_SCORE")]
        outcomes <- rep(NA, nrow(tempAssay))
        outcomes[tempAssay[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Active"] <- 1
        outcomes[tempAssay[,"PUBCHEM_ACTIVITY_OUTCOME"] == 1] <- 0
        outcomes[tempAssay[,"PUBCHEM_ACTIVITY_OUTCOME"] == "Inactive"] <- 0
        outcomes[tempAssay[,"PUBCHEM_ACTIVITY_OUTCOME"] == 2] <- 1
        tempAssay[,"PUBCHEM_ACTIVITY_OUTCOME"] <- outcomes
        colnames(tempAssay) <- c("cid", "activity", "score")
        if(sum(as.integer(tempAssay$cid) == as.character(tempAssay$cid)) < length(tempAssay$cid))
            stop("non-integer cid")
        if(sum(as.integer(tempAssay$score) == as.character(tempAssay$score)) < length(tempAssay$score))
            stop("non-integer score")
    }

    # handle duplicate cids
    if(length(unique(tempAssay$cid)) != length(tempAssay$cid)){
        if(duplicates == "drop"){
            warning("dropping duplicate cids")
            tempAssay <- tempAssay[! duplicated(tempAssay$cid),,drop=FALSE]
        } else if(! duplicates)
            stop("duplicate cids present")
    }

    # parse xmlFile
    xmlLines <- readLines(xmlFile)
    xmlLines <- paste(xmlLines, collapse="\n")
    xmlPointer <- xmlTreeParse(xmlLines, useInternalNodes=TRUE, addFinalizer=TRUE)
    targets <- xpathSApply(xmlPointer, "//x:PC-AssayTargetInfo_mol-id/text()", xmlValue, namespaces="x")
    targetTypes <- xpathSApply(xmlPointer,"//x:PC-AssayTargetInfo_molecule-type/@value", namespaces="x")
    type <- xpathSApply(xmlPointer, "//x:PC-AssayDescription_activity-outcome-method/@value", namespaces="x")[[1]]
    comments <- xpathSApply(xmlPointer, "//x:PC-AssayDescription_comment_E/text()", xmlValue, namespaces="x")
    scoring <- xpathSApply(xmlPointer, "//x:PC-ResultType_name/text()", xmlValue, namespaces="x")[[1]]
    free(xmlPointer)
    organism <- comments[grep("^Organism:\\W", comments)]
    organism <- gsub("^Organism:\\W(.*)$", "\\1", organism)[1]
    if(is.null(type)){
        type <- NA
    }
    if(is.null(targets)){
        targets <- NA
    }
    if(is.null(targetTypes)){
        targetTypes <- NA
    }
    if(is.null(scoring)){
        scoring <- NA
    }
    if(length(targets) != length(targetTypes)){
        stop(paste("error with aid", aid))
    }
    targetTypes <- targetTypes[! duplicated(targets)]
    targets <- targets[! duplicated(targets)]

    new("bioassay",
        aid=aid,
        source_id="PubChem BioAssay",
        assay_type=as.character(type),
        organism=as.character(organism),
        scoring=as.character(scoring),
        targets = as.character(targets),
        target_types = as.character(targetTypes),
        scores=tempAssay
    ) 
}

# delete an assay
dropBioassay <- function(database, aid){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    if(is.numeric(aid)){
        aid <- as.character(aid)
    } else if(! is.character(aid)){
        stop("input not class character")
    }
    if(length(aid) > 1){
        warning("too many inputs, only the first was kept")
        aid <- aid[1]
    }
    if(! grepl("^[a-zA-Z_0-9\\s]+$", aid, perl=T))
        stop("invalid input: must contain only alphanumerics and/or whitespace")

    con <- slot(database, "database")
    for(x in c("activity", "assays", "targets")){
        sql <- paste("DELETE FROM", x, "WHERE aid = $AID")
        dbBegin(con)
        dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aid))
        dbCommit(con)
    }
    invisible()
}

# adds a new assay
loadBioassay <- function(database, bioassay){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(class(bioassay) != "bioassay")
        stop("input not of class bioassay")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 

    con <- slot(database, "database")

    # check that aid is not already present 
    existingAid <- queryBioassayDB(database,
        paste("SELECT aid FROM assays",
        " WHERE aid = '", aid(bioassay), "' LIMIT 1", sep=""))$aid
    if(length(existingAid) != 0)
        stop("an assay with this aid is already present in the database")

    # get source id
    source_int <- queryBioassayDB(database,
        paste("SELECT source_id FROM sources",
        " WHERE description = '", source_id(bioassay), "' LIMIT 1", sep=""))$source_id
    if(length(source_int) < 1){
        stop("data source not found in database - please load it first with addDataSource")
    }

    # load assay details
    sql <- paste("INSERT INTO assays VALUES ('", source_int, "', $AID, $ASSAY_TYPE, $ORGANISM, $SCORING)", sep="")
    dbBegin(con)
    dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aid(bioassay), ASSAY_TYPE=assay_type(bioassay), ORGANISM=organism(bioassay), SCORING=scoring(bioassay)))
    dbCommit(con)

    # load targets
    if(length(targets(bioassay)) > 0){
        for(i in 1:length(targets(bioassay))){
            .addAssayTarget(database, aid(bioassay), targets(bioassay)[[i]], target_types(bioassay)[[i]])
        }    
    }

    # load activity scores
    .loadScores(database, aid(bioassay), scores(bioassay))

    invisible()
}

# adds a new target
.addAssayTarget <- function(database, aid, target, target_type){
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    con <- slot(database, "database")
    sql <- paste("INSERT INTO targets VALUES ($AID, $TARGET, $TARGET_TYPE)", sep="")
    dbBegin(con)
    dbGetPreparedQuery(con, sql, bind.data = data.frame(AID=aid, TARGET=target, TARGET_TYPE=target_type))
    dbCommit(con)
    invisible()
}

# loads activity scores
.loadScores <- function(database, aid, dataTable){
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    dataTable <- as.data.frame(dataTable)
    colnames(dataTable) <- c("CID", "ACTIVITY", "SCORE")
    con <- slot(database, "database")
    sql <- paste("INSERT INTO activity VALUES ('", aid, "', $CID, $ACTIVITY, $SCORE)", sep="")
    dbBegin(con)
    dbGetPreparedQuery(con, sql, bind.data = dataTable)
    dbCommit(con)
    invisible()
}

# creates database indicies 
addBioassayIndex <- function(database){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    message("Creating index: note this may take a long time for a large database")
    queryBioassayDB(database, "CREATE INDEX IF NOT EXISTS activity_aid ON activity (aid)")
    queryBioassayDB(database, "CREATE INDEX IF NOT EXISTS activity_cid ON activity (cid)")
    queryBioassayDB(database, "CREATE INDEX IF NOT EXISTS targets_aid ON targets (aid)")
    queryBioassayDB(database, "CREATE INDEX IF NOT EXISTS targets_target ON targets (target)")
    queryBioassayDB(database, "CREATE INDEX IF NOT EXISTS targetTranslations_target ON targetTranslations (target)")
    invisible()
}

# drops database indicies
dropBioassayIndex <- function(database){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    message("Removing index: note this may take a long time for a large database")
    queryBioassayDB(database, "DROP INDEX IF EXISTS activity_aid")
    queryBioassayDB(database, "DROP INDEX IF EXISTS activity_cid")
    queryBioassayDB(database, "DROP INDEX IF EXISTS targets_aid")
    queryBioassayDB(database, "DROP INDEX IF EXISTS targets_target")
    queryBioassayDB(database, "DROP INDEX IF EXISTS targetTranslations_target")
    invisible()
}

# add a new target id mapping 
loadIdMapping <- function(database, target, category, identifier){
    if(class(database) != "BioassayDB")
        stop("database not of class BioassayDB")
    if(! .writeable(database)){
        stop("database opened in read only mode")
    } 
    sql <- paste("INSERT INTO targetTranslations VALUES ('", target, "', '", category, "', '", identifier, "')", sep="")
    queryBioassayDB(database, sql)
    invisible()
}
