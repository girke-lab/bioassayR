\name{loadBioassay}
\alias{loadBioassay}
\title{
Add an assay to the database
}
\description{
Loads the results of a bioassay experiment (stored as a \code{bioassay} object) into the specified database.
The data source specified in the \code{bioassay} object be added to the database with
\code{addDataSource} before loading. If the assay identifier (aid) is already present in the database,
an error is returned and no additional data is loaded.
}
\usage{
loadBioassay(database, bioassay)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to load the data into.
}
  \item{bioassay}{
A \code{bioassay} object containing the data to load.
}
}
\author{
Tyler Backman
}
\examples{
## create sample database
myDatabaseFilename <- tempfile()
mydb <- newBioassayDB(myDatabaseFilename, indexed=FALSE)

## parse example assay data
extdata_dir <- system.file("extdata", package="bioassayR")
assayDescriptionFile <- file.path(extdata_dir, "exampleAssay.xml")
activityScoresFile <- file.path(extdata_dir, "exampleScores.csv")
myAssay <- parsePubChemBioassay("1000", activityScoresFile, assayDescriptionFile)

## load bioassay into database
addDataSource(mydb, description="PubChem BioAssay", version="unknown")
loadBioassay(mydb, myAssay)

## disconnect from and delete sample database
disconnectBioassayDB(mydb)
unlink(myDatabaseFilename)
}
\keyword{ utilities }
