\name{dropBioassay}
\alias{dropBioassay}
\title{
Delete an assay from a bioassayR database
}
\description{
Allows the user to delete all records from the database associated with a given assay identifier.
}
\usage{
dropBioassay(database, aid)
}
\arguments{
  \item{database}{
    A \code{BioassayDB} database to remove an assay from.
}
  \item{aid}{
    The assay identifier string (aid), matching an aid for an assay loaded into the database.
}
}
\author{
Tyler Backman
}
\examples{
## create sample database and load with data
myDatabaseFilename <- tempfile()
mydb <- newBioassayDB(myDatabaseFilename, indexed=FALSE)
extdata_dir <- system.file("extdata", package="bioassayR")
assayDescriptionFile <- file.path(extdata_dir, "exampleAssay.xml")
activityScoresFile <- file.path(extdata_dir, "exampleScores.csv")
myAssay <- parsePubChemBioassay("1000", activityScoresFile, assayDescriptionFile)
addDataSource(mydb, description="PubChem BioAssay", version="unknown")
loadBioassay(mydb, myAssay)

## delete the loaded assay
dropBioassay(mydb, "1000")

## disconnect from and delete sample database
disconnectBioassayDB(mydb)
unlink(myDatabaseFilename)
}
\keyword{ utilities }
