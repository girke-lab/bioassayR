\name{activityMatrix}
\alias{activityMatrix}
\title{
Output a matrix of compound bioactivity
}
\description{
Returns a complete \code{matrix} with the binary activity values for every compound, in every assay within a bioassayR database. 
}
\usage{
activityMatrix(database, maxAssayLimit = 100)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{maxAssayLimit}{
A numeric value showing the maximum number of assays that can be present in the database 
before performing this computation. This is a safety feature to prevent users from accidentally
initiating this computation across thousands of assays. As this function is computationally intensive,
it's typically used only on small databases. However, by supplying a larger value users
can increase this limit.
}
}
\value{
A complete \code{matrix} with the binary activity values for every compound, in every assay within a bioassayR database.
}
\author{
Tyler Backman
}
\examples{
## create a temp database and copy a few assays into it
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)
myDatabaseFilename <- tempfile()
mydb <- newBioassayDB(myDatabaseFilename, indexed=FALSE)
addDataSource(mydb, description="bioassayR_testdata", version="unknown")
loadBioassay(mydb, getAssay(sampleDB, 53224))
loadBioassay(mydb, getAssay(sampleDB, 53211))
loadBioassay(mydb, getAssay(sampleDB, 207758))

## compute a complete matrix from this temp database
myMatrix <- activityMatrix(mydb)
myMatrix

## disconnect from databases, and delete temp database
disconnectBioassayDB(mydb)
disconnectBioassayDB(sampleDB)
unlink(myDatabaseFilename)
}
\keyword{ utilities }
