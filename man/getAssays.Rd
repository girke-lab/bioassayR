\name{getAssays}
\alias{getAssays}
\title{
Retrieve multiple bioassays from a database
}
\description{
Retrieves a list of aids as a single \code{bioassaySet} matrix object
}
\usage{
getAssays(database, aids)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{aids}{
One or more assay identifier strings (aid), matching aid(s) for assays loaded into the database.
}
}
\value{
A \code{bioassaySet} object containing data from the specified assays.
}
\author{
Tyler William H Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve three assays
assays <- getAssays(sampleDB, c("673509","103","105"))
assays

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
