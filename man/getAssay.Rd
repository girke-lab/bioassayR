\name{getAssay}
\alias{getAssay}
\title{
Retrieve a bioassay
}
\description{
Retrieves a bioassay as a \code{bioassay} object from a bioassayR database by identifier.
}
\usage{
getAssay(database, aid)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{aid}{
The assay identifier string (aid), matching an aid for an assay loaded into the database.
}
}
\value{
A \code{bioassay} object containing the requested assay.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve an assay
assay <- getAssay(sampleDB, "673509")
assay

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
