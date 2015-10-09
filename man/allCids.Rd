\name{allCids}
\alias{allCids}
\title{
List compound cids in a \code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object
}
\description{
Returns a \code{vector} of small molecule cids contained within a
\code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object.
It can optionally only returned cids labeled as active.
}
\usage{
allCids(inputObject, activesOnly = FALSE)
}
\arguments{
  \item{inputObject}{
A \code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object to query.
}
  \item{activesOnly}{
logical. Should only active compounds be returned? Defaults to FALSE.
}
}
\value{
A \code{vector} of distinct small molecule cids. No particular order
is guranteed.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get all compound cids
myCids <- allCids(sampleDB)

## get only active compound cids
activeCids <- allCids(sampleDB, activesOnly = TRUE)

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
