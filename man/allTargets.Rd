\name{allTargets}
\alias{allTargets}
\title{
List distinct target(s) in a \code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object
}
\description{
Returns a \code{vector} of target ids contained within a
\code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object.
}
\usage{
allTargets(inputObject)
}
\arguments{
  \item{inputObject}{
A \code{BioassayDB}, \code{bioassay}, \code{bioassaySet}, or target matrix (\code{dgCMatrix}) object to query.
}
}
\value{
A \code{vector} of distinct target ids. No particular order
is guaranteed.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get all target ids
myTargets <- allTargets(sampleDB)

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
