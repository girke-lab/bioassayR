\name{assaySetTargets}
\alias{assaySetTargets}
\title{
Return targets of assays in a \code{bioassaySet} object
}
\description{
    This takes a \code{bioassaySet} of multiple assays and returns
    a vector of the targets of each, with the assay identifiers themselves (aids)
    as names. If a single assay contains multiple targets, these will all be listed.
}
\usage{
assaySetTargets(assays)
}
\arguments{
  \item{assays}{
A \code{bioassaySet} object with data from multiple assays, some of which may share a common target.
}
}
\value{
A \code{character} vector
of the targets of each, with the assay identifiers themselves (aids) as names
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

## get the targets for these assays
myTargets <- assaySetTargets(assays)
myTargets

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
