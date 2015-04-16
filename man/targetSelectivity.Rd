\name{targetSelectivity}
\alias{targetSelectivity}
\title{
Returns the target selectivity for a specified list of compounds (cids).
}
\description{
Queries a \code{BioassayDB} database and returns the target selectivity of the specified cids.
}
\usage{
targetSelectivity(database, cids, scoring = "total")
}
\arguments{
  \item{database}{
    A \code{BioassayDB} database to query.
}
  \item{cids}{
A string or integer vector containing query cids referring to a small molecules.
}
  \item{scoring}{
Must be one of two optional scoring methods "total" or "fraction". Fraction returns the target
selectivity for each compound as the fraction of screened distinct targets that showed activity in at least one assay. Total returns the total number of active distinct targets for each compound, and does not consider inactive targets in the calculation.
}
}
\value{
Returns an \code{numeric} vector containing the target selectivity for each query compound. Returned entires are named by their corresponding cid. 
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## make a vector with compounds of interest
compoundsOfInterest <- c(2244, 2662, 3033)

## get "total" active targets for each compound of interest
targetSelectivity(sampleDB, compoundsOfInterest, scoring="total")

## get fraction of active targets for each compound of interest
targetSelectivity(sampleDB, compoundsOfInterest, scoring="fraction")

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
