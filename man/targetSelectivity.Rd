\name{targetSelectivity}
\alias{targetSelectivity}
\title{
Returns the target selectivity for a specified list of compounds (cids).
}
\description{
Queries a \code{BioassayDB} database and returns the target selectivity of the specified cids.
}
\usage{
targetSelectivity(database, cids, scoring = "total", 
                  category=FALSE, multiTarget="keepOne")
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
If fractional activity is requested, active values take precedence over inactives: if a target is both active and inactive in different assays it will be regarded as active.
}
  \item{category}{
Include only once in selectivity counts any targets which share a common annotation of this category
(as used by the \code{translateTargetId} and \code{loadIdMapping} functions). For example,
with the PubChem BioAssay database one could use "UniProt", "kClust", or "domains" to get
selectivity by targets with unique UniProt identifiers, distinct amino acid sequences, or Pfam domains
respectively (the latter is also known as domain selectivity).
}
  \item{multiTarget}{
Decides how selectivity is counted with regard to multi-target assays. If \code{"drop"} these
assays are excluded entirely. If \code{"keepOne"} only the first target in the database is considered.
If \code{"all"} they are counted separately towards the total.
}
}
\value{
Returns an \code{numeric} vector containing the target selectivity for each query compound. Returned entires are named by their corresponding cid. 
}
\author{
Tyler Backman
}
\seealso{
\code{\link{translateTargetId}}
\code{\link{loadIdMapping}}
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
