\name{perTargetMatrix}
\alias{perTargetMatrix}
\title{
Collapse a \code{bioassaySet} object from multiple assays by combining assays with a common target
}
\description{
Creates a binary \code{sparseMatrix} object which has an activity value for each distinct target identifier rather than each distinct assay. If any assay for a given target vs compound combination shows active, this combination will be marked active in the resulting object.
}
\usage{
perTargetMatrix(assays)
}
\arguments{
  \item{assays}{
A \code{bioassaySet} object with data from multiple assays, some of which may share a common target.
}
}
\value{
A \code{sparseMatrix} which contains a value of 1 for each target vs compound combination which shows activity in at least one parent assay.
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

## collapse assays into perTargetMatrix
targetMatrix <- perTargetMatrix(assays)
targetMatrix

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
