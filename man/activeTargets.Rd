\name{activeTargets}
\alias{activeTargets}
\title{
Show targets against which a small molecule is active
}
\description{
Returns a \code{data.frame} of the protein targets (target_type for the protein must be 'protein'), which a given small molecule (specified by cid) shows activity against. For each target, a single row shows the total number of distinct screens it participated in, and the fraction of those in which it exhibits activity.
}
\usage{
activeTargets(database, cid)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{cid}{
A string or integer containing a cid referring to a small molecule.
}
}
\value{
A \code{data.frame} where the row names represent each target the specified compound shows activity against, and the columns show the total screens and the fraction in which the compound was active.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get targets that compound 2244 shows activity against
myTargets <- row.names(activeTargets(sampleDB, "2244"))

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
