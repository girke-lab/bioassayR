\name{getBioassaySetByCids}
\alias{getBioassaySetByCids}
\title{
Create \code{bioassaySet} sparse matrix object with activity data only for specified compounds
}
\description{
Takes a list of compounds, and creates a \code{bioassaySet} sparse matrix object with the activity data
for these compounds only, not including activity data from other compounds in the same assays.
}
\usage{
getBioassaySetByCids(database, cids)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{cids}{
One or more compounds IDs of interest.
}
}
\value{
A \code{bioassaySet} object containing data from the specified cids.
}
\author{
Tyler William H Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve activity data on 3 compounds 
activitySet <- getBioassaySetByCids(sampleDB, c("2244","3715","237"))
activitySet

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
