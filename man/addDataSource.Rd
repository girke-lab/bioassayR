\name{addDataSource}
\alias{addDataSource}
\title{
Add a new data source to a bioassayR database
}
\description{
This function adds a new data source (name/description and version) for tracking data within a bioassayR database. This can be used later to identify the source of any specific activity data within the database, or to limit analysis to data from specific source(s).
}
\usage{
addDataSource(database, description, version)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to add a new data source to.
}
  \item{description}{
A string containing a name or description of the new data source. This exact value will be used as a key for querying and loading data from this source.
}
  \item{version}{
A string with the version and/or date of the data source. This can be used to track the date in which a non-version data source was mirrored.
}
}
\author{
Tyler Backman
}
\examples{
## create a test database
library(bioassayR)
filename <- tempfile()
mydb <- newBioassayDB(filename, indexed=FALSE)

## add a new data source
addDataSource(mydb, description="bioassayR_sample", version="1.0")

## list data sources loaded
mydb

## close and delete database
disconnectBioassayDB(mydb)
unlink(filename)
}
\keyword{ utilities }
