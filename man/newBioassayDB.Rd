\name{newBioassayDB}
\alias{newBioassayDB}
\title{
Create a new bioassayR database
}
\description{
This function creates a new bioassayR database at the specified filesystem location, and returns a \code{BioassayDB} object connected to the new database.
}
\usage{
newBioassayDB(databasePath, writeable = TRUE, indexed = FALSE)
}
\arguments{
  \item{databasePath}{
Full path to the database file to be created.
}
  \item{writeable}{
logical. Should the database allow data to be modified and written to?
}
  \item{indexed}{
logical. Should a performance enhancing index be created? The default is false, as typically an index is added only after initial data is loaded. Data loading is much slower into an already indexed database.
}
}
\author{
Tyler Backman
}
\examples{
## get a temporary filename
library(bioassayR)
filename <- tempfile()

## create a new bioassayR database
mydb <- newBioassayDB(filename, indexed=FALSE)

## close and delete database
disconnectBioassayDB(mydb)
unlink(filename)
}
\keyword{ utilities }
