\name{connectBioassayDB}
\alias{connectBioassayDB}
\title{
Create a \code{BioassayDB} object connected to the specified
database file
}
\description{
This function returns a \code{BioassayDB} object for working with a pre-existing bioassayR database, already located on the users filesystem.
Users can download pre-built databases for use with this feature from \url{http://chemmine.ucr.edu/bioassayr}
}
\usage{
connectBioassayDB(databasePath, writeable = FALSE)
}
\arguments{
  \item{databasePath}{
Full path to the database file to be opened.
}
  \item{writeable}{
logical. Should the database allow data to be modified and written to?
}
}
\value{
\item{BioassayDB}{ for details see ?"BioassayDB-class"}
}
\author{
Tyler Backman
}
\examples{
## create a test database
library(bioassayR)
filename <- tempfile()
mydb <- newBioassayDB(filename, indexed=FALSE)
disconnectBioassayDB(mydb)

## connect to test database
mydb <- connectBioassayDB(filename)

## close and delete database
disconnectBioassayDB(mydb)
unlink(filename)
}
\keyword{ utilities }
