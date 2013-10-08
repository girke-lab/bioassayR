\name{disconnectBioassayDB}
\alias{disconnectBioassayDB}
\title{
Disconnect the database file from a \code{BioAssaySet} object
}
\description{
This function disconnects the underlying sqlite database from a \code{BioAssaySet} object. This is a critical step for writeable databases, but can be omitted for read only databases.
}
\usage{
disconnectBioassayDB(database)
}
\arguments{
  \item{database}{
A code{BioAssaySet} object to be disconnected.
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

## disconnect from database
mydb <- connectBioassayDB(filename)

## delete database file
unlink(filename)
}
\keyword{ utilities }
