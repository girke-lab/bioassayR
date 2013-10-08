\name{dropBioassayIndex}
\alias{dropBioassayIndex}
\title{
Remove index from a bioassayR database
}
\description{
Indexing a bioassayR database before performing queries will drastically improve query performance. However, it will also slow down loading large amounts of additional data. Therefore, 
it may be necessary to use this index to remove an index from a database before adding large quantities of data. Afterwards, the index can be re-generated using the \code{addBioassayIndex} function. 
}
\usage{
dropBioassayIndex(database)
}
\arguments{
  \item{database}{
A \code{BioAssaySet} database to have the index removed.
}
}
\author{
Tyler Backman
}
\examples{
## create test database
library(bioassayR)
filename <- tempfile()
mydb <- newBioassayDB(filename, indexed=TRUE)

## remove database index
dropBioassayIndex(mydb)

## load new data into database here 

## reactivate index
addBioassayIndex(mydb)

## close and delete test database
disconnectBioassayDB(mydb)
unlink(filename)
}
\keyword{ utilities }
