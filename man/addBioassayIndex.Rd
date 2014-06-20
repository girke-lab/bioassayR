\name{addBioassayIndex}
\alias{addBioassayIndex}
\title{
Index a bioassayR database
}
\description{
Indexing a bioassayR database before performing queries will drastically improve query performance. However, it will also slow down loading large amounts of additional data. Therefore, we recommend loading the majority of your data, using this function to index, and then performing queries.
}
\usage{
addBioassayIndex(database)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to be indexed.
}
}
\author{
Tyler Backman
}
\examples{
## create test database
library(bioassayR)
filename <- tempfile()
mydb <- newBioassayDB(filename, indexed=FALSE)

## load any data at this point

## add database index
addBioassayIndex(mydb)

# perform queries here

## close and delete test database
disconnectBioassayDB(mydb)
unlink(filename)
}
\keyword{ utilities }
