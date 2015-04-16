\name{queryBioassayDB}
\alias{queryBioassayDB}
\title{
Perform a SQL query on a bioassayR database
}
\description{
Provides extreme query flexibility by allowing the user to perform any SQLite query on a bioassayR database. This allows for analysis beyond that provided by the built in query functions.
}
\usage{
queryBioassayDB(object, query)
}
\arguments{
  \item{object}{
A \code{BioassayDB} object referring to a bioassayR database.
}
  \item{query}{
A string containing a valid SQLite query (see SQLite documentation for more details).
}
}
\value{
A \code{data.frame} containing the results of the specified query.
}
\references{
http://www.sqlite.org provides a complete reference for SQLite syntax that can be used with this function
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## inspect the structure of the database before forming a query
queryBioassayDB(sampleDB, "SELECT * FROM sqlite_master WHERE type='table'")

## find all activity data for compound cid 2244 
queryBioassayDB(sampleDB, "SELECT * FROM activity WHERE cid = '2244'")

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
