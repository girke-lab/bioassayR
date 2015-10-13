\name{selectiveAgainst}
\alias{selectiveAgainst}
\title{
Identify small molecules with selective binding against a target of interest
}
\description{
Allows the user to find compounds in the database that have been screened
against a large number of distinct targets, but show high binding selectivity
for a specific target of interest.
}
\usage{
selectiveAgainst(database, target, maxCompounds = 10, minimumTargets = 10)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{target}{
A string or integer containing a target_id referring to a target of interest.
}
  \item{maxCompounds}{
An integer representing the number of resulting compounds to return.
}
  \item{minimumTargets}{
An integer representing the minimum number of distinct targets a compound
must have been screened against to be included in the results.
}
}
\value{
A \code{data.frame} where the row names represent each compound showing binding specificity against the specified target.
The first column shows the number of distinct targets each compound shows activity against, and the second
column shows the total number of distinct targets it has been screened against.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## find target selective compounds active against a protein of interest
selectiveAgainst(sampleDB, target="166897622", maxCompounds=10,minimumTargets=20)

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
