\name{screenedAtLeast}
\alias{screenedAtLeast}
\title{
Return all compounds in the database screened at least 'minTargets' times.
}
\description{
Returns all compound cids screened against at least 'minTargets' distinct target identifiers. For a very large database (such as PubChem Bioassay) this function may take a long time to run.
}
\usage{
screenedAtLeast(database, minTargets)
}
\arguments{
  \item{database}{
    A \code{BioassayDB} database to query.
}
  \item{minTargets}{
    The minimum number of distinct targets for each returned cid.
}
}
\value{
Returns a \code{character} vector of all CIDs meeting the specified criteria.
}
\author{
Tyler Backman
}

\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get all CIDS screened against at least 2 distinct targets
highlyScreened <- screenedAtLeast(sampleDB, 2)
highlyScreened

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
