\name{screenedAtLeast}
\alias{screenedAtLeast}
\title{
Return all compounds in the database screened at least 'minTargets' times.
}
\description{
Returns all compound cids screened against at least 'minTargets' distinct target identifiers. For a very large database (such as PubChem Bioassay) this function may take a long time to run.
}
\usage{
screenedAtLeast(database, minTargets, inconclusives=TRUE)
}
\arguments{
  \item{database}{
    A \code{BioassayDB} database to query.
}
  \item{minTargets}{
    The minimum number of distinct targets for each returned cid.
}
  \item{inconclusives}{
    Logical. If \code{TRUE} (the default) screening results with inconclusive (\code{NA})
    are counted towards \code{minTargets}. If \code{FALSE} only active and inactive
    results are counted.
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

## get all CIDS screened against at least 2 distinct targets with conclusive results
highlyScreened <- screenedAtLeast(sampleDB, 2, inconclusives=FALSE)
highlyScreened

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
