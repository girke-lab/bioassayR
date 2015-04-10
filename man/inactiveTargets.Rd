\name{inactiveTargets}
\alias{inactiveTargets}
\title{
Takes a single cid and returns a table of the proteins it has been found inactive against.
}
\description{
Returns a \code{data.frame} of all targets a single cid (compound) has been found inactive against, and the number of times it has been found inactive in distinct assay experiments. If a compound has been found both active and inactive in different assays, it will be listed among these results.
}
\usage{
inactiveTargets(database, cid)
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{cid}{
A string or integer containing a cid referring to a small molecule.
}
}
\value{
A \code{data.frame} where the row names represent each target the specified compound shows inactivity against, and the column shows the number of assays in which it was found to be inactive.
}
\author{
Tyler Backman
}
\seealso{
\code{\link{activeTargets}}
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get targets that compound 2244 shows inactivity against
myCidInactiveTargets <- row.names(inactiveTargets(sampleDB, "2244"))

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
