\name{activeAgainst}
\alias{activeAgainst}
\title{
Show compounds active against a specified target
}
\description{
Returns a \code{data.frame} of small molecule cids which show activity against a specified target. Each row name represents a cid which shows activity, and the total screens and the percent active are shown in their respective columns.
}
\usage{
activeAgainst(database, target)
}
\arguments{
  \item{database}{
A \code{BioAssaySet} database to query.
}
  \item{target}{
A string or integer containing a target_id referring to a target of interest.
}
}
\value{
A \code{data.frame} where the row names represent each compound showing activity against the specified target. The second column shows the number of distinct assays in which this cid was screened against the target, and the first column shows the percentage of these which exhibited activity.
}
\author{
Tyler Backman
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## get cids of compounds which show activity against target 116516899
myCids <- row.names(activeAgainst(sampleDB, "166897622"))

## disconnect from database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
