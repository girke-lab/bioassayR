\name{scaleBioassaySet}
\alias{scaleBioassaySet}
\title{
Centers and standardizes the numeric activity scores for
a \code{bioassaySet} object (creates Z-scores)
}
\description{
Converts the numeric activity scores for a \code{bioassaySet} object
into per-assay Z-scores. Untested '0' values are not considered
in computing the value, only actives and inactives. In essence, this
is a special version of the R base \code{scale} function,
which ignores missing entries in a sparse matrix instead of using them
as zeros. A primary purpose of this function is to pass scaled
results to \code{perTargetMatrix}, in order to compute a numeric Z-score
compound vs. target matrix.
}
\usage{
scaleBioassaySet(bioassaySet, center=TRUE, scale=TRUE)
}
\arguments{
  \item{bioassaySet}{
A \code{bioassaySet} object with data from multiple assays. This should
be created with \code{getAssays} rather than \code{getBioassaySetByCids},
as the former includes the full assay data whereas the latter omits
scores for compounds other than those specified, and therefore will
not compute a coherent Z-score.
}
  \item{center}{
A logical value. If center is TRUE then centering is done by subtracting the assay means (omitting inconclusive NAs) from their corresponding scores, and if center is FALSE, no centering is done.
  }
  \item{scale}{
A logical value. Scaling is done by dividing the (centered) per-assay scores by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
  }
}
\value{
A \code{bioassaySet} object with standardized numeric scores, that can be
accessed with \code{scores(bioassaySet)}.
}
\author{
Tyler William H Backman
}
\seealso{
Functions: \code{getAssays}, \code{perTargetMatrix}
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve three assays
assays <- getAssays(sampleDB, c("347221","53211","624349"))

## disconnect from sample database
disconnectBioassayDB(sampleDB)

## compute and return standardized scores
scaledAssays <- scaleBioassaySet(assays)

## inspect scaled and unscaled scores
scores(assays)
scores(scaledAssays)

## NOTE: this example only returns non-NA Z-scores if tried with
## real data, not the test database used here
}
\keyword{ utilities }
