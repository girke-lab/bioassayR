\name{bioactivityFingerprint}
\alias{bioactivityFingerprint}
\title{
Create an \pkg{ChemmineR} \code{FPset} object that contains bioactivity
results for a given set of compounds and targets.
}
\description{
Returns a custom binary descriptor fingerprint for a given set of query cids and
target compounds, based on the activity data within a \code{bioassaySet} object.
}
\usage{
bioactivityFingerprint(bioassaySet, targets = FALSE, conflictResolver = "activesFirst")
}
\arguments{
\item{bioassaySet}{
  A \code{bioassaySet} object to generate fingeprints from. For a given compound set, this can be generated from a database using the \code{getBioassaySetByCids} function.
}
\item{targets}{
  An optional list of target id(s) to consider when creating the binary fingerprint. If a listed target is not in the bioassaySet, or has no active scores it will still be accepted, but create a fingerprint with all zeros for this location. The binary order of this list is preserved, so that direct comparison and combination of resulting \code{FPset} objects created with the same target list can be performed. If omitted, the target list for the \code{bioassaySet} object will be used, as returned by the \code{allTargets} function.
}
\item{conflictResolver}{
    Optionally allows users to choose how replicates (multiple assays sharing common compounds and targets) are resolved if they disagree. If 'activesFirst' any active score will take precedence over an inactive. If 'mode' the resulting score will be computed according to the statistical mode using \code{as.numeric(names(which.max(table(x))))}. Users can also optionally pass a function here which (for each cid/target pair) will receive a list of '2' (active) and '1' (inactive) values, and can then return any desired number as a summary to be included in the resulting table. For a large matrix, the default option 'activesFirst' offers the lowest computational overhead.
}
}
\value{
The returned object is a standard \pkg{ChemmineR} \code{FPset} object, and
can be used as described in the \pkg{ChemmineR} documentation. The order and number of binary bits for each compound can be set using the \code{targets} option, enabling the combination or comparison of multiple objects created with the same target list.
If a single compound has both active and inactive scores for the same target, it will
be resolved according to the \code{confictResolver} option.
}
\author{
Tyler William H Backman
}
\seealso{
Functions: \code{getBioassaySetByCids}, \code{getAssays}
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve all targets in database
targetList <- allTargets(sampleDB)

## get an activity fingerprint object for selected CIDs
queryCids <- c("2244", "3715", "2662", "3033", "133021", 
    "44563999", "44564000", "44564001", "44564002") 
myAssaySet <- getBioassaySetByCids(sampleDB, queryCids)
myFp <- bioactivityFingerprint(bioassaySet=myAssaySet)

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
