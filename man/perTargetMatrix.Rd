\name{perTargetMatrix}
\alias{perTargetMatrix}
\title{
Collapse a \code{bioassaySet} object from multiple assays by combining assays with a common target
}
\description{
Creates a \code{sparseMatrix} object which has an activity value for each distinct target identifier rather than each distinct assay. Users can optionally choose how replicates are resolved. By default active scores always take preference over inactives: if any assay for a given target vs compound combination shows active, this combination will be marked active in the resulting object.
Either binary activity categories or scalar numeric scores can be used. When used with numeric data, this will create a Z-score compound vs. target matrix similar to High Throughput Screening Fingerprints (HTSFPs).
This function is not designed for single assays with multiple targets, and if they are present only one of the targets will be considered.
}
\usage{
perTargetMatrix(assays, inactives = TRUE, assayTargets = FALSE, 
    targetOrder = FALSE, summarizeReplicates = "activesFirst", 
    useNumericScores = FALSE)
}
\arguments{
  \item{assays}{
A \code{bioassaySet} object with data from multiple assays, some of which may share a common target. If used with \code{useNumericScores = TRUE} this should be the output of \code{scaleBioassaySet} if a Z-score matrix is desired.
}
  \item{inactives}{
    A logical value. Include both active and inactive scores. If FALSE only active scores are returned. This is only used if \code{useNumericScores = FALSE}, with
    numeric scores inactives are always considered.
  }
  \item{assayTargets}{
      Provide a custom merge table of target identifiers for each assay. For example, if you have clustered the targets of many assays into bins you can here merge by common clusters instead of distinct targets. This must be vector of class \code{character} with names that correspond to your assay ids (aids) and values that correspond to desired targets or clusteres in the resulting matrix. Names and targets should be represented as a character, even if they are numeric names.
      Note that if an assay contains multiple targets, only the first is used.
  }
  \item{targetOrder}{
      An optional \code{character} vector of desired target names in order. This will become the row names in the resulting sparse matrix in exact order, making it possible to coherently bind together sparse matricies of different compounds. If a target is omitted from this list it will be dropped in the result, and if an extra target is included it will show up with all '0' (untested) entries. 
  }
  \item{summarizeReplicates}{
    Optionally allows users to choose how replicates (multiple assays sharing common compounds and targets) are resolved if they disagree. If 'activesFirst' any active score will take precedence over an inactive. If 'mode' the resulting score will be computed according to the statistical mode using \code{as.numeric(names(which.max(table(x))))}. Users can also optionally pass a function here which (for each cid/target pair) will receive a vector of '2' (active) and '1' (inactive) values (if  \code{useNumericScores = TRUE}), and can then return any desired number as a summary to be included in the resulting table. For a large matrix, the default option 'activesFirst' offers the lowest computational overhead. When used with \code{useNumericScores = TRUE} the option
    'activesFirst' will keep only the replicate with the greatest absolute value.
    To average across all replicates, one can pass the R function \code{mean}.
  }
  \item{useNumericScores}{
  A logical value. Use numeric score rather than binary data to create a scalar compound vs. target
  matrix. When used with the output of \code{scaleBioassaySet} this creates
  a Z-score compound vs. target matrix.
    \code{NaN} values are replaced with zeros, as these usually represent
    assays summarized with \code{scaleBioassaySet} where all compounds
    had an identical value.
  }
}
\value{
When used with \code{useNumericScores = FALSE} 
a \code{sparseMatrix} which contains a value of 2 for each target vs compound combination which shows activity in at least one parent assay, a value of 1 for inactive combinations, and a value of zero for untested or ambiguous values. 
Note that this is different from older versions of bioassayR (1.6 and older) which used to return a value of 1 for actives and did not have the option to process inactives.
When used with \code{useNumericScores = TRUE} the raw numeric scores are returned, with replicates summarized as specified with the \code{summarizeReplicates} option.
}
\author{
Tyler William H Backman
}
\seealso{
Functions: \code{scaleBioassaySet}, \code{getAssays}, \code{bioactivityFingerprint}
}
\examples{
## connect to a test database
extdata_dir <- system.file("extdata", package="bioassayR")
sampleDatabasePath <- file.path(extdata_dir, "sampleDatabase.sqlite")
sampleDB <- connectBioassayDB(sampleDatabasePath)

## retrieve three assays
assays <- getAssays(sampleDB, c("673509","103","105"))
assays

## collapse assays into perTargetMatrix
targetMatrix <- perTargetMatrix(assays)
targetMatrix

## disconnect from sample database
disconnectBioassayDB(sampleDB)
}
\keyword{ utilities }
