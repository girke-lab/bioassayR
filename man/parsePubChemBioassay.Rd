\name{parsePubChemBioassay}
\alias{parsePubChemBioassay}
\title{
Parse PubChem Bioassay Data
}
\description{
Parses a PubChem Bioassay experimental result from two required files (a csv file and an XML description) into a \code{bioassay} object.
}
\usage{
parsePubChemBioassay(aid, csvFile, xmlFile, duplicates = "drop",
    missingCid = "drop", scoreRegex = "inhibition")
}
\arguments{
  \item{aid}{
The assay identifier (aid) for the assay to be parsed.
}
  \item{csvFile}{
A CSV file for a given assay, as downloaded from PubChem Bioassay.
}
  \item{xmlFile}{
An XML description file for a given assay, as downloaded from PubChem Bioassay.
}
  \item{duplicates}{
Specifies how duplicate CIDs in the same assay are treated. If 'drop' is specified,
only the first of each duplicated cid is kept and a warning is returned. If 'FALSE'
processing will stop with an error if duplicates are present. If 'TRUE' duplicates
will be included without warning, which may cause erroneous results with other
bioassayR functions that assume a unique cid list for each assay.
}
  \item{missingCid}{
A value of either 'drop' or a logical value of FALSE. If 'FALSE' processing
will stop with an error for any input compounds with an empty cid string.
If 'drop' is specified, a warning will be issued and these compounds will be skipped.
}
  \item{scoreRegex}{
A regular expression (perl compatible, case insensitive) to be matched to the column names
in the CSV header, to identify relavent score rows. If any rows match this regex, the
first matching row will be used in place of the 'PUBCHEM_ACTIVITY_SCORE' and it's row
name will be stored as the assays scoring method. The default will identify
most PubChem Bioassays which contain protein target inhibition data.
}
}
\value{
A \code{bioassay} object containing the loaded data.
}
\references{
http://pubchem.ncbi.nlm.nih.gov NCBI PubChem
}
\author{
Tyler Backman
}
\examples{
## get sample data locations
extdata_dir <- system.file("extdata", package="bioassayR")
assayDescriptionFile <- file.path(extdata_dir, "exampleAssay.xml")
activityScoresFile <- file.path(extdata_dir, "exampleScores.csv")

## parse files
myAssay <- parsePubChemBioassay("1000", activityScoresFile, assayDescriptionFile)
myAssay
}
\keyword{ utilities }
