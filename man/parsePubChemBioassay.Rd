\name{parsePubChemBioassay}
\alias{parsePubChemBioassay}
\title{
Parse PubChem Bioassay Data
}
\description{
Parses a PubChem Bioassay experimental result from two required files (a csv file and an XML description) into a \code{bioassay} object.
}
\usage{
parsePubChemBioassay(aid, csvFile, xmlFile, duplicates = "drop")
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
