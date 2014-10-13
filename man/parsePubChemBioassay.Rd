\name{parsePubChemBioassay}
\alias{parsePubChemBioassay}
\title{
Parse PubChem Bioassay Data
}
\description{
Parses a PubChem Bioassay experimental result from two required files (a csv file and an XML description) into a \code{bioassay} object.
}
\usage{
parsePubChemBioassay(aid, csvFile, xmlFile)
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
