\name{samplebioassay}
\alias{samplebioassay}
\docType{data}
\title{
Sample activity data for use with bioassayR
}
\description{
This is sample bioactivity data, taken from assay identifier (aid) 1000 in the NCBI PubChem Bioassay database. These data are provided for testing the bioassayR library.
}
\usage{data(samplebioassay)}
\format{
  A data frame with activity scores for 4 distinct compounds.
  \describe{
    \item{\code{cid}}{unique compound identifer}
    \item{\code{activity}}{1=active, 0=inactive, NA=other}
    \item{\code{score}}{activity scores}
  }
}
\source{
http://pubchem.ncbi.nlm.nih.gov NCBI PubChem
}
\references{
http://pubchem.ncbi.nlm.nih.gov NCBI Pubchem
}
\examples{
## create a new bioassay object from these sample data 
data(samplebioassay)
myassay <- new("bioassay",aid="1000", source_id="PubChem BioAssay", targets="116516899", 
    target_types="protein", scores=samplebioassay)
myassay
}
\keyword{datasets}
