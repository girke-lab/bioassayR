\name{loadIdMapping}
\alias{loadIdMapping}
\title{
Load a target identifier mapping into a bioassayR database
}
\description{
Loads an identifier mapping for a bioassay target (stored in the database as an NCBI GI number) to another
protein target naming system. Common uses include UniProt identifiers, similarity clusters, and common names.
}
\usage{
loadIdMapping(database, target, category, identifier)
}
\arguments{
  \item{database}{
A writable \code{BioassayDB} database to insert data into.
}
  \item{target}{
A single protein target NCBI GI number.
}
  \item{category}{
The specified identifier type of the data being loaded, such as 'UniProt'.
}
  \item{identifier}{
A \code{character} object with the new identifier. This should be length one, and the function
should be re-ran to add multiple identifiers.
}
}
\references{
http://www.ncbi.nlm.nih.gov/protein NCBI Protein Database
http://www.uniprot.org UniProt Protein Database
}
\author{
Tyler Backman
}
\seealso{
\code{\link{translateTargetId}}
}
\examples{
## create sample database
myDatabaseFilename <- tempfile()
mydb <- newBioassayDB(myDatabaseFilename, indexed=FALSE)

## load a sample translation from GI 6686268 to UniProt P11712
loadIdMapping(mydb, "6686268", "UniProt", "P11712")

## get UniProt identifier(s) for GI Number 6686268
UniProtIds <- translateTargetId(mydb, "6686268", "UniProt")
UniProtIds

## disconnect from and delete sample database
disconnectBioassayDB(mydb)
unlink(myDatabaseFilename)
}
\keyword{ utilities }
