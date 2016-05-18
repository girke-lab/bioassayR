\name{translateTargetId}
\alias{translateTargetId}
\title{
Translate a protein target identifier to another identifier system
}
\description{
Returns a \code{character} vector of the protein target identifiers
using the specified category (classification system). This is most often
used to translate NCBI Protein GI numbers (as provided with the pre-build
PubChem Bioassay database) into UniProt identifiers.
}
\usage{
translateTargetId(database, target, category, fromCategory = "GI")
}
\arguments{
  \item{database}{
A \code{BioassayDB} database to query.
}
  \item{target}{
A protein target identifier to query (as set by the category option, with a default of "GI").
}
  \item{category}{
The specified identifier type to return, such as 'UniProt'.
}
  \item{fromCategory}{
The identifier type of the query (default "GI").
}
}
\value{
A \code{character} vector of the protein target identifiers
of the category specified for the target specified. An NA is returned
if no matching values exist in the database.
}
\references{
http://www.ncbi.nlm.nih.gov/protein NCBI Protein Database
http://www.uniprot.org UniProt Protein Database
}
\author{
Tyler Backman
}

\seealso{
\code{\link{loadIdMapping}}
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

## get GI identifier(s) for UniProt ID P11712
GIs <- translateTargetId(mydb, "P11712", "GI", "UniProt")
GIs

## disconnect from and delete sample database
disconnectBioassayDB(mydb)
unlink(myDatabaseFilename)
}
\keyword{ utilities }
