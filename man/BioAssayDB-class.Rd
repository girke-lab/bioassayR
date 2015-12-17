\name{BioassayDB-class}
\Rdversion{1.1}
\docType{class}
\alias{BioassayDB-class}
\alias{queryBioassayDB,BioassayDB-method}
\alias{show,BioassayDB-method}

\title{Class \code{"BioassayDB"}}
\description{
This class holds a connection to a bioassayR sqlite database.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{BioassayDB("datbasePath")}.
}
\section{Slots}{
  \describe{
    \item{\code{database}:}{Object of class \code{"SQLiteConnection"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{queryBioassayDB}{\code{signature(object = "BioassayDB")}: ... }
    \item{show}{\code{signature(object = "BioassayDB")}: ... }
	 }
}
\author{
Tyler Backman
}
\seealso{
Related classes: bioassaySet, bioassay. 
}
\examples{
showClass("BioassayDB")
}
\keyword{classes}
