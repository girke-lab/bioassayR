\name{BioAssaySet-class}
\Rdversion{1.1}
\docType{class}
\alias{BioAssaySet-class}
\alias{queryBioassayDB,BioAssaySet-method}
\alias{show,BioAssaySet-method}

\title{Class \code{"BioAssaySet"}}
\description{
This class holds a connection to a bioassayR sqlite database.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{BioAssaySet("datbasePath")}.
}
\section{Slots}{
  \describe{
    \item{\code{database}:}{Object of class \code{"SQLiteConnection"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{queryBioassayDB}{\code{signature(object = "BioAssaySet")}: ... }
    \item{show}{\code{signature(object = "BioAssaySet")}: ... }
	 }
}
\author{
Tyler Backman
}
\examples{
showClass("BioAssaySet")
}
\keyword{classes}
