\name{bioassaySet-class}
\Rdversion{1.1}
\docType{class}
\alias{bioassaySet-class}
\alias{show,bioassaySet-method}

\title{Class \code{"bioassaySet"}}
\description{
This class stores a large number of bioactivity scores from multiple assays and experiments as a single sparse matrix.
}
\section{Objects from the Class}{
Objects can be created with several functions including \code{getAssays} and \code{getBioassaySetByCids}.
}
\section{Slots}{
  \describe{
    \item{\code{activity}:}{Object of class \code{"dgCMatrix"} a binary sparse matrix of assays vs compounds where 0 represents inactive and 1 represents activity }
    \item{\code{scores}:}{Object of class \code{"dgCMatrix"} numeric activity scores }
    \item{\code{targets}:}{Object of class \code{"dgCMatrix"} a matrix of the targets for each aid listed in the activity and scores matrix }
    \item{\code{sources}:}{Object of class \code{"data.frame"} data sources for each assay }
    \item{\code{replicates}:}{Object of class \code{"factor"} a factor signaling which assays are replicates }
    \item{\code{source_id}:}{Object of class \code{"integer"} the source_id for each assay as an integer }
    \item{\code{assay_type}:}{Object of class \code{"character"} the experiment type for each assay }
    \item{\code{organism}:}{Object of class \code{"character"} scientific name of each target species }
    \item{\code{scoring}:}{Object of class \code{"character"} scoring method used in the scores matrix }
    \item{\code{target_types}:}{Object of class \code{"character"} type of target used in assay }
  }
}
\section{Methods}{
No methods defined with class "bioassaySet" in the signature.
}
\author{
Tyler William H Backman
}

\examples{
showClass("bioassaySet")
}
\keyword{classes}
