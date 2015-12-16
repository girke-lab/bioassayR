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
    \item{\code{activity}:}{Object of class \code{"dgCMatrix"} a sparse matrix of assays (rows) vs compounds (columns) where 0 represents untested, NA represents inconclusive, 1 represents inactive, and 2 represents activity }
    \item{\code{scores}:}{Object of class \code{"dgCMatrix"} numeric activity scores with the same dimensions as activity }
    \item{\code{targets}:}{Object of class \code{"dgCMatrix"} a binary matrix of the targets (columns) for each aid (rows) listed in the activity and scores matrix. A 1 represents a target for the given assay, and a 0 represents that the given target was not used in the assay. }
    \item{\code{sources}:}{Object of class \code{"data.frame"} data sources for each assay. There must be three columns titled 'source_id', 'description', and 'version.' Each row represents a data source for these data. The 'source_id' must be a numeric (integer) index that matches to those in the 'source_id' slot.}
    \item{\code{replicates}:}{Object of class \code{"factor"} a factor signaling which assays are replicates. This will have a length and order equal to the row names (assay aids) in the activity matrix, and will repeat a label for any pair of assays that are replicates. Non-replicate assays should have a unique value for each assay. If no replicates are present among the assay, this can optionally be a factor length zero.}
    \item{\code{source_id}:}{Object of class \code{"integer"} the source_id for each assay as an integer, in the same order as the row names of the activity matrix. }
    \item{\code{assay_type}:}{Object of class \code{"character"} the experiment type for each assay in the same order as the row names of the activity matrix. }
    \item{\code{organism}:}{Object of class \code{"character"} scientific name of each target species. }
    \item{\code{scoring}:}{Object of class \code{"character"} scoring method used in the scores matrix }
    \item{\code{target_types}:}{Object of class \code{"character"} type of target for each target id, where the names and order match the columns in the target matrix. }
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
