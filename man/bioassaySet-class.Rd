\name{bioassaySet-class}
\Rdversion{1.1}
\docType{class}
\alias{bioassaySet-class}
\alias{show,bioassaySet-method}
\alias{activity,bioassaySet-method}
\alias{activity<-,bioassaySet-method}
\alias{scores,bioassaySet-method}
\alias{scores<-,bioassaySet-method}
\alias{targets,bioassaySet-method}
\alias{targets<-,bioassaySet-method}
\alias{sources,bioassaySet-method}
\alias{sources<-,bioassaySet-method}
\alias{source_id,bioassaySet-method}
\alias{source_id<-,bioassaySet-method}
\alias{assay_type,bioassaySet-method}
\alias{assay_type<-,bioassaySet-method}
\alias{organism,bioassaySet-method}
\alias{organism<-,bioassaySet-method}
\alias{scoring,bioassaySet-method}
\alias{scoring<-,bioassaySet-method}
\alias{target_types,bioassaySet-method}
\alias{target_types<-,bioassaySet-method}
\alias{activity}
\alias{activity<-}
\alias{sources}
\alias{sources<-}

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
    \item{\code{source_id}:}{Object of class \code{"integer"} the source_id for each assay as an integer. The length should equal the number of rows in the activity matrix, with element names for each assay id (aid). }
    \item{\code{assay_type}:}{Object of class \code{"character"} the experiment type for each assay.  The length should equal the number of rows in the activity matrix, with element names for each assay id (aid). }
    \item{\code{organism}:}{Object of class \code{"character"} scientific name of each target species. The length should equal the number of rows in the activity matrix, with element names for each assay id (aid).}
    \item{\code{scoring}:}{Object of class \code{"character"} scoring method used in the scores matrix. The length should equal the number of rows in the activity matrix, with element names for each assay id (aid). }
    \item{\code{target_types}:}{Object of class \code{"character"} type of target for each target id, where the names and order match the columns in the target matrix. The length should equal the number of rows in the activity matrix, with element names for each assay id (aid). }
  }
}
\section{Methods}{
  \describe{
    \item{activity}{\code{signature(x = "bioassaySet")}: ... }
    \item{activity<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{scores}{\code{signature(x = "bioassaySet")}: ... }
    \item{scores<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{targets}{\code{signature(object = "bioassaySet")}: ... }
    \item{targets<-}{\code{signature(object = "bioassaySet")}: ... }
    \item{sources}{\code{signature(x = "bioassaySet")}: ... }
    \item{sources<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{source_id}{\code{signature(x = "bioassaySet")}: ... }
    \item{source_id<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{assay_type}{\code{signature(object = "bioassaySet")}: ... }
    \item{assay_type<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{organism}{\code{signature(x = "bioassaySet")}: ... }
    \item{organism<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{scoring}{\code{signature(x = "bioassaySet")}: ... }
    \item{scoring<-}{\code{signature(x = "bioassaySet")}: ... }
    \item{target_types}{\code{signature(x = "bioassaySet")}: ... }
    \item{target_types<-}{\code{signature(x = "bioassaySet")}: ... }
     }
}
\author{
Tyler William H Backman
}
\seealso{
Related classes: bioassay, bioAssayDB. 
}
\examples{
showClass("bioassaySet")
}
\keyword{classes}
