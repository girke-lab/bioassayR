\name{bioassay-class}
\Rdversion{1.1}
\docType{class}
\alias{bioassay}
\alias{bioassay-class}
\alias{aid,bioassay-method}
\alias{aid<-,bioassay-method}
\alias{assay_type,bioassay-method}
\alias{assay_type<-,bioassay-method}
\alias{organism,bioassay-method}
\alias{organism<-,bioassay-method}
\alias{scores,bioassay-method}
\alias{scores<-,bioassay-method}
\alias{scoring,bioassay-method}
\alias{scoring<-,bioassay-method}
\alias{show,bioassay-method}
\alias{source_id,bioassay-method}
\alias{source_id<-,bioassay-method}
\alias{target_types,bioassay-method}
\alias{target_types<-,bioassay-method}
\alias{targets,bioassay-method}
\alias{targets<-,bioassay-method}
\alias{aid}
\alias{aid<-}
\alias{assay_type}
\alias{assay_type<-}
\alias{organism}
\alias{organism<-}
\alias{scores}
\alias{scores<-}
\alias{scoring}
\alias{scoring<-}
\alias{show}
\alias{source_id}
\alias{source_id<-}
\alias{target_types}
\alias{target_types<-}
\alias{targets}
\alias{targets<-}

\title{Class \code{"bioassay"}}
\description{
This class represents the data from a bioassay experiment, where a number of small molecules
are screened against a defined target (such as a protein or living organism).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bioassay", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{aid}:}{Object of class \code{"character"} 
        containing the assay id. For assays sourced from NCBI PubChem, this should be a string containing the PubChem AID (assay identifier).
     }
    \item{\code{source_id}:}{Object of class \code{"character"}.
        This should match the description for a data source loaded via the \code{addDataSource()} function.
     }
    \item{\code{assay_type}:}{Object of class \code{"character"}.
        A string noting the type of bioactivity experiment, such as ``confirmatory'' to represent a confirmatory assay.
     }
    \item{\code{organism}:}{Object of class \code{"character"}.
        A string noting the scientific name of the assays target organism. }
    \item{\code{scoring}:}{Object of class \code{"character"}.
        A string noting the scoring method used for the bioactivity experiment. For example, IC50 or EC50. }
    \item{\code{targets}:}{Object of class \code{"character"}. A string or vector of strings
         containing the target identifier indicating the assay target. In the case of protein targeted assays sourced from NCBI PubChem, this should be a genbank ID. }
    \item{\code{target_types}:}{Object of class \code{"character"}.
        A string of text or vector of strings, representing (in the same order) the target types for each target. For example ``protein'' or ``cell.'' }
    \item{\code{scores}:}{Object of class \code{"data.frame"}
containing the bioactivity data to be loaded. This must be a 3 column data frame, with each row representing the bioactivity results of a single molecule. The first column represents the compound id (cid), which must be a unique value for each structurally distinct molecule. The second column is a binary value representing activity (1=active, 0=inactive, NA=inconclusive or untested) for the given assay. The last column represents a score, scored by the method specified with the \code{addBioassay()} function. Missing or non-applicable values in any column should be represented by a \code{NA} value. }
  }
}
\section{Methods}{
  \describe{
    \item{aid}{\code{signature(x = "bioassay")}: ... }
    \item{aid<-}{\code{signature(x = "bioassay")}: ... }
    \item{assay_type}{\code{signature(x = "bioassay")}: ... }
    \item{assay_type<-}{\code{signature(x = "bioassay")}: ... }
    \item{organism}{\code{signature(object = "bioassay")}: ... }
    \item{organism<-}{\code{signature(object = "bioassay")}: ... }
    \item{scores}{\code{signature(x = "bioassay")}: ... }
    \item{scores<-}{\code{signature(x = "bioassay")}: ... }
    \item{scoring}{\code{signature(x = "bioassay")}: ... }
    \item{scoring<-}{\code{signature(x = "bioassay")}: ... }
    \item{show}{\code{signature(object = "bioassay")}: ... }
    \item{source_id}{\code{signature(x = "bioassay")}: ... }
    \item{source_id<-}{\code{signature(x = "bioassay")}: ... }
    \item{target_types}{\code{signature(x = "bioassay")}: ... }
    \item{target_types<-}{\code{signature(x = "bioassay")}: ... }
    \item{targets}{\code{signature(x = "bioassay")}: ... }
    \item{targets<-}{\code{signature(x = "bioassay")}: ... }
	 }
}
\author{
Tyler Backman
}
\seealso{
Related classes: bioassaySet, bioAssayDB. 
}
\examples{
showClass("bioassay")

## create a new bioassay object from sample data
data(samplebioassay)
myassay <- new("bioassay",aid="1000", source_id="test", targets="116516899", 
    target_types="protein", scores=samplebioassay)
myassay
}
\keyword{classes}
