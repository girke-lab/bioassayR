## generics for BioassayDB

setGeneric("queryBioassayDB", function(object, query){
    standardGeneric("queryBioassayDB")
})

## generics for bioassay
setGeneric(name="aid", def=function(x) standardGeneric("aid"))
setGeneric(name="aid<-", def=function(x, value) standardGeneric("aid<-"))

setGeneric(name="source_id", def=function(x) standardGeneric("source_id"))
setGeneric(name="source_id<-", def=function(x, value) standardGeneric("source_id<-"))

setGeneric(name="assay_type", def=function(x) standardGeneric("assay_type"))
setGeneric(name="assay_type<-", def=function(x, value) standardGeneric("assay_type<-"))

setGeneric(name="scoring", def=function(x) standardGeneric("scoring"))
setGeneric(name="scoring<-", def=function(x, value) standardGeneric("scoring<-"))

setGeneric(name="targets", def=function(x) standardGeneric("targets"))
setGeneric(name="targets<-", def=function(x, value) standardGeneric("targets<-"))

setGeneric(name="target_types", def=function(x) standardGeneric("target_types"))
setGeneric(name="target_types<-", def=function(x, value) standardGeneric("target_types<-"))

setGeneric(name="scores", def=function(x) standardGeneric("scores"))
setGeneric(name="scores<-", def=function(x, value) standardGeneric("scores<-"))

## generics for bioassaySet

setGeneric(name="activity", def=function(x) standardGeneric("activity"))
setGeneric(name="activity<-", def=function(x, value) standardGeneric("activity<-"))

setGeneric(name="organism", def=function(x) standardGeneric("organism"))
setGeneric(name="organism<-", def=function(x, value) standardGeneric("organism<-"))

setGeneric(name="sources", def=function(x) standardGeneric("sources"))
setGeneric(name="sources<-", def=function(x, value) standardGeneric("sources<-"))

setGeneric(name="scoring", def=function(x) standardGeneric("scoring"))
setGeneric(name="scoring<-", def=function(x, value) standardGeneric("scoring<-"))

