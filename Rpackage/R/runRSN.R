#' Robust Spline Normalisation
#'
#' A wrapper around the lumi::rsn() function
#'
#'
#' @param eset an ExpressionSet object
#' @param refSet an ExpressionSet object
#' @param minValue TODO
#' @param ... TODO

#' The function can be run multicore by running if registered via the doMC or doParallel package
#'
#' @examples 
#' set.seed(1)
#' #TODO!
#'
#' @export
#' @importFrom stats quantile
#' @importFrom foreach getDoParRegistered
#' @importFrom Biobase exprs exprs<-
#' @importFrom lumi rsn
runRSN <-
function(eset,refSet,minValue=0,...){
    dots<-list(...)
    if (!is.null(dots[["myHiddenArgument"]])){
        x<-eset
        x.clean<-x[which(!is.na(x) & x>minValue & names(x)%in%names(refSet))]
        mat<-cbind(refSet[names(x.clean)], x.clean)
        x.clean<-rsn(x.lumi=mat ,targetArray=1)
        xNew<- x*NA
        xNew[rownames(x.clean)]<-x.clean[,2]
        xNew[ x<=minValue]<-0
        return(xNew)
    }
    ##Requires input on log2 scale!
    if (any(is.na(exprs(refSet)))) {stop("refSet contains NA in expression!")}
    if ( quantile(exprs(eset),0.98,na.rm=T)>30 | quantile(exprs(refSet),0.98,na.rm=T)>30 ) {stop("Data must be provided on log2 scale" );}
    refSet<-exprs(buildReferenceSet(refSet,refIsLogged2=TRUE))[,1]    
    x<-exprs(eset);
    `%backend%` <- `%do%`

    if (getDoParRegistered()) `%backend%` <- `%dopar%`
    
    cIdx<-NULL
    y<-foreach(cIdx = seq_len(ncol(eset)) ) %backend% {
        runRSN(x[,cIdx], refSet=refSet,  minValue=minValue, myHiddenArgument=TRUE)
    }
    y<-do.call(cbind,y)
    exprs(eset)<-y
    return(eset)
}
