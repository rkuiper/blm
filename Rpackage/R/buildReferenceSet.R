#' @importFrom Biobase ExpressionSet featureData experimentData
buildReferenceSet <-
function (eset, name, refIsLogged2 = FALSE) 
{
    if (missingArg(name)) 
        name <- "referenceSet"
    m<-exprs(eset)
    if (!refIsLogged2){
        m <- log2(m)
    }

    m <- rowMeans(m) ##Means on log2 scale

    if (!refIsLogged2) { #Make sure to return to original scale
        m<-2^m; 
    }
    m <- matrix(m, ncol = 1, dimnames = list(names(m), name))
    newEset <- ExpressionSet(m, featureData = featureData(eset), 
        experimentData = experimentData(eset), annotation = eset@annotation)
    newEset
}
