#' Bridge Linear Models
#'
#' By bridging linear models, a model can be described in terms of other features. 
#'
#' This function was developed specifically with the intention to convert a 
#' linear model known as SKY92 - which uses AffyMetrix probe IDs - into a model
#' that is described in terms of Ensembl Gene IDs. 
#' 
#' The AffyMetrix ProbeIDs reflect gene transcript expressions, while the Ensembl
#' gene IDs summarize the genes as a whole. Some genes might be targeted by
#' none, one or multiple probes, and similarly, some probes might correspond to
#' no, one or multiple Ensembl gene IDs. Although there might not be a one-to-one
#' relationship, there is most often a high correlation between features expressed
#' in terms of AffyMetrix Probe IDs and other representations such as Ensembl Gene
#' IDs. 
#'
#' The function uses the underlying covariance structure in the data to 
#' redistribute the weight of a model over new features.
#'
#' @param data A matrix with named features in the rows, and unique sample names in the columns
#' @param intercept A single numeric value indicating the model intercept
#' @param weights A numeric vector corresponding the model feature weights.
#' @param newFeatures A character vector, indicating which features to use for the bridged model.
#' @param pThres TODO
#' @param maxFeatPerUnknown TODO
#' @param seed TODO
#'
#' @examples 
#' set.seed(1)
#' #TODO!
#'
#' @references
#' Rowan Kuiper, Sonja Zweegman, Mark van Duin, Martin H. van Vliet, Erik H. van Beers, et.al; Prognostic and predictive performance of R-ISS with SKY92 in older patients with multiple myeloma: the HOVON-87/NMSG-18 trial. Blood Adv 2020; 4 (24): 6298–6309. doi:https://doi.org/10.1182/bloodadvances.2020002838
#' @export
#' @importFrom stats quantile
blm <-
function (data, intercept, weights, newFeatures, 
    pThres = 1e-05, maxFeatPerUnknown = Inf, seed=1) 
{
    res0 <- .bridgeLinearModel.internal(data, intercept, weights, newFeatures, pThres, maxFeatPerUnknown)
    newThr <- NULL
    if (!is.numeric(intercept) || length(intercept) != 1 || any(is.na(intercept))) {
        stop("Invalid value for argument 'intercept'")
    }
    if (!is.numeric(pThres) || length(pThres) != 1 || any(is.na(pThres)) || 
        pThres > 1 || pThres < 0) {
        stop("Invalid value for argument 'pThres'")
    }
    if (!is.numeric(maxFeatPerUnknown) || length(maxFeatPerUnknown) != 
        1 || any(is.na(maxFeatPerUnknown)) || maxFeatPerUnknown < 
        1) {
        stop("Invalid value for argument 'maxFeatPerUnknown'")
    }

    ##Run some cross validations to obtain an equal proportion >0 for before/after
    folds <- makeGroups(rownames(data), 5, seed = seed)
    pred <- rep(NA, nrow(data))
    names(pred) <- rownames(data)
    trueScore <- data[, names(weights)] %*% weights + intercept
    for (fIdx in seq_along(folds)) {
        cat(".")
        train <- scale(data[unlist(folds[-fIdx]), ])
        test <- scale(data[unlist(folds[fIdx]), ])
        res <- .bridgeLinearModel.internal(weights = weights, 
            intercept = intercept, data = train, newFeatures = newFeatures, 
            pThres = pThres, maxFeatPerUnknown = maxFeatPerUnknown)
        pred[folds[[fIdx]]] <- test[, names(res$weights)] %*% res$weights + res$intercept
    }
    newThr <- quantile(pred, mean(trueScore < 0))

    res0$intercept<-res0$intercept - newThr
    return(res0)
}
