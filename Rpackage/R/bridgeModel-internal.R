.bridgeLinearModel.internal <-
function (data, intercept, weights, newFeatures, pThres = 1e-05, maxFeatPerUnknown = Inf) 
{
    fittedModels <- .findMapping.internal(data, names(weights), 
        newFeatures, pThres = pThres, maxFeatPerUnknown = maxFeatPerUnknown)
    newProbes <- unique(c(unlist(lapply(fittedModels, function(x) {
        x$predictors
    }))))
    newWeights <- rep(0, length(newProbes))
    names(newWeights) <- newProbes
    newIntercept <- intercept
    for (probe in names(fittedModels)) {
        newWeights[fittedModels[[probe]]$predictors] <- newWeights[fittedModels[[probe]]$predictors] + 
            fittedModels[[probe]]$newWeights * weights[probe]
        newIntercept <- newIntercept + fittedModels[[probe]]$newIntercept * 
            weights[probe]
    }
    return(list(intercept = newIntercept, weights = newWeights))
}

#' @import limma
#' @import foreach
#' @importFrom stats prcomp model.matrix coef lm
.findMapping.internal <-
function (data, oldFeatures, newFeatures, 
    pThres = 1e-05, maxFeatPerUnknown = Inf) 
{
    oldFeatures <- unique(oldFeatures)
    newFeatures <- unique(newFeatures)
    if (length(intersect(oldFeatures, colnames(data))) != 
        length(oldFeatures)) {
        stop("Not all oldFeatures could be matched with the data")
    }
    if (length(unique(colnames(data))) != ncol(data)) {
        stop("Not all colnames are uniquely named")
    }
    if (nrow(data) < 35) {
        stop("At least 35 samples required in data.")
    }
    knownProbes <- intersect(oldFeatures, newFeatures)
    unknownProbes <- setdiff(oldFeatures, newFeatures)
    overlap <- intersect(colnames(data), newFeatures)
    if (length(overlap) == 0) {
        stop("No overlapping features found between source and targetdata")
    }
    df.source <- data.frame(data)
    x<-NULL
    fittedModels_unkown <- foreach(x = unknownProbes) %dopar% 
        {
            dm <- model.matrix(~1 + data[, x])
            colnames(dm) <- c("Intercept", "x")
            tt <- topTable(eBayes(lmFit(t(data[, overlap]), 
                dm)), coef = "x", number = Inf)
            pThres_cur <- max(pThres, tt[, "P.Value"][2])
            tt <- tt[tt[, "P.Value"] <= pThres_cur, , drop = FALSE]
            predictors <- rownames(tt)
            if (length(predictors) == 0) {
                newWeights <- numeric()
                newIntercept <- mean(df.source[, make.names(x)])
            }
            else {
                aPCA <- prcomp(data[, predictors])
                aPCA.x <- aPCA$x[, seq_len(min(c(ncol(aPCA$x), 
                  maxFeatPerUnknown, ceiling(nrow(data)/4)))), 
                  drop = FALSE]
                model <- lm(x ~ 1 + ., data.frame(x = df.source[, 
                  make.names(x)], aPCA.x))
                newWeights <- cbind(0, aPCA$rotation[, seq_len(ncol(aPCA.x)), 
                  drop = FALSE]) %*% coef(model)
                newIntercept <- coef(model)[1] - 1 * sum(aPCA$center * 
                  newWeights)
            }
            list(newWeights = newWeights, newIntercept = newIntercept, 
                predictors = predictors)
        }
    names(fittedModels_unkown) <- unknownProbes
    fittedModels_known <- sapply(knownProbes, function(x) {
        list(newWeights = matrix(1, dimnames = list(x, NULL)), 
            newIntercept = matrix(0, dimnames = list(NULL, "(Intercept)")), 
            predictors = knownProbes)
    }, simplify = FALSE)
    names(fittedModels_known) <- knownProbes
    c(fittedModels_unkown, fittedModels_known)
}
