makeGroups <-
function (v, ngroups, randomize = FALSE, seed = NULL) 
{
    if (!is.null(seed)) {
        set.seed(seed)
    }
    mygroups <- list()
    for (k in c(ngroups:1)) {
        if (randomize & length(v) > 1) {
            v <- sample(v)
        }
        mygroups[[k]] <- v[1:c(length(v)/k)]
        v <- v[-c(1:c(length(v)/k))]
    }
    mygroups
}
