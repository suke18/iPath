#' set up for the parallel computing for biocParallel.
#'
#' This function sets up the environment for parallel computing.
#' @param nproc number of processors
#' @param BPPARAM bpparameter from bpparam
#' @keywords internal
#' @return BAPPARAM settings
setUp_BPPARAM = function (nproc = 0, BPPARAM = NULL)
{
    if (is.null(BPPARAM)) {
        if (nproc != 0) {
            if (.Platform$OS.type == "windows") {
                result <- SnowParam(workers = nproc)
            }
            else {
                result <- MulticoreParam(workers = nproc)
            }
        }
        else {
            result <- bpparam()
        }
        return(result)
    }
    else {
        return(BPPARAM)
    }
}


#' remove genes with 0 sd
#'
#' This function helps remove non-informative genes.
#' @param Y is the expression matrix.
#' @importFrom stats sd
#' @keywords internal
#' @import stats
#' @return a processed matrix
rem_data = function(Y){
    rem_ids = which(apply(Y, 1, sd) == 0)
    if (length(rem_ids) == 0){
        return(Y)
    }else{
        return(Y[-rem_ids, ])
    }
}

