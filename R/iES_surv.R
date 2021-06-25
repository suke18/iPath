#' iES calculation Function
#'
#' This function allows to investigate on one specific pathway.
#' @import mclust
#' @import survival
#' @param iES_mat is iES matrix with rows corresponding to the pathway and columns corresponding to the patients.
#' @param cli clinical data associated to the gene expression data.
#' @param indVec binary vector indicating normal (0) and tumor (1).
#' @keywords iPath survival analysis for two groups of patients: perturbed and normal-like.
#' @return a matrix of survival analysis from coxph.
#' @param npatsThre the threshold of number of patients for survival analysis.
#' @export
#' @examples
#' data(PRAD_data)
#' data(GSDB_example)
#' iES_mat = iES_cal2(prad_exprs, GSDB = GSDB_example)
#' iES_surv(iES_mat, cli = prad_cli, indVec = prad_inds)
iES_surv = function(iES_mat, cli, indVec = NULL, npatsThre = 5){
    if (is.null(indVec) | length(indVec)!= ncol(iES_mat)){
        stop("input right binary vector indicating the patients phenotypes")
    }

    npaths = nrow(iES_mat)
    path_names = rownames(iES_mat)
    inds1 = which(indVec==0)
    inds2 = which(indVec==1)
    norm_Y = iES_mat[, inds1]
    tumor_Y = iES_mat[, inds2]
    tumor_pat_names = colnames(tumor_Y)
    tumor_com = intersect(tumor_pat_names, cli$bcr_patient_barcod)
    tumor_Y = tumor_Y[, which(tumor_pat_names %in% tumor_com)]

    iPath_res = sapply(1:npaths, function(i){
        norm_vec = norm_Y[i, ]
        tumor_vec = tumor_Y[i, ]
        tmp_m = Mclust(norm_vec, parameter = TRUE, modelNames = "V")
        id = which.max(tmp_m$parameters$pro)
        tmp_mean = tmp_m$parameters$mean[id]
        tmp_sd = sqrt(tmp_m$parameters$variance$sigmasq[id])
        # UP or DOWN Regulate
        if (tmp_mean < mean(tumor_vec)){
            thre =  tmp_mean + 2*tmp_sd
            perturb = names(tumor_vec)[which(tumor_vec >= thre)]
        }else{
            thre =  tmp_mean - 2*tmp_sd
            perturb = names(tumor_vec)[which(tumor_vec < thre)]
        }
        normlike = setdiff(tumor_com, perturb)
        nlen1 = length(perturb)
        nlen2 = length(normlike)
        if (nlen1 >= npatsThre & nlen2 >= npatsThre){
            tmp_cli = rbind(cli[which(cli$bcr_patient_barcode %in% perturb),],
                            cli[which(cli$bcr_patient_barcode %in% normlike),])
            tmp_cli$class = c(rep("perturb", nlen1), rep("normlike", nlen2))
            tmp_cox = coxph(Surv(times, patient.vital_status) ~ class, data=tmp_cli, method = "breslow")
            nperturb = length(perturb)
            c = summary(tmp_cox)$concordance[1]
            coef = tmp_cox$coefficients
            pval = summary(tmp_cox)$sctest["pvalue"]
            tmp_res = as.numeric(c(nperturb, c, coef, pval))
        }else{
            tmp_res = rep(NA, 4)
        }
        return(tmp_res)
    })
    iPath_res = t(iPath_res)
    dimnames(iPath_res) = list(path_names, c("nPerturb", "c-index", "coef", "pval"))
    return(iPath_res)
}

