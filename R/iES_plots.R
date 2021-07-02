########## The main demonstration plots #############

#' water fall plot
#'
#' This function allows you to express your love of cats.
#' @param iES_mat,gs_str is the iES_mat with tumor and normal and gs name.
#' @param indVec the binary indicator for normal(0) and tumor (1) patients.
#' @param title boolean true or false for including the title (gs_str) in the ggplot.
#' @import ggplot2
#' @import ggpubr
#' @return ggplot object containing the KM plot.
#' @keywords waterfall plot for normal and tumor sample
#' @examples
#' data(PRAD_data)
#' data(GSDB_example)
#' iES_mat = iES_cal2(prad_exprs, GSDB = GSDB_example)
#' water_fall(iES_mat, gs_str = "SimPathway1", indVec =prad_inds)
#' @export
water_fall = function(iES_mat, gs_str,indVec, title = TRUE){
    group_colors = c(tumor = "Brown", normal = "#56B4E9")
    inds1 = which(indVec==0)
    inds2 = which(indVec==1)
    iES_res = list(normal = iES_mat[, inds1], tumor = iES_mat[, inds2])

    gs_ind = which(rownames(iES_mat) == gs_str)
    normal = iES_res[[1]][gs_ind,]; tumor = iES_res[[2]][gs_ind, ]

    n_gap = round((length(tumor) + length(normal)) * 0.01)
    tmp_iES_mat = data.frame(value = c(tumor[order(tumor)], rep(0, n_gap), normal[order(normal)]),
                          type = c(rep("tumor", length(tumor)), rep("tumor", n_gap), rep("normal", length(normal))),
                          fill = c(rep("fill", length(tumor)), rep(NA, n_gap), rep("fill", length(normal))))
    if (title ==TRUE){
        p = ggplot(tmp_iES_mat, aes(x = 1:nrow(tmp_iES_mat), y = value)) +
            geom_area(aes(fill = type)) +
            theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            labs(x="Samples", y="Enrichment Score") +
            scale_fill_manual(values=group_colors) +
            ggtitle(gs_str)
    }else{
        p = ggplot(tmp_iES_mat, aes(x = 1:nrow(tmp_iES_mat), y = value)) +
            geom_area(aes(fill = type)) +
            theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            labs(x="Samples", y="Enrichment Score") +
            scale_fill_manual(values=group_colors)
    }
    return(p)
}

#' density fall plot
#'
#' This function allows you to express your love of cats.
#' @param iES_mat,gs_str is the iES_mat with tumor and normal and gs name.
#' @param indVec the binary indicator for normal(0) and tumor (1) patients.
#' @param title boolean true or false for including the title in the ggplot.
#' @import ggpubr
#' @keywords densityfall plot for normal and tumor sample.
#' @return ggplot object containing the KM plot.
#' @export
#' @examples
#' data(PRAD_data)
#' data(GSDB_example)
#' iES_mat = iES_cal2(prad_exprs, GSDB = GSDB_example)
#' density_fall(iES_mat, gs_str = "SimPathway1", indVec = prad_inds)
density_fall = function(iES_mat, gs_str,indVec, title = TRUE){
    inds1 = which(indVec==0)
    inds2 = which(indVec==1)
    iES_res = list(normal = iES_mat[, inds1], tumor = iES_mat[, inds2])

    gs_ind = which(rownames(iES_mat) == gs_str)
    normal = iES_res[[1]][gs_ind,]; tumor = iES_res[[2]][gs_ind, ]
    tmp_m = Mclust(normal, parameter = TRUE, modelNames = "V")
    id = which.max(tmp_m$parameters$pro)
    tmp_mean = tmp_m$parameters$mean[id]
    tmp_iES_mat = data.frame(value = c(normal, tumor),
                          type = c(rep("normal", length(normal)), rep("tumor", length(tumor))))
    if (title ==TRUE){
        p = ggdensity(tmp_iES_mat, x = "value",rug = TRUE,
                      color = "type", fill = "type",
                      palette = c("#56B4E9", "Brown"),
                      main = gs_str, legend.title = "") +
            geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
    }else{
        p = ggdensity(tmp_iES_mat, x = "value", rug = TRUE,
                      color = "type", fill = "type",
                      palette = c("#56B4E9", "Brown"),
                      legend.title = "")+
            geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
    }
    return(p)

}


#' iES survival for a certain pathway
#'
#' This function allows you to express your love of cats.
#' @import survminer
#' @import survival
#' @import mclust
#' @import ggpubr
#' @param iES_mat,gs_str is the GSDB iES_mat with tumor and normal and gs name.
#' @param cli clinical data corresponding to the expression data.
#' @param indVec the binary indicator for normal(0) and tumor (1) patients.
#' @param npatsThre the threshold of number of patients for survival analysis.
#' @param title boolean true or false for including the title (gs_str) in the ggplot.
#' @return ggplot object containing the KM plot.
#' @keywords densityfall plot for normal and tumor sample
#' @export
#' @examples
#' data(PRAD_data)
#' data(GSDB_example)
#' iES_mat = iES_cal2(prad_exprs, GSDB = GSDB_example)
#' iES_survPlot(iES_mat, cli = prad_cli, gs_str = "SimPathway1", indVec = prad_inds)
iES_survPlot = function(iES_mat, cli, gs_str, indVec = NULL, npatsThre=5, title = TRUE){
    npaths = nrow(iES_mat)
    path_names = rownames(iES_mat)
    inds1 = which(indVec==0)
    inds2 = which(indVec==1)
    iES_res = list(normal = iES_mat[, inds1], tumor = iES_mat[, inds2])

    norm_Y = iES_mat[, inds1]
    tumor_Y = iES_mat[, inds2]
    tumor_pat_names = colnames(tumor_Y)
    tumor_com = intersect(tumor_pat_names, cli$bcr_patient_barcod)
    tumor_Y = tumor_Y[, which(tumor_pat_names %in% tumor_com)]

    norm_vec = norm_Y[gs_str, ]
    tumor_vec = tumor_Y[gs_str, ]
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

    if (nlen1>=npatsThre & nlen2>=npatsThre){
        tmp_cli = rbind(cli[which(cli$bcr_patient_barcode %in% perturb),],
                        cli[which(cli$bcr_patient_barcode %in% normlike),])
        tmp_cli$class = c(rep("perturb", nlen1), rep("normlike", nlen2))
        tmp_cox = coxph(Surv(times, patient.vital_status) ~ class, data=tmp_cli, method = "breslow")
        nperturb = length(perturb)
        sfit <- survfit(Surv(times, patient.vital_status)~class, data = tmp_cli, conf.type="log-log")

        if (title ==TRUE){
            p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=FALSE,
                           legend.labs=c("perturbed", "normal-like"), legend.title="",
                           palette=c("darkred", "darksalmon"),
                           linetype = "strata",
                           title= gs_str, data = tmp_cli)
        }else{
            p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=FALSE,
                           legend.labs=c("perturbed", "normal-like"), legend.title="",
                           palette=c("darkred", "darksalmon"),
                           linetype = "strata", data = tmp_cli)
        }
    }else{
        p = NULL
    }
    return(p)
}

