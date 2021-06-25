
#' GSEA calculation
#'
#' This function calculates the GSEA enrichment score.
#' @param gene_list is a list of genes.
#' @param gene_set is a set of genes.
#' @param stats_vector a vector quantify the level of genes in the gene list.
#' @return the orignial GSEA score.
GSEA = function(gene_list, gene_set, stats_vector){
    # Input genelist must be ordered.
    tag_indicator = sign(match(gene_list, gene_set, nomatch = 0))
    no_tag_indicator = 1- tag_indicator
    N = length(gene_list); Nh = length(gene_set); Nm =  N - Nh

    sum_rank_tag = sum(stats_vector[tag_indicator ==1])
    if (sum_rank_tag == 0){
        ES = -1
    }else{
        norm_tag = 1.0/ sum_rank_tag; norm_no_tag = 1/Nm
        RES = cumsum(tag_indicator*stats_vector*norm_tag - no_tag_indicator*norm_no_tag)
        max.ES = max(RES)
        min.ES = min(RES)
        if (max.ES > - min.ES) {
            ES = signif(max.ES, digits=5)
        } else {
            ES = signif(min.ES, digits=5)
        }
    }
    return(ES)
}

#' iES calculation Function
#'
#' This function calculates the iES matrix which is the core of iPath.
#' @importFrom matrixStats rowSds
#' @importFrom Rcpp sourceCpp
#' @import BiocParallel
#' @param Y is the expression matrix.
#' @param GSDB is the gene set database.
#' @param BPPARAM parameters from the BiocParallel.
#' @param nPro number of processors (default = 0).
#' @return a matrix with rows corresponding to the pathways and columns corresponding to the patients.
#' @keywords iES statistics calculatioin.
#' @examples
#' data(PRAD_data)
#' data(GSDB_example)
#' iES_mat = iES_cal2(prad_exprs, GSDB = GSDB_example)
#' @export
iES_cal2 = function(Y, GSDB,  BPPARAM = NULL, nPro = 0){
    Y = as.matrix(Y)
    npats = ncol(Y)
    GSDB_paths = GSDB$genesets
    GSDB_paths_names = GSDB$geneset.names
    ngsets = length(GSDB_paths);
    message("start normalization ...")
    Y = rem_data(Y)
    gnames = rownames(Y)
    row_mean = rowMeans(Y)
    row_sd = rowSds(Y)
    tmp_norm = abs((Y-row_mean) / row_sd)

    order_array = apply(tmp_norm, 2, function(x) rev(order(x)))
    order_name = apply(order_array, 2, function(i) gnames[i])
    order_stats = sapply(1:npats, function(i) tmp_norm[, i][order_array[, i]])

    bp_fun = function(i) {
        one_stats = order_stats[, i]
        one_names = order_name[, i]
        names(one_stats) = one_names
        one_pat_vec = sapply(1:ngsets, function(j){
            one_match_pos = na.omit(match(GSDB_paths[[j]], one_names))
            return(caliES2(one_stats, one_match_pos))
        })
    }

    tmpParam = setUp_BPPARAM(nPro, BPPARAM = BPPARAM)
    tmpParam$progressbar = TRUE
    message("start calculating iES matrix ...")
    pats_iESs = bplapply(seq_len(npats), FUN = bp_fun, BPPARAM = tmpParam)
    iES_mat = do.call(cbind, pats_iESs)
    dimnames(iES_mat) = list(GSDB_paths_names, colnames(Y))
    return(iES_mat)
}


