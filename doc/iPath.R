## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----quickYo, eval = FALSE----------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("iPath")

## ----load_data, eval = TRUE,  message=FALSE, results='hide', include = TRUE----
library(iPath)
data(PRAD_data)
dim(prad_exprs)
data(GSDB_example)
head(prad_cli)

## ----cal_iES2, eval = TRUE, message=FALSE, results='hide', include = TRUE-----
iES_mat = iES_cal2(Y = prad_exprs, GSDB = GSDB_example)
iES_mat[1:2, 1:4]

## ----survival, eval = TRUE, message=FALSE, include = TRUE---------------------
surv_outcomes = iES_surv(iES_mat = iES_mat, cli = prad_cli, indVec = prad_inds)
head(surv_outcomes)

## ----waterfall, eval = TRUE, message=FALSE, include = TRUE--------------------
water_fall(iES_mat, gs_str = "SimPathway2", indVec = prad_inds)
density_fall(iES_mat, gs_str = "SimPathway2", indVec = prad_inds)

## ----onesurvival, eval = TRUE, message=FALSE, include = TRUE------------------
iES_survPlot(iES_mat = iES_mat, cli = prad_cli, gs_str = "SimPathway1", indVec = prad_inds, title = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

