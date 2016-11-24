# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .create_OUTP
#' @param peakset_1 list of confirmed peak intervals
#' @param peakset_2 list of discardedd peak intervals
#' @param pAdjustMethod adjusted pvalue for multiple comparison
#' @param .fdr parameter for false discovery rate
#' @param replicate.type type of replicate (Biological / Technical replicates)
#' @param outDir user can control where the exported BED file goes
#'
#' @return BED file return standard BED File
#' @export
#' @importFrom stats, p.adjust
#' @importFrom dplyr anti_join
#' @importFrom rtracklayer export.bed
#' @author Julaiti Shayiding
#' @example

.FDR.stats <- function(peakset_1, peakset_2, pAdjustMethod="BH", .fdr = 0.05
                       , replicate.type=c("Biological", "Technical"), outDir=getwd(), ...) {
  # check input param
  require(dplyr)
  #stopifnot(length(peakset)>0)
  pAdjustMethod = match.arg(pAdjustMethod)
  replicate.type = match.arg(replicate.type)
  stopifnot(is.numeric(.fdr))
  message("set purification on set of confirmed, discarded peaks")
  # peakset_1, peakset_2 must be casted to data.frame

  peakset_1 <- lapply(peakset_1, as.data.frame)
  peakset_2 <- lapply(peakset_2, as.data.frame)

  if (!dir.exists(outDir)) {
    dir.create(file.path(outDir))
  }
  if(replicate.type=="Biological") {
    .setPurf <- peakset_1
  } else {
    .setPurf <- Map(anti_join, peakset_1, peakset_2)
  }
  .setPurf <- lapply(.setPurf, function(x) as(x, "GRanges"))
  res <- lapply(.setPurf, function(ele_) {
    if(is.null(ele_$p.value)) {
      stop("p.value is required")
    } else {
      p <- ele_$p.value
      ele_$p.adj <- p.adjust(p, method = "BH")
      .filt <- split(ele_,
                     ifelse(ele_$p.adj <= .fdr,
                            "BH_Pass", "BH_Failed"))
    }
  })
  res <- lapply(res, unique)
  rslt <- lapply(names(res), function(ele_) {
    mapply(export.bed,
           res[[ele_]],
           paste0(outDir, ele_, ".", names(res[[ele_]]), ".bed"))
  })
  return(rslt)
}

#' @example
fdr.rslt <- .FDR.stats(.Confirmed.ERs, .Discarded.ERs, pAdjustMethod = "BH",
                       .fdr = 0.05, replicate.type = "Bio", outDir = "test/")
