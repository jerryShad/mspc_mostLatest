# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .create_OUTP
#' @param peakset GRangeList or list of data.frame objects
#' @param pAdjustMethod adjusted pvalue for multiple comparison
#' @param alpha threhold threshold value that keep the peaks below this value
#' @return GRangeList
#' @export
#' @importFrom stats, p.adjust
#' @importFrom rtracklayer export.bed
#' @author Julaiti Shayiding
#' @example

.FDR.stats <- function(peakset, pAdjustMethod="BH", threshold=0.05, ...) {
  # check input param
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1L]], "GRanges"))
  pAdjustMethod = match.arg(pAdjustMethod)
  stopifnot(is.numeric(threshold))
  res <- lapply(peakset, function(ele_) {
    if(is.null(ele_$p.value)) {
      stop("p.value is required")
    } else {
      o <- ele_$p.value
      p.adj <- p.adjust(p, method = "BH")
      ele_$p.adj <- p.adj
      .filt <- split(ele_,
                     ifelse(ele_$p.adj <= threshold,
                            "BH_Pass", "BH_Failed"))
      return(.filt)
    }
  })
  rslt <- lapply(names(res), function(ele_) {
    mapply(export.bed,
           res[[ele_]],
           paste0(ele_, "_", names(res[[ele_]]), ".bed"))
  })
  rslt <- lapply(rslt, unique)
  return(rslt)
}

#' example
## .FDR.stats(peakset = grs, pAdjustMethod = "BH", alpha = 0.05, ...)
