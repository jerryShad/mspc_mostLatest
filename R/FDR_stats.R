# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .create_OUTP
#' @param peak GRanges objects
#' @param pAdjustMethod adjusted pvalue for multiple comparison
#' @param alpha threhold threshold value that keep the peaks below this value
#' @return GRanges
#' @export
#' @importFrom stats, p.adjust
#' @author Julaiti Shayiding
#' @example

.FDR.stats <- function(peaks, pAdjustMethod="BH", alpha=0.05, ...) {
  # input param checking
  stopifnot(class(peaks)=="GRanges")
  stopifnot(is.numeric(alpha))
  pAdjustMethod = match.arg(pAdjustMethod)
  if(is.null(peaks$p.value)) {
    stop("required slot is missing")
  } else {
    p <- peaks$p.value
    p.adj <- p.adjust(p, method = pAdjustMethod)
    peaks$p.adj <- p.adj
    rslt <- split(peaks, ifelse(peaks$p.adj <= alpha,
                                "pass", "fail"))
    return(rslt)
  }
}

#' example
## gr <- GRanges( seqnames=Rle("chr1", 4),ranges=IRanges(c(3,33,54,91), c(23,42,71,107)),
##               rangeName=c("a1", "a4", "a7", "a11"), p.value=c(1e-22, 1e-6,1e-13, 1e-7))

## .FDR.stats(peaks = gr, pAdjustMethod = "BH", alpha = 0.05)
