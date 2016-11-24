##' MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##'
##' @title .Fisher.stats
##' @param .ovHit overlap hit index for sufficiently overlaped peaks
##' @param peakset list of peak interval
##' @return combine pvalue
##' @export
##' @importFrom metap sumlog
##' @importFrom XVectors extractList
##' @author Julaiti Shayiding
##' @example

.Fisher.stats <- function(hitTB, peakset, verbose=FALSE, ...) {
  # input param checking
  require(metap)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  message("retrieve pvalue of peaks")
  pval_List <- mapply(.get.pvalue, hitTB, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0.000000e+00)
    })
  }
  pval.TB <- as.data.frame(mapply(.helper.PV, pval_List))
  # pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  # Npresent <- rowSums( !is.na(pval.TB) )
  # comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  comb.pval <- suppressWarnings(
    .globSC <- apply(pval.TB[, 1:length(pval.TB)], 1, function(ro) {
      res <- sumlog(ro)$p
    })
  )
  comb.pval <- as.matrix(comb.pval)
  return(comb.pval)
}

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}
