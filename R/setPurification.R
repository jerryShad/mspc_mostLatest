# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .setPurification
#' @param peakList list of all confirmed and discarded peaks for each replicates
#' @param replicate.type type of replicate whether Biological or Technical
#' @return GRanges
#' @export
#' @importFrom dplyr anti_join
#' @author Julaiti Shayiding

.setPurification <- function(.grsConf, .grsDisc, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(inherits(grs[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  dfConf <- Map(as.data.frame, .grsConf)
  dfDisc <- Map(as.data.frame, .grsDisc)
  .setPurification <- ifelse(replicate.type=="Biological",
                             rslt <- dfConf,
                             rslt <- Map(anti_join, dfConf, dfDisc))
  .res <- lapply(.setPurification, function(x) as(x, "GRanges"))
  return(.res)
}

#' @example
#' res.STPR <- .setPurification(confirmedList, discardedList,"Biological")



