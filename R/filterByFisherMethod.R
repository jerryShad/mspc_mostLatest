## MSPC Project -- Bioconductor Packaaage for Multiple Sample Peak Calling
##
##' @title .filterByFisherMethod
##' @description
##' since we are getting combined p.value by Fisher methods, we need to set up second level filtering that
##' each peaks can be classified as confirmed or discarded . more correct description is needed.
##'
##' @param peakset
##' @param .ovHit
##' @param tau.s
##' @param comb.p
##' @param isFisherPass
##' @export
##' @importFrom XVector extractList
##' @author Julaiti Shayiding
##' @usage

.filterByFisherMethod <- function(peakset, .ovHit, tau.s, comb.p ,isFisherPass, ...) {
  if(isFisherPass) {
    lapply.func <- function(ele_) {
      keepMe <- sapply(comb.p, function(x) x<=tau.s)
      res <- ele_[keepMe]
    }
  } else {
    lapply.func <- function(ele_) {
      drop_ <- sapply(comb.p, function(x) x > tau.s)
      res <- ele_[drop_]
    }
  }
  .hitIdx <- lapply(.ovHit, lapply.func)
  .expandAsGR <- Map(unlist,
                     mapply(extractList, peakset, .hitIdx))
  .expandAsGR[[1L]] <- unique(.expandAsGR[[1L]])
  .expandAsGR <- setNames(.expandAsGR, names(peakset))
  return (.expandAsGR)
}

#' @example
# .Confirmed.ERs <- .filterByFisherMethod(total.ERs, keepList, tau.s=1.0E-08, comb.p, isFisherPass = TRUE)
# .Discarded.ERs <- .filterByFisherMethod(total.ERs, keepList, tau.s=1.0E-08, comb.p, isFisherPass = FALSE)


