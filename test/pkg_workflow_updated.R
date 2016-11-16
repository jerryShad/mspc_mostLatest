## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##' @description
##' all updated function code and workflow of MSPC Packages
##'
##==============================================================================================
##' @details read bed files as GRanges objects from file directory
##' function implementation

readPeakFiles <- function(peakFolder, pvalueBase=1L,verbose=FALSE, ...) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  stopifnot(length(peakFolder)>0)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      .gr <- as(import.bed(ele_), "GRanges")
      if(is.null(.gr$p.value)) {
        .gr <- .pvalueConversion(.gr, 1L)
      }
    }), tools::file_path_sans_ext(basename(files))
  )
  return(f.read)
}

##' @examples
myData <- readPeakFiles(peakFolder = "test/testData/", pvalueBase = 1L)

#-----------------------------------------------------------------------------------------------
#' @description pvalueConversion
.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$p.value)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  return(x)
}

##==============================================================================================
##' @details splitting GRanges
##' function implementation

.denoise.ERs <- function(grs, tau.w= 1.0E-04, .fileName="", outDir=getwd(), verbose=FALSE, ...) {
  # check input param
  stopifnot(class(grs[[1L]])=="GRanges")
  stopifnot(length(grs)>0)
  stopifnot(is.numeric(tau.w))
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if(!dir.exists(outDir)) {
    dir.create(file.path(outDir))
    setwd(file.path(outDir))
  }
  res <- lapply(seq_along(grs), function(x) {
    .gr <- grs[[x]]
    .grNM <- names(grs)[x]
    .drop <- .gr[.gr$p.value > tau.w]
    export.bed(.drop, sprintf("%s/%s.%s.bed", outDir, .fileName, .grNM))
    .keep <- .gr[.gr$p.value <= tau.w]
    return(.keep)
  })
  rslt <- setNames(res, names(grs))
  return(rslt)
}

#' @example
total.ERs <- .denoise.ERs(myData, tau.w = 1.0E-04, .fileName = "noisePeak", outDir = "test/")

##==============================================================================================
##' @details peakOverlapping
##' function implementation
## I am on the right track, keep coding *_*

.peakOverlapping <- function(peakset, FUN=which.max, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  res <- list()
  for(i in seq_along(peakset)) {
    que <- peakset[[i]]
    queHit <- as(findOverlaps(que), "List")
    supHit <- lapply(peakset[- i], function(ele_) {
      ans <- as(findOverlaps(que, ele_), "List")
      out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
      out.idx0 <- out.idx0[!is.na(out.idx0),]
      ans <- ans[out.idx0]
      ans
    })
    res[[i]] = DataFrame(c(list(que=queHit), sup=supHit))
    names(res[[i]]) = c(names(peakset[i]),names(peakset[- i]))
  }
  rslt <- lapply(res, function(x) as.matrix(x[names(res[[1]])]))
  rslt <- DataFrame(rbind(rslt[[1]],
                          unique(do.call("rbind", rslt[2: length(rslt)]))))
  rslt <- lapply(rslt, function(x) as(x, "CompressedIntegerList"))
  return(rslt)
}

#' @example
#' Hit <- .peakOverlapping(peakset = total.ERs, FUN = which.max)

##==========================================================================================================
MSPC.Analyzer <- function(peakset, ovHit, replicate.type=c("Biological","Technical"), tau.s=1.0E-08, ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(tau.s))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  length(peakset)-1,
                  length(peakset))

  cnt.ovHit <- as.matrix(Reduce('+', lapply(ovHit, lengths)))
  keepMe <- cnt.ovHit >= min.c
  ##
  dropList <- lapply(ovHit, function(ele_) ele_[!keep_me])
  init.discardPeaks <- Map(unlist,
                           mapply(extractList, peakset, dropList))
  ##
  keepList <- lapply(ovHit, function(x) x[keepMe])
  pval_List <- mapply(.get.pvalue, keepList, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
  }
  pval.TB <- as.data.frame(Map(.helper.PV, pval_List))
  #pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  #Npresent <- rowSums( !is.na(pval.TB) )
  #comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  comb.pval <- suppressWarnings(
    .globSC <- apply(pval.TB, 1, function(ro) {
      require(metap)
      res <- sumlog(ro)$p
    })
  )
  comb.pval <- as.matrix(comb.pval)
  ##--------------------------------------------------------------------------------
  ## TODO BEGIN : FIXME: to make more compatible
  Confirmed_idx <- lapply(keepList, function(elm) {
    saved <- sapply(comb.pval, function(x) x <= tau.s)
    res <- elm[saved]
  })
  .Confirmed.ERs <- Map(unlist,
                        mapply(extractList, peakset, Confirmed_idx))
  .Confirmed.ERs[[1L]] <- unique(.Confirmed.ERs[[1L]])
  .Confirmed.ERs <- setnames(.Confirmed.ERs, names(total.ERs))
  # TODO: export .confirmed.ERs as BED file into desired directory

  ##================================================================================
  Discarded_idx <- lapply(keepList, function(elm) {
    droped <- sapply(comb.pval, function(x) x > tau.s)
    res <- elm[droped]
  })
  .Fisher.discPeaks <- Map(unlist,
                           mapply(extractList, peakset, Discarded_idx))
  .Discarded.ERs <- suppressWarnings(mapply(c, .init.discPeaks, .Fisher.discPeaks))
  .Discarded.ERs[[1L]] <- unique(.Discarded.ERs[[1L]])
  .Discarded.ERs <- setNames(.Discarded.ERs, names(.Fisher.discPeaks))
  # TODO : export .discarded.ERs as BED File into desired directory

  ##-------------------------------------------------------------------------------
  ## Note: to use anti_join function, needed to be casted as data.frame
  .setPurification <- ifelse(replicate.type=="Biological",
                             res <- .Confirmed.ERs,
                             res <- Map(anti_join, .Confirmed.ERs, .Discarded.ERs))
  .setPurification <- lapply(.setPurification, function(x) as(x, "GRanges"))

  #--------------------------------------------------------------------------------
  # implementation for BH correction test
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
  ## TODO END
}

.get.pvalue <- function(hit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, hit)
  return(res)
}

##================================================================================================
.setPurification <- ifelse(replicate.type=="Biological",
                           res <- .Confirmed.ERs,
                           res <- Map(anti_join, .Confirmed.ERs, .Discarded.ERs))

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

.create_OUTP <- function(peaks, pAdjustMethod="BH", alpha=0.05, ...) {
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

#' @example
.BH_output <- Map(.create_OUTP, .setPurif)
.BH_output <- lapply(.BH_output, unique)
##======================================================================================================
## generate list of desired output BED files

both <- do.call("rbind", c(.Confirmed.ERs, .Discarded.ERs))
cn <- c("letter", "confirmed", "seq")
DF <- cbind(read.table(text = chartr("_", ".", rownames(both)), sep = ".", col.names = cn), both)

DF <- transform(DF, stringency = ifelse(p.value <= tau.s, "Stringent", "Weak"))

.res.out <- by(DF, DF[c("letter", "stringency", "confirmed")],
           function(x) export.bed(x[-(1:3)],
                                 sprintf("%s_%s_%s.csv", x$letter[1], x$stringency[1], x$confirmed[1])))

## alternative : I could apply below solution for list of confirmed peaks

res <- lapply(names(.Confired.ERs), function(nm) {
  mapply(export.bed,
         .Confirmed.ERs[[nm]],
         paste0(nm, "_", names(.res_accepted[[nm]]), ".bed"))
})

##=======================================================================================================

