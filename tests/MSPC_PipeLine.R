## MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
##
## updated workflow of MSPC Packages
##' @author Julaiti Shayiding

##--------------------------------------------------------------------------------------------------
#' @description pvalueConversion
.pvalueConversion <- function(x, pvalueBase = 10L, ...) {
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

##--------------------------------------------------------------------------------------------------

readPeakFiles <- function(peakFolder, pvalueBase=1L) {
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
        .gr <- .pvalueConversion(.gr, pvalueBase)
      }
    }), tools::file_path_sans_ext(basename(files))
  )
  return(f.read)
}

##' @examples
myData <- readPeakFiles(peakFolder = "test/dat_Musera/", pvalueBase = 10L)
options(scipen = 0)
##---------------------------------------------------------------------------------------------------

##' @details splitting GRanges
##' function implementation

.denoise.ERs <- function(grs, tau.w= 1.0E-04, .fileName="", outDir=getwd(), verbose=FALSE, ...) {
  # check input param
  stopifnot(inherits(grs[[1L]], "GRanges"))
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
    export.bed(.drop, sprintf("%s/%s.%s.bed", outDir, .grNM,.fileName))
    .keep <- .gr[.gr$p.value <= tau.w]
    return(.keep)
  })
  rslt <- setNames(res, names(grs))
  return(rslt)
}

#' @example
total.ERs <- .denoise.ERs(myData, tau.w = 1.0E-04, .fileName = "noisePeak", outDir = "test/")
options(scipen = 0)
##----------------------------------------------------------------------------------------------------

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
      out.idx0 <- out.idx0[!is.na(out.idx0)]
      ans <- ans[out.idx0]
    })
    res[[i]] <- DataFrame(c(list(que=queHit), sup=supHit))
    names(res[[i]]) <- c(names(peakset[i]),names(peakset[- i]))
  }
  rslt <- lapply(res, function(x) as.matrix(x[names(res[[1L]])]))
  rslt <- DataFrame(rbind(rslt[[1L]],
                          unique(do.call("rbind", rslt[2: length(rslt)]))))
  rslt <- lapply(rslt, function(x) as(x, "CompressedIntegerList"))
  return(rslt)
}

#' @example
#'
Hit <- .peakOverlapping(peakset = total.ERs, FUN = which.max)

##---------------------------------------------------------------------------------------------------

filterByOverlapHit <- function(.ovHit, peakset, replicate.type=c("Biological", "Technical"),
                   isSuffOverlap=TRUE, verbose=FALSE, ...) {
  # check input param
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  len <- length(peakset)-1,
                  len <- length(peakset))
  cnt.ovHit <- as.matrix(Reduce('+', lapply(.ovHit, lengths)))
  if(isSuffOverlap) {
    message("sufficient overlapped peaks are detected")
    keepHit <- lapply(.ovHit, function(ele_) {
      keepMe <- sapply(cnt.ovHit, function(x) x >= min.c)
      res <- ele_[keepMe]
    })
    return(keepHit)
  } else {
    message("peaks are discarded due to insufficient overlap")
    dropHit <- lapply(.ovHit, function(ele_) {
      droped <- sapply(cnt.ovHit, function(x) x < min.c)
      res <- ele_[droped]
    })
    rslt <- Map(unlist,
                mapply(extractList, peakset, dropHit))
    return(rslt)
  }
}

#' @example
keepList <- filterByOverlapHit(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=TRUE)
initDisc.ERs <- filterByOverlapHit(Hit, peakset = total.ERs, replicate.type = "Biological", isSuffOverlap=FALSE)

##----------------------------------------------------------------------------------------------------

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
  message("list of Fisher' combined p-value")
  comb.pval <- as.matrix(comb.pval)
  return(comb.pval)
}

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

#' @example
comb.p <- .Fisher.stats(keepList, total.ERs)

##-----------------------------------------------------------------------------------------------------

.filterByFisherMethod <- function(peakset, .ovHit, tau.s, comb.p ,isFisherPass=TRUE, ...) {
  # check input param
  stopifnot(class(peakset[[1L]])=="GRanges")
  stopifnot(is.numeric(tau.s))
  if(missing(comb.p)) {
    comb.p <- NA
    warning(paste("Fisher' combined pvalue", comb.p, "is missing"))
  }
  if(isFisherPass) {
    message("peaks are saved due to valid Fisher' combined p-value")
    lapply.func <- function(ele_) {
      keepMe <- sapply(comb.p, function(x) x<= tau.s)
      res <- ele_[keepMe]
    }
  } else {
    message("peaks are discarded due to obtaied Fisher combined p-value is not significant")
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
.Confirmed.ERs <- .filterByFisherMethod(total.ERs, keepList, tau.s=1.0E-08, comb.p, isFisherPass = TRUE)
.FisherDisc.ERs <- .filterByFisherMethod(total.ERs, keepList, tau.s=1.0E-08, comb.p, isFisherPass = FALSE)

##----------------------------------------------------------------------------------------------------
##' @description need to merge initDiscard.ERs and .FisherDiscard.ERs

.Discarded.ERs <- suppressWarnings(mapply(c, initDisc.ERs, .FisherDisc.ERs))

##----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------

.FDR.stats <- function(peakset_1, peakset_2, pAdjustMethod="BH", .fdr = 0.05
                       , replicate.type=c("Biological", "Technical"), outDir=getwd(), ...) {
  # check input param
  require(dplyr)
  if(missing(peakset_1) | missing(peakset_2)) {
    stop("required argument is missing")
  }
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

options(scipen = 0)
fdr.rslt <- .FDR.stats(.Confirmed.ERs, .Discarded.ERs, pAdjustMethod = "BH",
                       .fdr = 0.05, replicate.type = "Biological", outDir = "test/")
##===============================================================================================================

options(scipen = 0)
confirmedERs <- lapply(.Confirmed.ERs, as.data.frame)
discardedERs <- lapply(.Discarded.ERs, as.data.frame)

library(tidyverse)
library(magrittr)

create_output <- function(output_path, list1, list2, tau.s) {
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path))
  }

  #mList_1 <- lapply(list1, as.data.frame)
  #mList_2 <- lapply(list2, as.data.frame)

  names(list1) <- paste("confirmed", names(list1), sep = ".")
  names(list2) <- paste("discarded", names(list2), sep = ".")
  both <- do.call(rbind, c(list1, list2))
  both %<>% rownames_to_column(var = "cn")
  both %<>% separate(cn, c("original_list", "letters", "seq"), sep = "\\.")
  both %<>% mutate(peakStringency = ifelse(p.value <= tau.s , "Stringent", "Weak"))

  list_of_dfs <- both %>% split(list(.$letters, .$peakStringency, .$original_list))
  csv_names <- paste0(output_path, names(list_of_dfs), ".bed")
  return(mapply(export.bed, list_of_dfs, csv_names))
}

#' @example
options(scipen = 0)
.Confirmed.ERs
.Discarded.ERs

options(scipen = 0)
create_output(list1 = confirmedERs, list2 = discardedERs, tau.s = 1.0E-08, output_path ="test/")

##-----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------
## wrapper function for getting stringent / weak peaks accordingly
## below R script can generate expected format BED files in desired way

func.5 <- function(tau.s, ...) {
  # prepare data to be exported as BED file
  conf <- lapply(.Confirmed.ERs, as.data.frame)
  disc <- lapply(.Discarded.ERs, as.data.frame)

  both <- do.call("rbind", c(conf, disc))
  cn <- c("letter", "seq")
  DF <- cbind(read.table(text = chartr("_", ".", rownames(both)), sep = ".",
                         col.names = cn), both)
  DF <- transform(DF, stringency = ifelse(p.value <= tau.s, "Stringent", "Weak"))

  res <- by(DF, DF[c("letter", "stringency",".confirmed")],
            function(x) export.bed(x[-(1:length(.Confirmed.ERs))],
                                   sprintf("%s_%s_%s.bed", x$letter[1], x$stringency[1], x$confirmed[1])))
  return(res)
}

## Alternative solution for exporting nested list as BED file
res <- lapply(names(.res_accepted), function(nm)
  mapply(write.csv,
         .res_accepted[[nm]],
         paste0(nm, "_", names(.res_accepted[[nm]]), ".csv")))

##------------------------------------------------------------------------------------------------------
## Generating grouped bar plot, pie chart for exported bed files
#' @example mini data

savedDF <- list(
  bar.saved = data.frame(start=sample(100, 15), stop=sample(150, 15), score=sample(36, 15)),
  cat.saved = data.frame(start=sample(100, 20), stop=sample(100,20), score=sample(45,20)),
  foo.saved = data.frame(start=sample(125, 24), stop=sample(140, 24), score=sample(32, 24))
)

dropedDF <- list(
  bar.droped = data.frame(start=sample(60, 12), stop=sample(90,12), score=sample(35,12)),
  cat.droped = data.frame(start=sample(75, 18), stop=sample(84,18), score=sample(28,18)),
  foo.droped = data.frame(start=sample(54, 14), stop=sample(72,14), score=sample(25,14))
)

comb <- do.call("rbind", c(savedDF, dropedDF))
cn <- c("letter", "saved","seq")
DF <- cbind(read.table(text = chartr("_", ".", rownames(comb)), sep = ".", col.names = cn), comb)
DF <- transform(DF, updown = ifelse(score>= 12, "stringent", "weak"))
by(DF, DF[c("letter", "saved", "updown")],
   function(x) write.csv(x[-(1:3)],
                         sprintf("%s_%s_%s.csv", x$letter[1], x$updown[1], x$saved[1])))

##-----------------------------------------------------------------------------
## first solution for getting grouped bar plot
library(dplyr)
library(ggplot2)

plot_data <- DF %>%
  group_by(letter, saved, updown) %>%
  tally

ggplot(plot_data, aes(x = saved, y = n, fill = saved)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ letter + updown, ncol = 2)

##------------------------------------------------------------------------------
## or alternative solution for getting bar plot:
ggplot(plot_data, aes(x = letter, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(~updown+saved, ncol = 2)

##------------------------------------------------------------------------------
## solution for getting pie chart for file bar

ggplot(plot_data, aes(x = 1, y = percentage, fill = letter)) +
  geom_bar(stat = "identity", width =1) +
  facet_wrap(~updown+saved, ncol = 2) +
  coord_polar(theta = "y") +
  theme_void()

#--------------------------------------------------------------------------------
## alternative solution for getting pie char in different angle
library(dplyr)
library(tidyr)
library(ggplot2)

plot_data <- DF %>%
  unite(interaction, saved, updown, sep = "-") %>%
  group_by(letter, interaction) %>%
  tally %>%
  mutate(percentage = n/sum(n)) %>%
  filter(letter == "bar")

ggplot(plot_data, aes(x = 1, y = percentage, fill = interaction)) +
  geom_bar(stat = "identity", width =1) +
  coord_polar(theta = "y") +
  theme_void()

#---------------------------------------------------------------------------------
