#' @rdname getage
#' @title Calculates age
#' @description Computes ages at specific dates
#' @param dob date(s) of birth
#' @param date date(s) for which to calculate age
#' @return \code{getage} returns the age from the date of birth and a certain date
#' @examples
#' getage("1977-01-31", "2014-01-30")
#' @export
getage <- function(dob, date)  { 
  dob  <- as.POSIXlt(dob)
  date <- as.POSIXlt(date)
  age  <- date$year - dob$year
  # if month of DoB has not been reached yet, then a year younger
  idx <- which(date$mon < dob$mon)
  if(length(idx) > 0) age[idx]  <- age[idx] - 1
  # if same month as DoB but day has not been reached, then a year younger
  idx <- which(date$mon == dob$mon & date$mday < dob$mday)
  if(length(idx) > 0) age[idx]  <- age[idx] - 1
  return(age)
}

#' @rdname vfstats
#' @title Statistical analyses for visual fields data
#' @description 
#' \itemize{
#'   \item\code{vfaggregate} computes summary statistics of visual field data
#'   \item\code{vfmean} computes the mean statistics of visual field data. It is
#'     a wrapper for vfaggregate but only to compute means
#'   \item\code{vfretestdist} computes the conditional distribution from test-retest data
#' }
#' @details
#' \itemize{
#'   \item\code{vfaggregate} this is a restricted version of \code{\link{aggregate}}
#'     that only allows to use part of the key hierarchically, and operates on all
#'     data frames of the \code{VisualField} object. The restriction is that only
#'     aggregates that are allowed are `\code{newkey = c("id", "eye")}` and
#'     `\code{newkey = c("id", "eye", "date")}`. It returns the aggregated value for all
#'     numeric columns grouped and ordered by the new key (id and eye, or id, eye,
#'     and date). If the aggregate grouping is by \code{eye} and the function, then
#'     the \code{date} returned is the average.
#' }
#' @param by aggregate by \code{date}, that is by id, eye, and date (default) or by
#'   \code{eye}, that is by id and eye
#' @param fun a function to compute the summary statistics which can be applied to
#'   all data subsets. The default is `\code{mean}`
#' @return \code{vfaggregate} and \code{vfmean} return a vf data frame with aggregate values
#' @examples
#' # aggregate by date
#' vfaggregate(vfpwgRetest24d2, by = "date")           # compute the mean
#' vfaggregate(vfpwgRetest24d2, by = "date", fun = sd) # compute standard deviation
#' # aggregate by eye
#' vfaggregate(vfpwgRetest24d2, by = "eye")           # compute the mean
#' vfaggregate(vfpwgRetest24d2, by = "eye", fun = sd) # compute standard deviation
#' @export
vfaggregate <- function(vf, by = "date", fun = mean, ...) {
  if(by == "eye") {
    form <- . ~ id + eye
    nkey <- 1:2
  } else if(by == "date") {
    form <- . ~ id + eye + date
    nkey <- 1:3
  } else stop("type of summary statistics not allowed")
  # we need to do different things depending on what type of data we are dealing with
  # find which column of the non-key columns is numeric, which is character
  cnames <- names(vf)
  cclass <- sapply(cnames[5:length(cnames)], function(col) class(vf[,col]))
  numer <- which(cclass == "numeric" | cclass == "integer") + 4 # numeric columns
  other <- setdiff(5:length(cclass), numer) # all other columns
  # aggregate numeric columns
  vfa <- aggregate(form, data = vf[,c(nkey, numer)], fun, na.action = na.pass, ...)
  # time is no longer relevant and we set it to midnight
  vfa$time <- "00:00:00"
  # add column for date. if aggregate is for eye and fun is mean then find average
  if(by == "eye") {
    vfa$date <- as.Date(sapply(1:nrow(vfa), function(i) {
      return(as.character(mean.Date(vf$date[vf$id == vfa$id[i] & vf$eye == vfa$eye[i]])))
    }))
  }
  # sort
  vfa <- vfa[order(vfa$id, vfa$eye, vfa$date),]
  # add other columns. If they have the same value the new key, then keep
  # that value otherwise, blank
  for(col in cnames[other]) {
    vfa[,col] <- sapply(1:nrow(vfa), function(i) {
      if(by == "eye") idx <- which(vf$id == vfa$id[i] & vf$eye == vfa$eye[i])
      else idx <- which(vf$id == vfa$id[i] & vf$eye == vfa$eye[i] & vf$date == vfa$date[i]) # by == "date"
      if(length(unique(vf[idx,col])) == 1) return(vf[idx[1],col])
      return(NA)
    })
  }
  # return the aggregated visual field data with columns in the correct order
  return(vfa[,cnames])
}

#' @rdname vfstats
#' @examples
#' # mean by date
#' vfmean(vfpwgRetest24d2, by = "date")
#' # mean by eye
#' vfmean(vfpwgRetest24d2, by = "eye")
#' @export
vfmean <- function(vf, by = "date", ...) return(vfaggregate(vf, by, fun = mean, ...))

#' @rdname vfstats
#' @param vf a table with visual fields data. Data is rounded, which leaves
#' sensitivity data unchanged, but it is necessary for the nature of the
#' algorithm if the data passed are TD or PD values or summary stats such as
#' averages. Beware of the locations in the blind spot, which very likely need
#' to be removed
#' @param nbase number of visual fields to be used as baseline
#' @param nfollow number of visual fields to be used as follow up
#' @param alpha significance level to derive the conditional retest intervals.
#' Default value is \code{0.1}
#' @param ... arguments to be passed to or from methods. A useful one to try
#' is type of quantile calculation `\code{type}` use in \code{\link{quantile}}
#' @return \code{vfretestdist} returns a list with the following elements:
#' \itemize{
#'   \item\code{x} with all the test values (x-axis)
#'   \item\code{y} the distribution of retest dB values conditional to each
#'   test value in \code{x}. It is a list with as many entries as \code{x}
#'   \item\code{n} number of retest values conditional to each value in \code{x}.
#'   It is a list with as many entries as \code{x}
#'   \item\code{ymed} median for each value in \code{x}. It is a list with as
#'   many entries as \code{x}
#'   \item\code{ylow} quantile value for significance \code{1 - alpha / 2}
#'   for each value in \code{x}. It is a list with as many entries as \code{x}
#'   \item\code{yup} quantile value for significance \code{alpha / 2}
#'   for each value in \code{x}. It is a list with as many entries as \code{x}
#' }
#' Together \code{ylow} and \code{yup} represent the lower and upper limit of the
#' \code{(1 - alpha)\%} confidence intervals at each value \code{x}.
#' @examples
#' # get the retest sensitivity data after removing the blind spot
#' retest <- vfretestdist(vfpwgRetest24d2, nbase = 1, nfollow = 1)
#' 
#' plot(0, 0, typ = "n", xlim = c(0, 40), ylim = c(0,40),
#'      xlab = "test in dB", ylab = "retest in dB", asp = 1)
#' for(i in 1:length(retest$x)) {
#'   points(rep(retest$x[i], length(retest$y[[i]])), retest$y[[i]],
#'          pch = 20, col = "lightgray", cex = 0.75)
#' }
#' lines(c(0,40), c(0,40), col = "black")
#' lines(retest$x, retest$ymed, col = "red")
#' lines(retest$x, retest$ylow, col = "red", lty = 2)
#' lines(retest$x, retest$yup, col = "red", lty = 2)
#' @export
vfretestdist <- function(vf, nbase = 1, nfollow = 1, alpha = 0.1, ...) {
  # get all subject id and eye to process
  idall <- paste(vf$id, vf$eye, sep = "_")
  idu <- unique(idall)
  neyes <- length(idu)     # total number of eyes to process
  vf <- round(vf[,getvfcols()]) # remove the first 4 columns corresponding to the key
  if(length(getlocmap()$bs) > 0) vf <- vf[,-getlocmap()$bs] # remove blind spot if necessary
  # check if it can be done
  if(any(sapply(idu, function(idu) length(which(idall == idu))) < nbase + nfollow))
    stop("not enough baseline and follow up data for the retest analysis")
  # gather all possible values
  retest  <- NULL
  retest$x <- sort(unique(c(as.matrix(vf))))
  retest$y <- as.list(rep(NA, length(retest$x)))
  # construct all possible permutations from the series
  combs <- combinations((nbase + nfollow), nbase)
  # start analysis for each possible combination
  for(i in 1:nrow(combs)) {
    # do the analysis for each subject-eye
    for(j in 1:neyes) {
      # first get the baseline data vs the followup data
      idx <- which(idall == idu[j])
      vfiter <- vf[idx[1:(nbase + nfollow)],]
      dbbase <- vfiter[combs[i,],]
      dbfollow <- as.numeric(colMeans(vf[idx[-combs[i,]],]))
      # analyse baseline
      idxbase <- NULL
      for(k in 1:ncol(dbbase)) {
        if(all(dbbase[,k] == dbbase[1,k])) idxbase <- c(idxbase, k)
      }
      # mount conditional distribution
      if(length(idxbase) > 0) {
        for(k in 1:length(idxbase)) {
          idxcp <- which(retest$x == dbbase[1,idxbase[k]])
          if(length(retest$y[[idxcp]]) == 1 && is.na(retest$y[[idxcp]])) {
            retest$y[[idxcp]] <- dbfollow[idxbase[k]]
          } else {
            retest$y[[idxcp]] <- c(retest$y[[idxcp]], dbfollow[idxbase[k]])
          }
        }
      }
    }
  }
  # remove whatever has no data in it
  idx <- which(is.na(retest$y))
  if(length( idx ) > 0) {
    retest$x      <- retest$x[-idx]
    retest$y[idx] <- NULL
  }
  # sort the results and obtain percentiles at level alpha
  for(i in 1:length(retest$y)) {
    retest$y[[i]] <- sort(retest$y[[i]])
    qq <- as.numeric(quantile(retest$y[[i]],
                              probs = c(0.5, alpha / 2, 1 - alpha / 2), ...))
    retest$n[i] <- length(retest$y[[i]])
    retest$ymed[i] <- qq[1]
    retest$ylow[i] <- qq[2]
    retest$yup[i]  <- qq[3]
  }
  return(retest)
}