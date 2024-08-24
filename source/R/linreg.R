#' @rdname linreg
#' @title Global and pointwise linear regression analyses
#' @description Functions that compute global and pointwise linear regression analyses:
#' \itemize{
#'   \item\code{glr} performs global linear regression analysis
#'   \item\code{plr} performs pointwise linear regression (PLR) analysis
#'   \item\code{poplr} performs PoPLR analysis as in O'Leary et al (see reference)
#' }
#' @details
#' \itemize{
#'   \item\code{poplr} there is a small difference between this implementation of
#'     PoPLR and that proposed by O'Leary et al. The combined S statistic in the
#'     paper used a natural logarithm. Here we not only use a logarithm of base 10
#'     but we also divide by the number of locations. This way the S statistic has
#'     a more direct interpretation as the average number of leading zeros in the
#'     p-values for pointwise (simple) linear regression. That is, if S = 2, then
#'     the p-values have on average 2 leading zeros, if S = 3, then 3 leading zeros,
#'     and so on
#' }
#' @return
#' \itemize{
#'   \item{\code{glr} and \code{plr} return a list with the following
#'     \itemize{
#'       \item\code{id} patient ID
#'       \item\code{eye} patient eye
#'       \item\code{testSlope} slope for \code{glr} or list of slopes for \code{plr}
#'         to test as null hypotheses
#'       \item\code{nvisits} number of visits
#'       \item\code{years} years from baseline. Used for the pointwise linear
#'         regression analysis
#'       \item\code{data} data analyzed. For \code{glr}, it is the values of the
#'         global indes analyzed. For \code{plr}, each column is a location of the
#'         visual field used for the analysis. Each row is a visit (as many as years)
#'       \item\code{pred} predicted values. Each column is a location of the visual
#'         field used for the analysis. Each row is a visit (as many as years)
#'       \item\code{sl} slopes estimated at each location for pointwise (simple)
#'         linear regression
#'       \item\code{int} intercept estimated at each location for pointwise (simple)
#'         linear regression
#'       \item\code{tval} t-values obtained for the left-tailed-t-tests for the slopes
#'         obtained in the pointwise (simple) linear regression at each location
#'       \item\code{pval} p-values obtained for the left-tailed t-tests for the slopes
#'         obtained
#'     }
#'   }
#'   \item{\code{poplr} returns a list with the following additional fields
#'     \itemize{
#'       \item\code{csl} the modified Fisher's S-statistic for the left-tailed permutation test
#'       \item\code{cslp} the p-value for the left-tailed permutation test
#'       \item\code{csr} the modifed Fisher's S-statistic for the right-tailed permutation test
#'       \item\code{csrp} the p-value for the right-tailed permutation test
#'       \item\code{pstats} a list with the poinwise slopes (`\code{sl}`), intercepts
#'         (`\code{int}`), standard errors (`\code{se}`), and p-values (`\code{pval}`) obtained
#'         for the series at each location analyzed and for all \code{nperm} permutations
#'         (in `\code{permutations}`)
#'       \item\code{cstats} a list with all combined stats:
#'         \itemize{
#'           \item\code{csl, csr} the combined Fisher S-statistics for the left- and right-tailed
#'             permutation tests respectively
#'           \item\code{cslp, csrp} the corresponding p-values for the permutation tests
#'           \item\code{cslall, csrall} the combined Fisher S-statistics for all permutations
#'         }
#'     }
#'   }
#' }
#' @param g a data.frame with date on the first column and the value of the
#'        global index on the second column
#' @param testSlope slope, or slopes, to test as null hypothesis. Default is 0.
#'   if a single value, then the same null hypothesis is used for all locations.
#'   If a vector of values, then (for \code{plr} and \code{poplr}) each
#'   location of the visual field will have a different null hypothesis. The length of
#'   testSlope must be 1 or equal to the number of locations to be used in the PLR or
#'   PoPLR analysis
#' @references
#' N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in
#' glaucoma: estimating the overall significance of deterioration with permutation
#' analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology
#' and Visual Science, 53, 2012
#' @examples
#' vf <- vffilter(vfpwgRetest24d2, id == 1) # select one patient
#' res <- glr(getgl(vf)[,c("date", "tmd")]) # linear regression with mean deviation (MD)
#' res <- plr(gettd(vf))   # pointwise linear regression (PLR) with TD values
#' res <- poplr(gettd(vf)) # Permutation of PLR with TD values
#' @export
glr <- function(g, testSlope = 0) {
  years <- as.numeric(g[,1] - g[1,1]) / 365.25 # it should be difference in years from baseline date
  y <- g[,2]
  y <- y[order(years)] # sort just in case
  years <- years[order(years)]
  nvisits <- length(y)
  precision <- 1e-6
  X <- matrix(c(rep(1, length(years)), years), nvisits, 2)
  ssvf <- (nvisits - 1 ) * var(y)
  ssyears <- (nvisits - 1) * var(years)
  kvyears <- (nvisits - 2) * ssyears
  reg <- t(solve(t(X) %*% X) %*% t(X) %*% y)
  int <- reg[1]
  sl <- reg[2]
  v <- (ssvf - ssyears * sl^2) / kvyears
  v[v < 0] <- 0
  se <- sqrt(v)
  if(sd(y) <= precision) {
    int <- as.numeric(y[1])
    sl <- 0
    se <- precision
  }
  tval <- (sl - testSlope) / se
  pval <- pt(tval, nvisits - 2)
  pred <- sl * years + int
  return(list(id = g$id[1], eye = g$eye[1], testSlope = testSlope,
              nvisits = nvisits, dates = g$date, years = years, data = y, pred = pred,
              sl = sl, int = int, se = se, tval = tval, pval = pval))
}

#' @rdname linreg
#' @param vf visual fields sensitivity data
#' @export
plr <- function(vf, testSlope = 0) {
  if(nrow(unique(data.frame(vf$id, vf$eye))) != 1)
    stop("all visual fields must belong to the same subject id and eye")
  vf <- vfsort(vf) # sort just in case
  nvisits <- nrow(vf)
  years <- as.numeric(vf$date - vf$date[1]) / 365.25 # it should be difference in years from baseline date
  bs <- getlocmap()$bs
  y <- vf[,getvfcols()]
  y[,bs] <- NA # ignore blind spot locations in the analysis
  precision <- 1e-6
  X <- matrix(c(rep(1, length(years)), years), nvisits, 2)
  ssvf <- (nvisits - 1 ) * apply(y, 2, var)
  ssyears <- (nvisits - 1) * var(years)
  kvyears <- (nvisits - 2) * ssyears
  reg <- t(solve(t(X) %*% X) %*% t(X) %*% as.matrix(y))
  int <- t(reg[,1])
  sl <- t(reg[,2])
  v <- (ssvf - ssyears * sl^2) / kvyears
  v[v < 0] <- 0
  se <- sqrt(v)
  # get the locations for which sensitivity did not change
  invariantloc <- which(apply(y, 2, sd) <= precision)
  # locations with non-changing series in sensitivity: slope is zero,
  # intercept is not defined, and standard error is nominally very small
  int[invariantloc] <- as.numeric(y[1,invariantloc])
  sl[invariantloc] <- 0
  se[invariantloc] <- precision
  tval <- (sl - testSlope) / se
  pval <- pt(tval, nvisits - 2)
  # convert all to data frames and assign the corresponding column names, then return
  sl <- as.data.frame(sl)
  int <- as.data.frame(int)
  se <- as.data.frame(se)
  tval <- as.data.frame(tval)
  pval <- as.data.frame(pval)
  # predicted values
  pred <- sapply(as.list(rbind(int, sl)), function(beta) {beta[1] + beta[2] * years})
  return(list(id = vf$id[1], eye = vf$eye[1], testSlope = testSlope,
              nvisits = nvisits, dates = vf$date, years = years, data = vf[,getvfcols()], pred = pred,
              sl = sl, int = int, se = se, tval = tval, pval = pval))
}

#' @rdname linreg
#' @param nperm number of permutations. If the number of visits is 7 or less, then
#' \code{nperm = factorial(nrow(vf))}. For series greater than 8 visits, default is
#' factorial(7). For series up to 7 visits, it is the factorial of the number of visits
#' (with less than 7 visits, the number of possible permutations is small and results
#' can be unreliable. For instance, for 5 visits, the number of possible permutations is
#' only 120.)
#' @param trunc truncation value for the Truncated Product Method (see reference)
#' @export
poplr <- function(vf, testSlope = 0, nperm = factorial(7), trunc = 1) {
  if(nrow(unique(data.frame(vf$id, vf$eye))) != 1)
    stop("all visual fields must belong to the same subject id and eye")
  vf <- vfsort(vf) # sort just in case
  nvisits <- nrow(vf)
  years <- as.numeric(vf$date - vf$date[1]) / 365.25 # it should be difference in years from baseline date
  bs <- getlocmap()$bs
  y <- vf[,getvfcols()]
  y[,bs] <- NA # ignore blind spot locations in the analysis
  if(nvisits < 7) warning("random permutation analysis may be imprecise with less than 7 visual fields")
  if(nvisits < 8) {
    porder <- matrix(unlist(permn(nvisits)), ncol = nvisits, byrow = TRUE)
    # is number of permutations is smaller than nrow(porder) do random sampling
    if(nperm < nrow(porder))
      porder <- rbind(porder[1,], porder[sample(nrow(porder), nperm - 1),])
    else nperm <- nrow(porder)
  } else {
    if(nperm > 10000)
      stop("I'm sorry Dave, I'm afraid I can't do that. I think you know what the problem is just as well as I do.")
    porder <- t(replicate(factorial(8), c(1:nvisits)[sample(nvisits)]))
    porder <- rbind(c(1:nvisits), porder)
    porder <- unique(porder)[1:nperm,]
    if(nrow(porder) != nperm)
      stop("something went wrong and did not get the number of permutations you wanted")
  }
  # get the p-values from pointwise linear regression for series and all permutations
  pstats <- poplrpvals(y, years, porder, testSlope)
  # ... and compute the combined S statistic, after removing the blind spot
  pval <- pstats$permutations$pval
  if(length(bs) > 0) pval <- pval[,-bs]
  cstats <- poplrsstats(pval, trunc = trunc)
  # predicted values
  pred <- sapply(as.list(rbind(pstats$int, pstats$sl)), function(beta) {beta[1] + beta[2] * years})
  return(list(id = vf$id[1], eye = vf$eye[1], testSlope = testSlope,
              nvisits = nvisits, dates = vf$date, years = years, data = vf[,getvfcols()], pred = pred,
              sl = pstats$sl, int = pstats$int, se = pstats$se, tval = pstats$tval,
              pval = pstats$pval, nperm = nperm,
              csl = cstats$csl, cslp = cstats$cslp,
              csr = cstats$csr, csrp = cstats$csrp,
              pstats = pstats, cstats = cstats))
}

###################################################################################
# INTERNAL FUNCTIONS: routines to ease load on poplr code
###################################################################################
# Internal functions, computes p-values from simple linear regression in all locations (columns)
# in vf for the series and for all permutations in porder
#' @noRd
poplrpvals <- function(vf, years, porder, testSlope = 0) {
  if(!is.numeric(testSlope)) stop("testSlope must be numeric")
  if(!(length(testSlope) == 1 || length(testSlope) == ncol(vf)))
    stop("testSlope must be either a numeric scalar or a vector with the same length as columns there are in vf")
  colNames <- names(vf)
  vf <- as.matrix(vf)
  # number of permutations, locations, and tests
  nperm <- nrow(porder)
  nloc  <- ncol(vf)
  nvisits <- nrow(vf)
  precision <- 1e-6
  sl <- matrix(c( NA ), nrow = nperm, ncol = nloc)
  int <- matrix(c( NA ), nrow = nperm, ncol = nloc)
  se <- matrix(c( NA ), nrow = nperm, ncol = nloc)
  # add defaults for slope hypothesis tests when slr analysis is to be performed
  if(length(testSlope) == 1) testSlope <- rep(testSlope, nloc)
  # point-wise linear regression over time permutation-invarian values
  syears <- sum(years)
  myears <- mean(years)
  ssyears <- (nvisits - 1) * var(years)
  kvyears <- (nvisits - 2) * ssyears
  mvf <- apply(vf, 2, mean)
  ssvf <- (nvisits - 1 ) * apply(vf, 2, var)
  # compute slopes per location, ...
  sl <- sapply(1:nloc, function(loc)
    (matrix(years[porder], nrow(porder), ncol(porder)) %*% vf[,loc] - syears * mean(vf[,loc])) / ssyears)
  # ..., and then intercepts, ...
  int <- matrix(rep(mvf, nperm), nrow(sl), ncol(sl), byrow = TRUE) - myears * sl
  # ..., and then compute standard errors.
  varslope <- (matrix(rep(ssvf, nperm), nrow(sl), ncol(sl), byrow = TRUE) - ssyears * sl^2) / kvyears
  varslope[varslope < 0] <- 0
  se <- sqrt(varslope)
  # get the locations for which sensitivity did not change
  invariantloc <- which(apply(vf, 2, sd) <= precision)
  # locations with non-changing series in sensitivity: slope is zero,
  # intercept is not defined, and standard error is nominally very small
  sl[,invariantloc] <- 0
  int[,invariantloc] <- vf[1,invariantloc]
  se[,invariantloc] <- precision
  # Get t-values and the corresponding p-values
  tval <- (sl - t(matrix(rep(testSlope, nperm), nloc, nperm))) / se
  pval <- pt(tval, nvisits - 2)
  pval[pval < precision] <- precision
  pval[pval > (1 - precision)] <- 1 - precision
  # convert all to data frames and assign the corresponding column names, then return
  sl   <- as.data.frame(sl)
  int  <- as.data.frame(int)
  se   <- as.data.frame(se)
  tval <- as.data.frame(tval)
  pval <- as.data.frame(pval)
  names(sl) <- colNames
  names(int) <- colNames
  names(se) <- colNames
  names(tval) <- colNames
  names(pval) <- colNames
  return(list(sl = sl[1,], int = int[1,], se = se[1,], tval = tval[1,], pval = pval[1,],
              permutations = list(sl = sl, int = int, se = se, tval = tval, pval = pval)))
}

# Internal function: computes the modified Fisher S, applying the Truncated Product Method,
# if requested, from the p-values obtained for the series at each location and for all
# permutations. It returns the p-value based on the observed Fisher S statistic and the
# distribution obtained from the series permutations. It does so for a left-tailed hypothesis
# test and for the right-tailed hypothesis test
#' @noRd
poplrsstats <- function(pval, trunc = 1) {
  ##############
  # input checks
  ##############
  # truncation must be between zero and one
  if(trunc <= 0 | trunc > 1)
    stop("truncation must be between 0 and 1")
  # init
  nperm <- nrow(pval)
  # Apply the Truncated Product Method if required (i.e. trunc between 0 and 1)
  # left-tail analysis
  tpl <- apply(pval, 1, min)
  tpl[tpl < trunc] <- trunc
  kl <- matrix(rep(1, nrow(pval) * ncol(pval)), nrow(pval), ncol(pval))
  kl[pval > tpl] <- 0
  # right-tail analysis
  pvalr <- 1 - pval
  tpr <- apply(pvalr, 1, min)
  tpr[tpr < trunc] <- trunc
  kr <- matrix(rep(1, nrow(pvalr) * ncol(pvalr)), nrow(pvalr), ncol(pvalr))
  kr[pvalr > tpr] <- 0
  # combine p-value test statistics with a modified Fisher S statistic
  csl <- -rowSums(kl * log(pval))
  csr <- -rowSums(kr * log(1 - pval))
  # observed and permutation test statistics
  cslp <- 1 - rank(csl) / nperm
  csrp <- 1 - rank(csr) / nperm
  return(list(csl = csl[1], cslp = cslp[1], cslall = csl, cslpall = cslp,
              csr = csr[1], csrp = csrp[1], csrall = csr, csrpall = csrp))
}