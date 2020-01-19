#' @rdname vfloaders
#' @title Loaders from perimeters
#' @description Functions to load from commercial perimeters
#' @details The XML loader for the Humphrery Field Analyser (HFA) by Carl Zeiss Meditec
#' is essentially a XML parser that reads in the XML generated with the scientific
#' export license. The DICOMM loader is also a parser to read HFA data generated in a
#' DICOMM file. The loader for the Octopus perimeter by Haag-Streit is a csv reader
#' from files generated with the Eyesuite software
#' @return visual field data
#' @param file file to load (it is a XML extracted with the export license for
#' the HFA loader), a CSV outputed by the Eyesuite software for the Octopus
#' perimeter and a CSV with two columns for the batch HFA loader. The two columns
#' for the batch HFA loader must be named `\code{file}` and `\code{type}` and must
#' have the full file name (file path + name) of each XML file to be loaded and the
#' corresponding patient type, respectively
#' @param type type of patient. It can be `\code{ctr}` (for control or healthy
#' subject-eye) or `\code{pwg}` (for patient with glaucoma) or other
#' @param repeated function to apply if there are repeated values in a particular location
#' @export
loadhfaxml <- function(file, type = "pwg", repeated = mean) {
  vf <- td <- tdp <- pd <- pdp <- g <- gp <- NA
  dat <- xmlToList(xmlParse(file))$PATIENT
  id <- dat$PATIENT_ID
  dob <- as.Date(dat$BIRTH_DATE)
  dat <- dat$STUDY # next level
  date <- as.Date(dat$VISIT_DATE)
  age <- getage(dob, date)
  dat <- dat$SERIES # next level
  eye <- switch(dat$SITE, "0" = "OS", "1" = "OD")
  dat <- dat$FIELD_EXAM # next level
  time <- dat$EXAM_TIME
  # get the location map
  dat <- dat$STATIC_TEST # next level
  duration <- dat$EXAM_DURATION
  if(dat$FALSE_POSITIVE_METHOD == "1") {
    fpr <- as.numeric(dat$FALSE_POSITIVE_PERCENT) / 100
  } else if(dat$FALSE_POSITIVE_METHOD == "0") {
    fpr <- as.numeric(dat$FALSE_POSITIVES$ERRORS) / as.numeric(dat$FALSE_POSITIVES$TRIALS)
  } else {
    fpr <- as.numeric(NA)
  }
  if(dat$FALSE_NEGATIVE_METHOD == "1") {
    fnr <- as.numeric(dat$FALSE_NEGATIVE_PERCENT) / 100
  } else if(dat$FALSE_NEGATIVE_METHOD == "0") {
    fnr <- as.numeric(dat$FALSE_NEGATIVES$ERRORS) / as.numeric(dat$FALSE_NEGATIVES$TRIALS)
  } else {
    fnr <- as.numeric(NA)
  }
  fl <- as.numeric(dat$FIXATION_CHECK$ERRORS) / as.numeric(dat$FIXATION_CHECK$TRIALS)
  info <- data.frame(id = id, eye = eye, date = date, time = time, age = age,
                     type = type, fpr = fpr, fnr = fnr, fl = fl, duration = duration,
                     stringsAsFactors = FALSE)
  dat <- dat$THRESHOLD_TEST # next level
  nloc <- as.numeric(dat$NUM_THRESHOLD_POINTS)
  # return coordinates and value and ensure they are sort as they should
  s <- lapply(1:nloc, function(i) {
    x <- as.numeric(dat$THRESHOLD_SITE_LIST[[i]]$X)
    y <- as.numeric(dat$THRESHOLD_SITE_LIST[[i]]$Y)
    idx <- grep("THRESHOLD", names(dat$THRESHOLD_SITE_LIST[[i]]))
    val <- do.call(repeated, list(x = as.numeric(dat$THRESHOLD_SITE_LIST[[i]][idx])))
    return(list(x, y, val))
  })
  s <- matrix(unlist(s), nloc, 3, byrow = TRUE)
  # if eye = "OS", x is negative
  if(eye == "OS") s[,1] <- -s[,1]
  s <- s[order(s[,1]),]
  s <- s[order(-s[,2]),]
  vf <- cbind(info, t(s[,3]))
  names(vf)[getvfcols()] <- paste0("l", 1:nloc)
  # need to find blindspots locations
  dat <- dat$STATPAC # next level
  nlocd <- as.numeric(dat$NUM_TOTAL_DEV_VALUE_POINTS)
  locd <- lapply(1:nlocd, function(i) {
    x <- as.numeric(dat$TOTAL_DEVIATION_VALUE_LIST[[i]]$X)
    y <- as.numeric(dat$TOTAL_DEVIATION_VALUE_LIST[[i]]$Y)
    return(list(x, y))
  })
  locd <- matrix(unlist(locd), nlocd, 2, byrow = TRUE)
  # if eye = "OS", x is negative
  if(eye == "OS") locd[,1] <- -locd[,1]
  locd <- locd[order(locd[,1]),]
  locd <- locd[order(-locd[,2]),]
  bs <- which(!(paste(s[,1], s[,2], sep = "_") %in% paste(locd[,1], locd[,2], sep = "_")))
  cutoffs <- c(50, 5, 2, 1, 0.5) # cutoff lookup for TD and PD probability levels
  # total deviation values
  if("NUM_TOTAL_DEV_VALUE_POINTS" %in% names(dat)) {
    td <- xmldevvals(dat$TOTAL_DEVIATION_VALUE_LIST, as.numeric(dat$NUM_TOTAL_DEV_VALUE_POINTS),
                     eye, bs)
    td <- cbind(info, t(td))
    names(td)[getvfcols()] <- paste0("l", 1:nloc)
  }
  # total deviation probability values
  if("NUM_TOTAL_DEV_PROBABILITY_POINTS" %in% names(dat)) {
    tdp <- cutoffs[1 + xmldevvals(dat$TOTAL_DEVIATION_PROBABILITY_LIST,
                                  as.numeric(dat$NUM_TOTAL_DEV_PROBABILITY_POINTS), eye, bs)]
    tdp <- cbind(info, t(tdp))
    names(tdp)[getvfcols()] <- paste0("l", 1:nloc)
  }
  # pattern deviation values
  if("NUM_PATTERN_DEV_VALUE_POINTS" %in% names(dat)) {
    pd <- xmldevvals(dat$PATTERN_DEVIATION_VALUE_LIST, as.numeric(dat$NUM_PATTERN_DEV_VALUE_POINTS),
                     eye, bs)
    pd <- cbind(info, t(pd))
    names(pd)[getvfcols()] <- paste0("l", 1:nloc)
  }
  # pattern deviation probability values
  if("NUM_PATTERN_DEV_PROBABILITY_POINTS" %in% names(dat)) {
    pdp <- cutoffs[1 + xmldevvals(dat$PATTERN_DEVIATION_PROBABILITY_LIST,
                                  as.numeric(dat$NUM_PATTERN_DEV_PROBABILITY_POINTS),
                                  eye, bs)]
    pdp <- cbind(info, t(pdp))
    names(pdp)[getvfcols()] <- paste0("l", 1:nloc)
  }
  cutoffs <- c(50, NA, 5, 2, 1, 0.5) # cutoff lookup for probability levels of global indices
  if("GLOBAL_INDICES" %in% names(dat)) {
    vals <- vf[,getvfcols()]
    if(length(bs) > 1) vals <- vals[,-bs]
    g <- data.frame(msens = apply(vals, 1, mean),
                    ssens = apply(vals, 1, sd),
                    tmd = as.numeric(dat$GLOBAL_INDICES$MD),
                    tsd = NA,
                    pmd = NA,
                    psd = as.numeric(dat$GLOBAL_INDICES$PSD),
                    vfi = as.numeric(dat$GLOBAL_INDICES$VFI))
    gp <- data.frame(msens = NA,
                     ssens = NA,
                     tmd = cutoffs[1 + as.numeric(dat$GLOBAL_INDICES$MD_PROBABILITY)],
                     tsd = NA,
                     pmd = NA,
                     psd = cutoffs[1 + as.numeric(dat$GLOBAL_INDICES$PSD_PROBABILITY)],
                     vfi = NA)
    g  <- cbind(info, g)
    gp <- cbind(info, gp)
  }
  return(list(vf = vf, td = td, tdp = tdp, pd = pd, pdp = pdp, g = g, gp = gp))
}

#' @rdname vfloaders
#' @export
loadhfadicom <- function(file, type = "pwg", repeated = mean) {
  vf <- td <- tdp <- pd <- pdp <- g <- gp <- NA
  # load and arrange data for processing
  dat <- readDICOM(file)$hdr
  dat <- as.data.frame(eval(parse(text = paste0("dat$`", file, "`"))),
                       stringsAsFactors = FALSE)
  # extract the groups we are interested on
  groups <- sort(unique(dat$group))
  # get test grid
  gridtxt <- grep("Test Pattern", dicomelement(dat, groups[2], "CodeMeaning"), value = TRUE)
  # depending on the grid, the bs is in different locations
  if(grepl("24-2", gridtxt)) bs <- visualFields::locmaps$p24d2$bs
  if(grepl("30-2", gridtxt)) bs <- visualFields::locmaps$p30d2$bs
  if(grepl("10-2", gridtxt)) bs <- visualFields::locmaps$p10d2$bs
  if(grepl("30-1", gridtxt)) bs <- visualFields::locmaps$p30d1$bs
  if(grepl("60-4", gridtxt)) bs <- visualFields::locmaps$p60d4$bs
  # get id, test date, eye, age, time, and test duration
  id   <- dicomelement(dat, groups[3], "PatientID")
  dob  <- as.Date(dicomelement(dat, groups[3], "PatientsBirthDate"), "%Y%m%d")
  date <- as.Date(dicomelement(dat, groups[8], "PerformedProcedureStepStartDate"), "%Y%m%d")
  age <- getage(dob, date)
  eye <- switch(dicomelement(dat, groups[5], "Laterality"), "L" = "OS", "R" = "OD")
  time <- round(as.numeric(dicomelement(dat, groups[8], "PerformedProcedureStepStartTime")))
  time <- substr(gsub('(?=(?:.{2})+$)', ":", time, perl = TRUE), 2, 9)
  duration <- formatDuration(as.numeric(dicomelement(dat, groups[7], "VisualFieldTestDuration")))
  # get false positives, false negatives, and fixation losses
  if(dicomelement(dat, groups[7], "FalsePositivesEstimateFlag") == "YES") {
    fpr <- as.numeric(dicomelement(dat, groups[7], "FalsePositives")) /
           as.numeric(dicomelement(dat, groups[7], "PositiveCatchTrials"))
  } else {
    fpr <- as.numeric(NA)
  }
  if(dicomelement(dat, groups[7], "FalseNegativesEstimateFlag") == "YES") {
    fnr <- as.numeric(dicomelement(dat, groups[7], "FalseNegatives")) /
           as.numeric(dicomelement(dat, groups[7], "NegativeCatchTrials"))
  } else {
    fnr <- as.numeric(NA)
  }
  fl <- as.numeric(dicomelement(dat, groups[7], "PatientNotProperlyFixatedQuantity")) /
        as.numeric(dicomelement(dat, groups[7], "FixationCheckedQuantity"))
  info <- data.frame(id = id, eye = eye, date = date, time = time, age = age,
                     type = type, fpr = fpr, fnr = fnr, fl = fl, duration = duration,
                     stringsAsFactors = FALSE)
  # get sensitivity, TD, PD and probability values
  s <- data.frame(
    x    = as.numeric(dicomelement(dat, groups[7], "VisualFieldTestPointXCoordinate")),
    y    = as.numeric(dicomelement(dat, groups[7], "VisualFieldTestPointYCoordinate")),
    val  = as.numeric(dicomelement(dat, groups[7], "SensitivityValue")),
    td   = as.numeric(dicomelement(dat, groups[7], "AgeCorrectedSensitivityDeviationValue")),
    tdp  = as.numeric(dicomelement(dat, groups[7], "AgeCorrectedSensitivityDeviationProbabilityValue")),
    pd   = as.numeric(dicomelement(dat, groups[7], "GeneralizedDefectCorrectedSensitivityDeviationValue")),
    pdp  = as.numeric(dicomelement(dat, groups[7], "GeneralizedDefectCorrectedSensitivityDeviationProbabilityValue")),
    seen = dicomelement(dat, groups[7], "StimulusResults")
  )
  if(length(bs) > 1) s$td[bs] <- s$tdp[bs] <- s$pd[bs] <- s$pdp[bs] <- NA
  s$val[s$seen == "NOT SEEN"] <- -2
  s$seen <- NULL
  if(eye == "OS") s$x <- -s$y
  s <- s[order(s$x),]
  s <- s[order(-s$y),]
  vf  <- cbind(info, t(s$val))
  td  <- cbind(info, t(s$td))
  tdp <- cbind(info, t(s$tdp))
  pd  <- cbind(info, t(s$pd))
  pdp <- cbind(info, t(s$pdp))
  names(vf)[getvfcols()] <- names(td)[getvfcols()] <- names(tdp)[getvfcols()] <- 
    names(pd)[getvfcols()] <- names(pdp)[getvfcols()] <- paste0("l", 1:nrow(s))
  cutoffs <- c(50, NA, 5, 2, 1, 0.5) # cutoff lookup for probability levels of global indices
  if("GLOBAL_INDICES" %in% names(dat)) {
    vals <- vf[,getvfcols()]
    if(length(bs) > 1) vals <- vals[,-bs]
    g <- data.frame(msens = apply(vals, 1, mean),
                    ssens = apply(vals, 1, sd),
                    tmd = as.numeric(dicomelement(dat, groups[7], "GlobalDeviationFromNormal")),
                    tsd = NA,
                    pmd = NA,
                    psd = as.numeric(dicomelement(dat, groups[7], "LocalizedDeviationfromNormal")),
                    vfi = NA)
    gp <- data.frame(msens = NA,
                     ssens = NA,
                     tmd = as.numeric(dicomelement(dat, groups[7], "GlobalDeviationProbability")),
                     tsd = NA,
                     pmd = NA,
                     psd = as.numeric(dicomelement(dat, groups[7], "LocalizedDeviationProbability")),
                     vfi = NA)
    g  <- cbind(info, g)
    gp <- cbind(info, gp)
  }
  return(list(vf = vf, td = td, tdp = tdp, pd = pd, pdp = pdp, g = g, gp = gp))
}

#' @rdname vfloaders
#' @param dateFormat format to be used for date. Its default value
#' is \code{\%d.\%m.\%Y}
#' @export
loadoctopus <- function(file, type = "pwg", repeated = mean, dateFormat = "%d.%m.%Y") {
}

#' @rdname vfloaders
#' @export
loadhfaxmlbatch <- function(file, repeated = mean) {
  csvforxml <- read.csv(file, stringsAsFactors = FALSE)
  vflist <- loadhfaxml(csvforxml[1,1], csvforxml[1,2], repeated = repeated)
  for(i in 2:nrow(csvforxml)) {
    vflist1 <- loadhfaxml(csvforxml[i,1], csvforxml[i,2], repeated = repeated)
    vflist$vf <- vfjoin(vflist$vf, vflist1$vf)
    # join only whatever else is available
    if(is.data.frame(vflist1$td))  vflist$td  <- vfjoin(vflist$td,  vflist1$td)
    if(is.data.frame(vflist1$tdp)) vflist$tdp <- vfjoin(vflist$tdp, vflist1$tdp)
    if(is.data.frame(vflist1$pd))  vflist$pd  <- vfjoin(vflist$pd,  vflist1$pd)
    if(is.data.frame(vflist1$pdp)) vflist$pdp <- vfjoin(vflist$pdp, vflist1$pdp)
    if(is.data.frame(vflist1$g))   vflist$g   <- vfjoin(vflist$g,   vflist1$g)
    if(is.data.frame(vflist1$gp))  vflist$gp  <- vfjoin(vflist$gp,  vflist1$gp)
  }
  return(vflist)
}

#' @rdname vfloaders
#' @export
loadhfadicombatch <- function(file, repeated = mean) {
  csvfordicom <- read.csv(file, stringsAsFactors = FALSE)
  vflist <- loadhfadicom(csvfordicom[1,1], csvfordicom[1,2], repeated = repeated)
  for(i in 2:nrow(csvfordicom)) {
    vflist1 <- loadhfaxml(csvfordicom[i,1], csvfordicom[i,2], repeated = repeated)
    vflist$vf <- vfjoin(vflist$vf, vflist1$vf)
    # join only whatever else is available
    if(is.data.frame(vflist1$td))  vflist$td  <- vfjoin(vflist$td,  vflist1$td)
    if(is.data.frame(vflist1$tdp)) vflist$tdp <- vfjoin(vflist$tdp, vflist1$tdp)
    if(is.data.frame(vflist1$pd))  vflist$pd  <- vfjoin(vflist$pd,  vflist1$pd)
    if(is.data.frame(vflist1$pdp)) vflist$pdp <- vfjoin(vflist$pdp, vflist1$pdp)
    if(is.data.frame(vflist1$g))   vflist$g   <- vfjoin(vflist$g,   vflist1$g)
    if(is.data.frame(vflist1$gp))  vflist$gp  <- vfjoin(vflist$gp,  vflist1$gp)
  }
  return(vflist)
}

###################################################################################
# INTERNAL FUNCTIONS: routines to ease load on XML loader
###################################################################################
# helps parse XML TD and PD values and their corresponding probability levels
#' @noRd
xmldevvals <- function(dat, nloc, eye, bs) {
  # return coordinates and value
  res <- lapply(1:nloc, function(i) {
    x <- as.numeric(dat[[i]][1])
    y <- as.numeric(dat[[i]][2])
    val <- as.numeric(dat[[i]][3])
    return(list(x, y, val))
  })
  res <- matrix(unlist(res), nloc, 3, byrow = TRUE)
  # if eye = "OS", x is negative
  if(eye == "OS") res[,1] <- -res[,1]
  res <- res[order(res[,1]),]
  res <- res[order(-res[,2]),]
  res <- res[,3]
  if(length(bs) > 0) {
    idx <- 1:(nloc + length(bs))
    aux <- rep(NA, nloc + length(bs))
    aux[idx[-bs]] <- res
    res <- aux
  }
  return(res)
}

# get value from an element of a group in a DICOM scheme
#' @noRd
dicomelement <- function(dat, group, name) {
  datgr <- dat[which(dat$group == group),]
  return(datgr$value[datgr$name == name])
}

# format duration from seconds to hh:mm:ss
#' @noRd
formatDuration <- function(ss) {
  # convert from seconds to hh:mm:ss
  hh <- floor(ss / 360)
  ss <- ss - 360 * hh
  mm <- floor(ss / 60)
  ss <- ss - 60 * mm
  # convert to strings
  hh <- as.character(hh)
  mm <- as.character(mm)
  ss <- as.character(ss)
  if(nchar(hh) == 1) hh <- paste0("0", hh)
  if(nchar(mm) == 1) mm <- paste0("0", mm)
  if(nchar(ss) == 1) ss <- paste0("0", ss)
  return(paste(hh, mm, ss, sep = ":"))
}