#' @rdname vf
#' @title Visual field dataset
#' @description The main object of the visualFields package is a table with 
#' a specific format and fields that are mandatory for their management and
#' processing (mainly statistical analysis). Each record (row) in the table
#' contains data for a single visual field test. The mandatory fields specify
#' subject (by its ID code), eye, and test date and time. There are required
#' fields statistical and reliability analyses (e.g., age for the determination
#' of total-deviation and pattern-deviation values, and for global indices and
#' fpr, fnr, fl for the proportion of false positives, false negative, and
#' fixation losses). The rest of mandatory fields are sensitivity or deviation
#' data for each visual field test location. (The number of fields for
#' tested locations varies with the location map, 54 for the 24-2, 76 for the
#' 30-2, 68 for the 10-2, etc.). Check section \code{Structure of visual fields data}
#' below for details about the required  structure of the table contatining the
#' visual fields datasets.
#' 
#' The following functions carry out analysis on visual fields data:
#' \itemize{
#'   \item\code{vfdesc} descriptive summary of a visual field dataset
#'   \item\code{vfsort} sort visual field data
#'   \item\code{vfisvalid} check if a table with visual field data is properly
#'     formatted and valid for analysis
#'   \item\code{vfread} read a csv file with visual field data
#'   \item\code{vfwrite} write a csv file with visual field data
#'   \item\code{vfjoin} joins two visual field datasets
#'   \item\code{vffilter} filters elements from a visual field dataset with
#'     matching conditions. This function is just a wrapper for \code{dplyr}'s
#'     function \code{\link{filter}}
#'   \item\code{vfselect} select visual field data by index or the first or last
#'     \code{n} visits per subject and eye
#'   \item\code{gettd} computes total-deviation (TD) values and probability values
#'   \item\code{gettdp} computes total-deviation (TD) probability values
#'   \item\code{getpd} computes pattern-deviation (PD) values
#'   \item\code{getpdp} computes pattern-deviation (PD) probability values
#'   \item\code{getgh} computes the general height (GH) from the TD tables
#'   \item\code{getgl} computes visual fields global indices
#'   \item\code{getglp} computes computes visual fields global indices probability values
#' }
#' @details
#' \itemize{
#'   \item\code{vfselect} when selecting the last or first few visual fields per
#'     subject and eye, if that subject and eye has fewer than \code{n} visits,
#'     then all visits are returned
#' }
#' @section Structure of visual fields data:
#' Visual fields data is the central object used in visualFields. It is a table of
#' visual field data collected with the same perimeter, background and stimulus
#' paradigm (e.g., static automated perimetry or frequency-doubling  perimetry),
#' stimulus size (e.g., Goldmann size III), grid of visual field test locations
#' (e.g., 24-2), and psychophysical testing strategy (e.g., SITA standard).
#' Normative values can be obtained from appropriate datasets with data for healthy
#' eyes and these normative values can then be used to generate statistical analyses
#' and visualizations of data for patients with retinal or visual anomalies.
#'
#' Each record correspond to a specific test for an eye of a subject taken on a
#' specific date at a specific time. Visual field data must have the following
#' columns
#' \itemize{
#'   \item\code{id} an id uniquely identifying a subject. This field is mandatory
#'   \item\code{eye} should be "OD" for right eye or "OS" for left eye. This field is
#'     mandatory
#'   \item\code{date} test date. This field is mandatory
#'   \item\code{time} test time. This field is mandatory
#'   \item\code{age} age of the patient on the test date. This field is required to
#'     obtain total-deviation, pattern-deviation values, and other age-dependent
#'     local and global indices
#'   \item\code{type} type of subject, Could be a healthy subject (ctr for control)
#'     or a patient with glaucoma (pwg) or a patient with idiopatic intraocular
#'     hypertension (iih) or other. This field is no required for management or
#'     statistical analysis.
#'   \item\code{fpr} false positive rate. This field is no required for management or
#'     statistical analysis.
#'   \item\code{fnr} false negative rate. This field is no required for management or
#'     statistical analysis.
#'   \item\code{fl} fixation losses. This field is no required for management or
#'     statistical analysis.
#'   \item\code{l1..ln} sensitivity, total-deviation, or pattern-deviation values for
#'     each location. For analysis with visualFields there should be as many columns
#'     as coordinates in the location map set in the visualFields environment. These
#'     fields are mandatory.
#' }
#' @examples
#' # get dataset description from visual field table
#' vfdesc(vfctrSunyiu24d2)
#' # sort dataset
#' vfsort(vfctrSunyiu24d2[c(5, 4, 10, 50, 30),])
#' # check if a visualField is valid
#' vf <- vfctrSunyiu24d2
#' vfisvalid(vf) # valid visual field data
#' vf$id[5] <- NA
#' vfisvalid(vf) # invalid visual field data
#' # write and read visual field data
#' vf <- vfctrSunyiu24d2
#' tf <- tempfile("vf")
#' vfwrite(vf, file = tf) # save current locmap in a temp file
#' head(vfread(tf)) # read the temp file
#' # join visual fields datasets
#' vfjoin(vfctrSunyiu24d2, vfpwgRetest24d2)
#' # visual field subselection
#' vffilter(vf, id == 1) # fields corresponding to a single subject
#' vffilter(vf, id == 1 & eye == "OD") # fields for a single subject's right eye
#' unique(vffilter(vf, eye == "OS")$eye) # only left eyes
#' vffilter(vfjoin(vfctrSunyiu24d2, vfpwgRetest24d2), type == "ctr") # get only controls
#' vffilter(vfjoin(vfctrSunyiu24d2, vfpwgRetest24d2), type == "pwg") # get only patients
#' # select visual fields by index
#' vfselect(vfctrSunyiu24d2, sel = c(1:4, 150))
#' # select last few visual fields per subject and eye
#' vfselect(vfpwgRetest24d2, sel = "last")
#' # select first few visual fields per subject and eye
#' vfselect(vfpwgRetest24d2, sel = "first")
#' vfselect(vfpwgRetest24d2, sel = "first", n = 5) # get the last 5 visits
#' # compute visual field statistics
#' vf  <- vfpwgSunyiu24d2
#' td  <- gettd(vf)  # get TD values
#' tdp <- gettdp(td) # get TD probability values
#' pd  <- getpd(td)  # get PD values
#' pdp <- getpdp(pd) # get PD probability values
#' gh  <- getgh(td)  # get the general height
#' g   <- getgl(vf)  # get global indices
#' gp  <- getglp(g)  # get global indices probability values
#' @param vf visual field dataset
#' @return \code{vfdesc} returns descriptive statistics of a visual field dataset
#' @export
vfdesc <- function(vf) {
  summary(vf) # PLACEHOLDER. TO BE REPLACED
}

#' @rdname vf
#' @param decreasing sort decreasing or increasing?
#' Default is increasing, that is \code{decreasing = FALSE}
#' @param ... arguments to be passed to or from methods
#' @return \code{vfsort} returns a sorted visual field dataset
#' @export
vfsort <- function(vf, decreasing = FALSE)
  return(vf[order(vf$id, vf$eye, vf$date, vf$time, decreasing = decreasing),])

#' @rdname vf
#' @param vf visual field data
#' @return \code{vfisvalid} returns \code{TRUE} or \code{FALSE}
#' @export
vfisvalid <- function(vf) {
  # check mandatory fields exist and have the correct format
  mandatory <- c("id", "eye", "date", "time", "age")
  eyecodes  <- c("OD", "OS", "OU")
  missingField <- !all(sapply(mandatory, function(field) {
    if(!(field %in% names(vf))) {
      warning(paste("Missing mandatory fields:", field, call. = FALSE))
      return(FALSE)
    }
    return(TRUE)
  }))
  if(missingField) return(FALSE)
  # no NAs allowed in mandatory fields
  nacols <- !all(sapply(mandatory, function(field) {
    if(any(is.na(vf[,field]))) {
      warning(paste("The mandatory field", field, "contains NAs"), call. = FALSE)
      return(FALSE)
    }
    return(TRUE)
  }))
  if(nacols) return(FALSE)
  # check eye has only allowed values "OD" "OS", or "OU"
  if(!all(unique(vf$eye) %in% eyecodes)) {
    warning(paste("Wrong eye code. They must be one of the following:",
                  paste(eyecodes, collapse = ", ")), call. = FALSE)
    return(FALSE)
  }
  # check date does have Date class
  if(class(vf$date) != "Date") {
    warning("Field 'date' must be sucessfully converted to 'Date' class", call. = FALSE)
    return(FALSE)
  }
  # check data structure for all locations is correct
  if((ncol(vf) - getlocini() + 1) != length(getvfcols())) {
    warning("Unexpected number of columns with visual field data", call. = FALSE)
    return(FALSE)
  }
  # check that all data columns are numeric (or are all their values NA meaning, which
  # may happen for locations to be excluded from statistical analysis due to their
  # proximity to the blind spot)
  if(!all(sapply(getvfcols(), function(loc) all(is.na(vf[,loc])) || is.numeric(vf[,loc])))) {
    warning("Columns with visual field data are non-numeric", call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#' @rdname vf
#' @param file the name of the csv file from where to read the data
#' @param dateformat format to be used for date. Its default value
#'   is \code{\%Y-\%m-\%d}
#' @param eyecodes codification for right and left eye, respectively.
#'   By default in visualField uses `\code{OD}` and `\code{OS}` for
#'   right and left eye respectively, but it is common to receive csv
#'   files with the codes `\code{R}` and `\code{L}`. The code `\code{OU}`
#'   for both eyes is also allowed
#' \code{eyecodes} should be equal to `\code{c("OD", "OS")}` or `\code{c("R", "L")}`.
#'   By default it is `\code{eyecodes = c("OD", "OS", "OU")}`
#' @param ... arguments to be passed to or from methods
#' @return \code{vfread} returns a visual field dataset
#' @export
vfread <- function(file, dateformat = "%Y-%m-%d", eyecodes = c("OD", "OS", "OU"), ...) {
  vf <- read.csv(file, stringsAsFactors = FALSE, ...)
  # reformat eye
  vf$eye[vf$eye == eyecodes[1]] <- "OD"
  vf$eye[vf$eye == eyecodes[2]] <- "OS"
  vf$eye[vf$eye == eyecodes[3]] <- "OU"
  # reformat date
  vf$date <- as.Date(vf$date, dateformat)
  if(!vfisvalid(vf)) warning("visual field dataset read with warnings. Check the loaded data")
  return(vf)
}

#' @rdname vf
#' @param file the name of the csv file where to write the data
#' @return \code{vfwrite} No return value
#' @export
vfwrite <- function(vf, file, dateformat = "%Y-%m-%d", eyecodes = c("OD", "OS", "OU"), ...) {
  # change date format and eye codes
  vf$date <- format(vf$date, dateformat)
  vf$eye[vf$eye == "OD"] <- eyecodes[1]
  vf$eye[vf$eye == "OS"] <- eyecodes[2]
  vf$eye[vf$eye == "OU"] <- eyecodes[3]
  write.csv(vf, file, row.names = FALSE, ...)
}

#' @rdname vf
#' @param vf1,vf2 the two visual field data objects to join or merge
#' @return \code{vfjoin} returns a visual field dataset
#' @export
vfjoin <- function(vf1, vf2) {
  # join rows of info tables together
  vf <- rbind(vf1, vf2)
  # check key
  if(any(duplicated(data.frame(vf$id, vf$eye, vf$date, vf$time))))
    stop("Cannot join visual field data if result yields duplicated keys")
  return(vfsort(vf))
}

#' @rdname vf
#' @return \code{vffilter} returns a visual field dataset
#' @export
vffilter <- function(vf, ...) return(vfsort(return(filter(vf, ...))))

#' @rdname vf
#' @param sel it can be two things, an array of indices to select from visual field data
#'   or a string with the values `\code{first}` or `\code{last}` indicating that only the
#'   first  few n visits per subject `\code{id}` and `\code{eye}` are to be  selected.
#'   Default is `\code{last}`.
#' @param n number of visits to select. Default value is 1, but it is ignored if
#'   \code{sel} is an index array
#' @return \code{vfselect} returns a visual field dataset
#' @export
vfselect <- function(vf, sel = "last", n = 1) {
  if(is.numeric(sel)) {
    if(!is.vector(sel))
      stop("wrong format of selection: sel ought to be a number of vector of numbers")
    if(min(sel) < 1 || max(sel) > nrow(vf))
      stop("index out of bounds")
    vf <- vf[unique(sel),]
  } else {
    if(!(is.atomic(sel) && length(sel)) && !is.character(sel))
      stop("wrong type of selection")
    if(sel != "last" & sel != "first")
      stop("wrong type of selection")
    # sort the data increasing or decreasing to select the first or last n visits
    if(sel == "last") {
      decreasing <- TRUE
    } else decreasing <- FALSE
    vf <- vfsort(vf, decreasing = decreasing)
    # for each unique subject id and eye return the first n visits (or all if there
    # are fewer visits than n)
    uid <- unique(data.frame(id = vf$id, eye = vf$eye))
    vf <- do.call(rbind.data.frame, lapply(1:nrow(uid), function(ii) {
      idx <- which(vf$id == uid$id[ii] & vf$eye == uid$eye[ii])
      if(length(idx) < n) return(vf[idx,])
      else return(vf[idx[1:n],])
    }))
    return(vfsort(vf))
  }
  return(vfsort(vf))
}

#' @rdname vf
#' @return \code{gettd} returns a visual field dataset with total deviation values
#' @export
gettd <- function(vf) {
  if(!vfisvalid(vf)) stop("cannot compute TD values")
  return(getnv()$tdfun(vf))
}

#' @rdname vf
#' @param td total-deviation (TD) values
#' @return \code{gettdp} returns a visual field dataset with total deviation probability values
#' @export
gettdp <- function(td) {
  if(!vfisvalid(td)) stop("cannot compute TD probability values")
  return(getnv()$tdpfun(td))
}

#' @rdname vf
#' @return \code{getpd} returns a visual field dataset with pattern deviation values
#' @export
getpd <- function(td) {
  if(!vfisvalid(td)) stop("cannot compute PD values")
  return(getnv()$pdfun(td))
}

#' @rdname vf
#' @param pd pattern-deviation (PD) values
#' @return \code{getpdp} returns a visual field dataset with pattern deviation probability values
#' @export
getpdp <- function(pd) {
  if(!vfisvalid(pd)) stop("cannot compute PD probability values")
  return(getnv()$pdpfun(pd))
}

#' @rdname vf
#' @return \code{getgh} returns the general height of visual fields tests
#' @export
getgh <- function(td) {
  if(!vfisvalid(td)) stop("cannot compute the general height (GH)")
  return(getnv()$ghfun(td))
}

#' @rdname vf
#' @return \code{getgl} returns visual fields global indices
#' @export
getgl <- function(vf) {
  if(!vfisvalid(vf)) stop("cannot compute global indices")
  td  <- gettd(vf)
  tdp <- gettdp(td)
  pd  <- getpd(td)
  pdp <- getpdp(pd)
  gh  <- getgh(td)
  return(getnv()$glfun(vf, td, pd, tdp, pdp, gh))
}

#' @rdname vf
#' @param g global indices
#' @return \code{getglp} returns probability values of visual fields global indices
#' @export
getglp <- function(g)
  return(getnv()$glpfun(g))