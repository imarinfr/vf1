#' @rdname vfctrIowaPC26
#' @title Central visual field
#' @description Locations of the visual field tested have eccentricities up to 26
#' degrees and were obtained with a custom static automated perimetry.
#' Data are from 98 eyes of 98 ocular healthy subjects. Each subject underwent two
#' visual field tests, one of the central visual field (64 locations within 26 degrees
#' of fixation) and one of the peripheral visual field (64 locations with eccentricity
#' from 26 to up to 81 degrees)
#' @format See section \code{Structure of visual fields data} in \code{\link{vfdesc}}
#' @details Data are for locations within the central 26 degrees. The data for locations
#' with eccentricity from 26 to up to 81 degrees are in \code{\link{vfctrIowaPeri}}.
#' This dataset of healthy eyes was used to generate the normative values
#' \code{iowa_PC26_pw}, and \code{iowa_PC26_pw_cps} included in \code{\link{normvals}}.
#' @seealso \code{\link{vfpwgSunyiu24d2}}, \code{\link{vfctrIowaPeri}},
#'          \code{\link{vfctrSunyiu10d2}}, \code{\link{vfctrSunyiu24d2}},
#'          \code{\link{vfpwgRetest24d2}}
#' @references
#' I. Marin-Franch, P. H. Artes, L. X. Chong, A. Turpin, and M. Wall.
#' \emph{Data obtained with an open-source static automated perimetry test
#' of the full visual field in healthy adults}. Data in Brief, 21:75–82, 2018.
#' @keywords dataset
"vfctrIowaPC26"
