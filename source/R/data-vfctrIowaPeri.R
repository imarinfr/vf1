#' @rdname vfctrIowaPeri
#' @title Peripheral visual field
#' @description Locations of the visual field tested have eccentricities from 26
#' to up to 81 degrees and were obtained with a custom static automated perimetry.
#' Data are from 98 eyes of 98 ocular healthy subjects. Each subject underwent two
#' visual field tests, one of the central visual field (64 locations within 26 degrees
#' of fixation) and one of the peripheral visual field (64 locations with eccentricity
#' from 26 to up to 81 degrees)
#' @format See section \code{Structure of visual fields data} in \code{\link{vfdesc}}
#' @details Data are for locations with eccentricity from 26 to up to 81 degrees.
#' The dataset for locations within the central 26 degrees are in \code{\link{vfctrIowaPC26}}.
#' This dataset of healthy eyes was used to generate the normative values
#' \code{iowa_Peri_pw}, and \code{iowa_Peri_pw_cps} included in \code{\link{normvals}}.
#' @seealso \code{\link{vfpwgSunyiu24d2}}, \code{\link{vfctrIowaPC26}},
#'          \code{\link{vfctrSunyiu10d2}}, \code{\link{vfctrSunyiu24d2}},
#'          \code{\link{vfpwgRetest24d2}}
#' @references
#' I. Marin-Franch, P. H. Artes, L. X. Chong, A. Turpin, and M. Wall.
#' \emph{Data obtained with an open-source static automated perimetry test
#' of the full visual field in healthy adults}. Data in Brief, 21:75–82, 2018.
#' @keywords dataset
"vfctrIowaPeri"
