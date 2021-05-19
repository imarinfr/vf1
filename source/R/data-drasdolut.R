#' @rdname drasdolut
#' @title Precomputed X and Y displacement of ganglion cell bodies for any given X and Y location on the retina
#' @description It contains a first list with two LUTs for the X and Y displacement of ganglion cell bodies
#' for arbitrary locations in the retina (in mm assuming 24 mm axial length). The other two elements of the list
#' contain precomputed vectors of degrees and mm on the retina for the same schematic eye, used for conversions.
#' These are used by the function vf2gc().
#' @format A large list containing
#' \describe{
#'   \item{Drasdo_LUT}{a list of four elements: xlut and ylut are 2d matrices 
#'   containing X and Y ganglion cell positions for any given location. 
#'   Xv and Yv are vectors defining the corresponding locations for the matrices along the X and Y axis.}
#'   \item{Degs}{A vector of degrees from the fovea, using a schematic eye. Corresponds to distances on the retina stored in MM}
#'   \item{MM}{A vector of MM distance from the fovea, using a schematic eye. Corresponds to distances in degrees stored in Degs}
#' }
#' @references
#' G. Montesano, G. Ometto, R. E. Hogg, L. M. Rossetti, D. F. Garway-Heath,
#' and D. P. Crabb. \emph{Revisiting the Drasdo Model: Implications for 
#' Structure-Function Analysis of the Macular Region}. 
#' Translational Vision Science and Technology,
#' 9(10):15, 2020
#' 
#' N. Drasdo, C. L. Millican, C. R. Katholi, and C. A. Curcio. \emph{The
#' length of Henle fibers in the human retina and a model of ganglion
#' receptive field density in the visual field}. Vision Research,
#' 47:2901–2911, 2007
#' @keywords ganglion cell displacement Drasdo
"drasdolut"