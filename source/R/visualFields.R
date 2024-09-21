#' @rdname visualFields
#' @title visualFields: statistical methods for visual fields
#' @description visualFields is a collection of tools for analyzing the field of vision.
#' It provides a framework for development and use of innovative methods
#' for visualization, statistical analysis, and clinical interpretation
#' of visual-field loss and its change over time. It is intended to be a
#' tool for collaborative research.
#' @details The development version of visualFields 1.x, can be found in
#' \url{https://github.com/imarinfr/vf1}. For developers who want to collaborate
#' extending, updating, and patching visualFields, all necessary imports are to
#' be added to the source file \code{visualFields.R}. visualField developers can
#' use the source codes here as examples on how to craft new source code and keep
#' documentation that is consistent with the rest of the package, roxygen2, and CRAN.
#'
#' The previous version of visualFields, 0.6, is still available for use in
#' \url{https://github.com/imarinfr/vf0}, but is no longer maintained.
#'
#' This work was supported by the NIH grant number \bold{R01EY007716} and 
#' the Veterans Administration grant number \bold{I01 RX-001821-01A1}.
#' @seealso \code{OPI}: the Open Perimetry Initiative
#' \url{https://opi.lei.org.au/} and \url{https://www.optocom.es/opi/}
#' @references
#' Marín-Franch I & Swanson WH. \emph{The visualFields package: A tool for
#' analysis and visualization of visual fields}. Journal of Vision, 2013,
#' 13(4):10, 1-12
#' 
#' Turpin A, Artes PH, & McKendrick AM. \emph{The Open Perimetry Interface:
#' An enabling tool for clinical visual psychophysics}. Journal of Vision,
#' 2012, 12(11):22, 21–25
#' @import utils graphics
#' @importFrom stats sd lm quantile predict aggregate approx na.pass optimize pt var mad
#' @importFrom Hmisc wtd.mean wtd.var wtd.quantile
#' @importFrom dplyr filter
#' @importFrom gtools combinations
#' @importFrom combinat permn
#' @importFrom XML xmlParse xmlToList
#' @importFrom oro.dicom readDICOMFile
#' @importFrom grDevices colorRampPalette rgb col2rgb chull pdf dev.off
#' @importFrom polyclip polyclip polyoffset
#' @importFrom deldir deldir tile.list
#' @importFrom plotrix draw.circle draw.ellipse
#' @importFrom tools toTitleCase
#' @importFrom rlang sym
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel column selectInput div
#'             mainPanel br tabsetPanel tabPanel plotOutput htmlOutput actionButton icon
#'             reactiveVal observeEvent updateSelectInput renderText renderPlot shinyApp
#' @importFrom shinyjs useShinyjs
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom htmlTable htmlTable
#' @importFrom boot boot
#' @importFrom pracma inv interp2
"_PACKAGE"

#' @rdname vfenv
#' @title Settings in the visualField environment
#' @description Functions to set and get settings in the visualField environment
#' @details
#' \itemize{
#'   \item\code{setdefaults} sets the default location map, normative value and
#'     graphical parameters visualFields environment
#'   \item\code{setnv} sets normative values in the visualFields environment
#'   \item\code{getnv} gets current normative values from the visualFields
#'     environment
#'   \item\code{setlocmap} sets a location map in the visualFields environment
#'   \item\code{getlocmap} gets the current location map from the visualFields
#'     environment
#'   \item\code{setgpar} sets graphical parameters in the visualFields environment
#'   \item\code{getgpar} gets current graphical parameters from the visualFields
#'     environment
#'   \item\code{setlocini} sets the column where visual field data start in the
#'     visualFields environment
#'   \item\code{getlocini} gets the column where visual field data starts from
#'     the visualFields environment
#'   \item\code{getlocini} gets the column where visual field data starts from
#'     the visualFields environment
#'   \item\code{getvfcols} gets all the columns with visual field data
#' }
#' @examples
#' # get and set normative values
#' getnv()$info$name                # print name of set normative values
#' setnv(normvals$iowa_PC26_pw_cps) # set pointwise normative values
#' getnv()$info$name                # print name of set normative values
#' setdefaults()                    # return back to defaults
#' # get and set a location map
#' getlocmap()$name         # name of set normative values
#' setlocmap(locmaps$p30d2) # set the 30-2 location map
#' getlocmap()$name         # name of set normative values
#' setdefaults()            # return back to defaults
#' # get and set a graphical parameters
#' getgpar()$tess$xlim  # limits of x axis
#' setgpar(gpars$pPeri) # set graphical parameters for the Peripheral test
#' getgpar()$tess$xlim  # limits of x axis
#' setdefaults()        # return back to defaults
#' # get and set initial column for visual field data
#' getlocini()
#' getvfcols() # get columns with visual fields data
#' setlocini(15)
#' getvfcols() # get columns with visual fields data
#' setdefaults() # return back to defaults
#' @return \code{setdefaults}: No return value
#' @export
setdefaults <- function() {
  setlocini(11)                             # set starting location for visual fields data
  setlocmap(visualFields::locmaps$p24d2)    # set default location map
  setnv(visualFields::normvals$sunyiu_24d2) # set default normative values
  setgpar(visualFields::gpars$p24d2)        # set default graphical parameters
}

#' @rdname vfenv
#' @return \code{getnv}: Returns the normative value currently in used by visualFields
#' @export
getnv <- function() return(.vfenv$nv)

#' @rdname vfenv
#' @param nv normative values to to set in the visualFields environment
#' @return \code{setnv}: No return value
#' @export
setnv <- function(nv) {
  if(is.null(nv)) stop("normative values input is null")
  assign("nv", nv, envir = .vfenv)
}

#' @rdname vfenv
#' @return \code{getgpar}: Returns the graphical parameters currently in used by visualFields
#' @export
getgpar <- function() return(.vfenv$gpar)

#' @rdname vfenv
#' @param gpar structure with all graphical parameters
#' @return \code{setgpar}: No return value
#' @export
setgpar <- function(gpar) {
  if(is.null(gpar)) stop("graphical parameters input is null")
  assign("gpar", gpar, envir = .vfenv)
}

#' @rdname vfenv
#' @return \code{getlocmap}: Returns the location map currently in used by visualFields
#' @export
getlocmap <- function() return(.vfenv$locmap)

#' @rdname vfenv
#' @param locmap location map to to set in the visualFields environment
#' @return \code{setlocmap}: No return value
#' @export
setlocmap <- function(locmap) {
  if(is.null(locmap)) stop("location map input is null")
  assign("locmap", locmap, envir = .vfenv)
}

#' @rdname vfenv
#' @return \code{getlocini}: Returns the column where visual field data starts
#' @export
getlocini <- function() return(.vfenv$locini)

#' @rdname vfenv
#' @param locini column from where to start reading the visual field data
#' @return \code{setlocini}: No return value
#' @export
setlocini <- function(locini = 11) {
  if(is.null(locini)) stop("locini input is null")
  assign("locini", locini, envir = .vfenv)
}

#' @rdname vfenv
#' @return \code{getvfcols}: Returns the columns with visual field data
#' @export
getvfcols <- function() return(getlocini() - 1 + 1:nrow(getlocmap()$coord))

#' @rdname vfenv
#' @export
.vfenv <- new.env(parent = globalenv(), size = 3) # create visualFields environment

.onAttach <- function(libname, pkgname) setdefaults()