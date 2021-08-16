#' @rdname locmap
#' @title Locmap management
#' @description Functions to handle location maps, which are lists with x and y
#' coordinates and other importan information about the visual field test
#' locations. Check section \code{Structure of location maps} below for details
#' @details
#' \itemize{
#'     \item\code{locread} reads a csv file with location map data
#'     \item\code{locwrite} writes a csv file with location map data
#' }
#' @section Structure of location maps:
#' Each element in the list \code{locmaps} is a location map that
#' contains the following fields
#' \itemize{
#'   \item\code{name} descriptive name
#'   \item\code{desc} brief description
#'   \item\code{coord} coordinates of the visual field locations
#'   \item\code{bs} if not empty, the locations that ought to be removed
#'     for statistical analysis due to their proximity to the blind spot
#' }
#' @examples
#' # write and read location map
#' tf <- tempfile("locmap")
#' locwrite(getlocmap(), file = tf) # save current locmap in a temp file
#' print(locread(tf, name = "name", desc = "desc", bs = c(1, 2))) # read the temp file
#' @param file the name of the file which the data are to be read from
#' @param name to give the location map
#' @param desc brief description for the location map
#' @param bs locations that should be excluded from statistical analysis because
#'   of their proximity to the blind spot
#' @param ... arguments to be passed to or from methods
#' @return \code{locread} a list with information about a location map
#' @export
locread <- function(file, name = "", desc = "", bs = numeric(), ...)
  return(list(name  = name,
              desc  = desc,
              coord = read.csv(file, stringsAsFactors = FALSE, ...),
              bs    = bs))

#' @rdname locmap
#' @param locmap location map from which to get coordinates to export as csv file
#' @return \code{locwrite} No return value
#' @export
locwrite <- function(locmap, file, ...)
  write.csv(locmap$coord, file, row.names = FALSE, ...)