#' @rdname vfplots
#' @title Plots for visual fields data
#' @description Graphical tools for visualization and statistical analysis of
#' visual fields.
#' @details
#' The following functions generate plots using visual fields data
#' \itemize{
#'   \item\code{vfgpar} generates simple graphical parameters
#'   \item\code{vftess} generates a structure to handle the visual field tessellation.
#'     Check section \code{Tesselation in visualFields} below for further details
#'   \item\code{vfcolscheme} generates the structures to handle the color scheme
#'     Check section \code{Color schemes in visualFields} below for further details
#'   \item\code{vfprogcolscheme} generates the structures to handle the color scheme
#'     for progression analysis. Check section \code{Color schemes in visualFields}
#'     below for further details
#'   \item\code{vfplot} plots a single test for visual field data
#'   \item\code{vfplotsens} plots a single test for visual field sensitivity data
#'     with a grayscale where darker means greater sensitivity loss
#'   \item\code{vfplotdev} plots a single test for visual field total or pattern
#'     deviation data with probability scales represented in color
#'   \item\code{vfplotplr} plots the results of pointwise linear regression for
#'     a series of visual fields for an eye from a subject
#'   \item\code{vflegoplot} the legoplot shows the differences between the average
#'     values of visual field tests taken as baseline and those at the end of
#'     follow up
#'   \item\code{vflegoplotsens} the legoplot for visual field sensitivity data with
#'     a grayscale where darker means greater sensitivity loss
#'   \item\code{vflegoplotdev} the legoplot for visual field total or pattern
#'     deviation data with probability scales represented in color
#' }
#' @section Structure of graphical parameters:
#' Graphical parameters for visualFields must be a list containing
#' \itemize{
#'   \item\code{coord} print x and y coordinates. They could be different from the
#'     the real visual field location testing coordinates in complex visual field
#'     grids to help readability and improve visualization of statistical results
#'   \item\code{tess} tesselation for the visual field maps. Check section
#'     \code{Tesselation in visualFields}
#'   \item\code{colmap} color map representing the probability scale. Check section
#'     \code{Color schemes in visualFields}
#' }
#' A default graphical parameters can be generated with \code{generategpar}
#' @section Tesselation in visualFields:
#' A tesselation in visualFields must be defined with a list containing
#' \itemize{
#'   \item\code{xlim}, \item\code{ylim} 2-dimensional vectors containing the minimum
#'     and maximum x and y values
#'   \item\code{floor} the value to be assinged to any sensitivity value lower than
#'     \code{floor}
#'   \item\code{tiles} a list of as many tiles defining the tesselation as visual field
#'     test locations. Each element of the list is a table with x and y coordinates defining
#'     a polygon containing the corresponding test location. Each polygon is thus the tile
#'     for each visual field test location
#'   \item\code{hull} a table with x and y coordinates defining the outer hull of the
#'     tessellation
#' }
#' A default tessellation can be generated with \code{vftess}
#' @section Color schemes in visualFields:
#' A color scheme in visualFields must be defined with a list containing
#' \itemize{
#'   \item\code{map} a table mapping probabilities levels with colors defined
#'     in hexadecimal base
#'   \item\code{fun} a function that takes sensitivity values and deviation
#'     probability levels and returns the corresponding color code.
#' }
#' A default color scheme can be generated with \code{vfcolscheme}
#' @param coord print x and y coordinates. Check section
#'   \code{Structure of graphical parameters} for details
#' @param tess tesselation for the visual field maps. Check section
#'   \code{Tesselation in visualFields} for details
#' @param probs probability scale to use for TD and PD values. It is a numeric vector
#'   of probabilities with values in [0,1]. The values 0 and 1 must be included.
#'   Although not technically necessary, it would be best if it is the same as for
#'   the normative values used
#' @param cols corresponding colors for each of the probability levels
#' @param ltprobs,ltcols color map for progression with the alternative hypothesis
#'   lower than (LT)
#' @param gtprobs,gtcols color map for progression with the alternative hypothesis
#'   lower than (GT)
#' @param neprobs,necols color map for progression with the alternative hypothesis
#'   not equal (NE)
#' @param bprobs,bcols color map for progression with blth alternative
#'   hypotheses LT and GT (B for both)
#' @examples
#' # generate a structure with default graphical parameters for the 30-2 map
#' vfgpar(locmaps$p30d2$coord)
#' @return \code{vfgpar} returns a list with graphical parameters to be used for vfplots
#' @export
vfgpar <- function(coord, tess = vftess(coord),
                   probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
                   cols  = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4),
                             "#F7F0EB", # within normal limits
                             colorRampPalette(c("#00FF00", "#008000"))(4)),
                   floor = 0,
                   ltprobs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 1),
                   ltcols  = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB", "#008000"),
                   gtprobs = c(0, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
                   gtcols  = c("#000000", "#FF0000", "#F7F0EB", colorRampPalette(c("#00FF00", "#008000"))(4)),
                   neprobs = c(0, 0.0025, 0.005, 0.01, 0.25, 0.975, 0.99, 0.995, 0.9975, 1),
                   necols  = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
                               colorRampPalette(c("#FFFF00", "#FF0000"))(4)),
                   bprobs  = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
                   bcols   = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
                               colorRampPalette(c("#00FF00", "#008000"))(4)))
  return(list(coord = coord, tess = tess, colmap = vfcolscheme(probs, cols, floor),
              progcolmap = list(
                lt = vfprogcolscheme(ltprobs, ltcols),
                gt = vfprogcolscheme(gtprobs, gtcols),
                ne = vfprogcolscheme(neprobs, necols),
                b  = vfprogcolscheme(bprobs, bcols))))

#' @rdname vfplots
#' @param floor Flooring value, typically in dB. Default is 0
#' @param delta Distance over which the boundary should be shifted. See for \code{\link{polyclip}}
#' @examples
#' # generate a structure with default tesselation for the 30-2 map
#' vftess(locmaps$p30d2$coord)
#' @return \code{vftess} returns a list with the \code{xlim}, \code{ylim}, tessellation tiles and an outer hull
#' to be used for vfplots
#' @export
vftess <- function(coord, floor = 0, delta = 3) {
  # get and expand the convex hull
  hull <- coord[chull(coord),]
  hull <- as.data.frame(polyoffset(hull, delta, jointype = "round")[[1]])
  # get tiles
  tiles <- lapply(tile.list(deldir(coord)), function(tt) data.frame(x = tt$x, y = tt$y))
  # and intersect with the convex hull to obtain the
  tiles <- lapply(tiles, function(tt) data.frame(polyclip(tt, hull)))
  return(list(xlim  = c(min(hull$x), max(hull$x)), ylim  = c(min(hull$y), max(hull$y)),
              floor = floor, tiles = tiles, hull  = hull))
}

#' @rdname vfplots
#' @examples
#' # default color scheme
#' vfcolscheme()
#' @return \code{vfcolscheme} returns a list with a lookup table and a function that define the color scheme
#' to be used for vfplots
#' @export
vfcolscheme <- function(probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
                        cols  = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4),
                                  "#F7F0EB", # within normal limits
                                  colorRampPalette(c("#00FF00", "#008000"))(4)),
                        floor = 0) {
  if(any(probs < 0) || any(probs > 1))
    stop("probability values must be between 0 and 1")
  if(!any(probs == 0) || !any(probs == 1))
    stop("probability values must include values 0 and 1")
  probs <- sort(probs)
  map <- data.frame(probs = probs, cols = cols, stringsAsFactors = FALSE) # fucking hate R defaulting to factors
  fun <- as.function(alist(vf = , devp = , {
    vf[which(is.na(devp))]   <- Inf
    devp[which(is.na(devp))] <- Inf # assign infinite to locations to ignore (e.g., blind spot)
    cols <- rep(NA, length(devp))
    for(i in nrow(map):1) cols[devp <= map$probs[i]] <- map$cols[i]
    cols[vf < floor] <- "#000000" # not seen are plotted in black
    return(cols)
  }))
  return(list(map = map, fun = fun))
}

#' @rdname vfplots
#' @examples
#' # default color scheme for progression
#' vfprogcolscheme()
#' @return \code{vfprogcolscheme} returns the default \code{vfcolscheme} to be used for vfplots
#' @export
vfprogcolscheme <- function(probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 1),
                            cols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4),
                                     "#F7F0EB", "#008000")) {
  if(any(probs < 0) || any(probs > 1))
    stop("probability values must be between 0 and 1")
  if(!any(probs == 0) || !any(probs == 1))
    stop("probability values must include values 0 and 1")
  probs <- sort(probs)
  map <- data.frame(probs = probs, cols = cols, stringsAsFactors = FALSE) # fucking hate R defaulting to factors
  fun <- as.function(alist(vals = , {
    vals[which(is.na(vals))] <- Inf # Assign infinite to locations to ignore (e.g., blind spot)
    cols <- rep(NA, length(vals))
    for(i in nrow(map):1) cols[vals <= map$probs[i]] <- map$cols[i]
    return(cols)
  }))
  return(list(map = map, fun = fun))
}

#' @rdname vfplots
#' @param vf the visual fields data to plot
#' @param type the type of data to plot: sensitivities (`\code{s}`),
#' total deviation values (`\code{td}`), pattern deviation values (`\code{pd}`),
#' a hybrid plot that shows sensitivity grayscale with TD values and corresponding
#' probability levels (`\code{tds}`), or PD values and corresponding probability
#' levels (`\code{pds}`). Default is `\code{td}`.
#' @param ... other graphical arguments. See \code{\link{plot}}
#' @examples
#' # plot visual field values for the last field in the series for the first
#' # subject in the dataset vfpwgSunyiu24d2
#' # grayscale with sensitivity values
#' vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "s")
#' # TD values
#' vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "td")
#' # PD values
#' vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "pd")
#' # hybrid sensitivities and TD values
#' vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "tds")
#' # hybrid sensitivities and PD values
#' vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "pds")
#' @return \code{vfplot} No return value
#' @export
vfplot <- function(vf, type = "td", ...) {
  if(nrow(vf) != 1) stop("can plot only 1 visual field at a time")
  if(!vfisvalid(vf)) stop("incorrect visual field data structure. Cannot plot")
  nv   <- getnv()   # get normative values
  gpar <- getgpar() # get graphical parameters
  locs <- getlocini():ncol(vf)
  # left or right eye
  if(vf$eye == "OS") gpar$tess$xlim <- gpar$tess$xlim[2:1]
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(mar = c(0, 0, 0, 0), ...)
  # maximum values for grayscales used in type "s", "tds", and "pds"
  maxdb <- nv$agem$model(vf$age)
  # for locations in the blind spot, we input the values of the previous locations 
  maxdb[which(is.na(maxdb))] <- maxdb[which(is.na(maxdb)) - 1]
  # dispatch to raw sensitivity (grayscale) plots, TD and PD (color) plots, or the hybrid
  # sensitivity with TD and PD plot
  if(type == "s") # sensitivities
    vfplotsens(gpar, vf[,locs], maxdb, ...)
  else { # TD, PD or hybrid plots
    if(type %in% c("td", "tds")) {
      dev  <- gettd(vf)
      devp <- gettdp(dev)
    } else if(type %in% c("pd", "pds")) {
      dev  <- getpd(gettd(vf))
      devp <- getpdp(dev)
    } else stop("wrong type of plot requested. Must be 's', 'td', 'pd', 'td', or 'pds'")
    if(type %in% c("td", "pd"))
      vfplotdev(gpar, vf[,locs], dev[,locs], devp[,locs], ...)
    else
      vfplotsdev(gpar, vf[,locs], maxdb, dev[,locs], devp[,locs], ...)
  }
}

#' @rdname vfplots
#' @param gpar graphical parameters
#' @param maxval maximum value, typically in dB, for the generation of a grayscale
#' @param digits digits to round values to plot. Default is 0
#' @return \code{vfplotsens} No return value
#' @export
vfplotsens <- function(gpar, vf, maxval, digits = 0, ...) {
  # background gray shades
  fcol <- (vf - gpar$tess$floor) / (maxval - gpar$tess$floor)
  fcol[fcol > 1] <- 1
  fcol[fcol < 0] <- 0
  # foreground text gray shades
  tcol <- rep(0.3, length(vf))
  tcol[fcol < 0.5] <- 0.7
  # blind spot
  fcol[getlocmap()$bs] <- 1
  tcol[getlocmap()$bs] <- 0.3
  # convert to hexadecimal color
  fcol <- rgb(fcol, fcol, fcol)
  tcol <- rgb(tcol, tcol, tcol)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE,
       axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles))
    polygon(gpar$tess$tiles[[i]], col = fcol[i], border = NA)
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")
  txt <- round(vf, digits)
  txt[is.na(vf)] <- ""
  txt[vf < gpar$tess$floor] <- paste0("<", gpar$tess$floor)
  text(gpar$coord$x, gpar$coord$y, txt, col = tcol, ...)
}

#' @rdname vfplots
#' @param dev deviation (TD or PD) values
#' @param devp deviation (TD or PD) probability values
#' @return \code{vfplotdev} No return value
#' @export
vfplotdev <- function(gpar, vf, dev, devp, digits = 0, ...) {
  # background colors and foreground text gray shades
  cols <- gpar$colmap$fun(vf, devp)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE,
       axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles)) {
    tiles  <- gpar$tess$tiles[[i]] # get tiles
    otiles <- polyoffset(tiles, -0.75, jointype = "round")[[1]] # shrink tiles to plot white region
    polygon(tiles, col = cols[i], border = "lightgray")  # plot tiles with grayscales
    polygon(otiles, col = "white", border = NA) # plot tiles with white background to show text
  }
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)  
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")   
  txt <- round(dev, digits)
  txt[is.na(dev)] <- ""
  text(gpar$coord$x, gpar$coord$y, txt, col = rgb(0.3, 0.3, 0.3), ...)
}

#' @rdname vfplots
#' @return \code{vfplotsdev} No return value
#' @export
vfplotsdev <- function(gpar, vf, maxval, dev, devp, digits = 0, ...) {
  # background gray shades
  fcol <- (vf - gpar$tess$floor) / (maxval - gpar$tess$floor)
  fcol[fcol > 1] <- 1
  fcol[fcol < 0] <- 0
  # foreground text gray shades
  tcol <- rep(0.3, length(vf))
  tcol[fcol < 0.5] <- 0.7
  # blind spot
  fcol[getlocmap()$bs] <- 1
  tcol[getlocmap()$bs] <- 0.3
  # convert to hexadecimal color
  fcol <- rgb(fcol, fcol, fcol)
  tcol <- rgb(tcol, tcol, tcol)
  # background colors and foreground text gray shades
  cols <- gpar$colmap$fun(vf, devp)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE,
       axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles)) {
    tiles  <- gpar$tess$tiles[[i]] # get tiles
    otiles <- polyoffset(tiles, -0.75, jointype = "round")[[1]] # shrink tiles to plot white region
    polygon(tiles, col = cols[i], border = "lightgray")  # plot tiles with grayscales
    polygon(otiles, col = fcol[i], border = NA) # plot tiles with white background to show text
  }
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)  
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")   
  txt <- round(dev, digits)
  txt[is.na(dev)] <- ""
  text(gpar$coord$x, gpar$coord$y, txt, col = tcol, ...)
}

#' @rdname vfplots
#' @param alternative alternative hypothesis used in progression analyses.
#'   Allowed values are `\code{LT}` (as in "lower than", default),
#'   `\code{GT}` (as in "greater than"), `\code{NE}` (as in "not equal"),
#'   and `\code{both}` (both `\code{LT}` and `\code{GT}`)
#' @param xoffs,yoffs offset x and y where to print the slope values. That is,
#' the distance from the center of each Voronoy polygons in degrees of visual angle
#' @examples
#' # plot results from pointwise linear regression for the series of
#' # visual fields for the right eye in the dataset vfpwgSunyiu24d2
#' # with sensitivity values
#' vfplotplr(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "s")
#' # TD values
#' vfplotplr(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "td")
#' # PD values
#' vfplotplr(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "pd")
#' @return \code{vfplotplr} No return value
#' @export
vfplotplr <- function(vf, type = "td", alternative = "LT", xoffs = 0, yoffs = 0, ...) {
  res <- plr(vf, type) # if more than 1 ID/eye then it crashes as it should
  gpar <- getgpar() # get graphical parameters
  # left or right eye
  if(vf$eye[1] == "OS") gpar$tess$xlim <- gpar$tess$xlim[2:1]
  # format intercept and slope
  sl <- format(round(res$sl, 1), nsmall = 1L)
  sl[sl == "NA"] <- ""
  if(alternative == "LT") {
    cols <- getgpar()$progcolmap$lt$fun(res$pval)
  } else if(alternative == "GT") {
    cols <- getgpar()$progcolmap$gt$fun(res$pval)
  } else if(alternative == "NE") {
    cols <- getgpar()$progcolmap$ne$fun(res$pval)
  } else if(alternative == "both") {
    cols <- getgpar()$progcolmap$b$fun(res$pval)
  } else stop("wrong alternative")
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(mar = c(0, 0, 0, 0), ...)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE, axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles)) {
    tiles  <- gpar$tess$tiles[[i]] # get tiles
    otiles <- polyoffset(tiles, -0.75, jointype = "round")[[1]] # shrink tiles to plot white region
    polygon(tiles, col = cols[i], border = "lightgray")  # plot tiles with grayscales
    polygon(otiles, col = "white", border = NA) # plot tiles with white background to show text
  }
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")
  text(gpar$coord$x + xoffs, gpar$coord$y + yoffs, sl, col = rgb(0.3, 0.3, 0.3), ...)
}

#' @rdname vfplots
#' @param grp number of baseline (first) and last visual fields to group.
#' Default is `\code{3}`
#' @examples
#' # legoplot for the series of visual fields for the right eye
#' # of the subject in the dataset vfpwgSunyiu24d2
#' # with sensitivity values
#' vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "s")
#' # TD values
#' vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "td")
#' # PD values
#' vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "pd")
#' @return \code{vflegoplot} No return value
#' @export
vflegoplot <- function(vf, type = "td", grp = 3, ...) {
  if(nrow(unique(data.frame(vf$id, vf$eye))) != 1)
    stop("all visual fields must belong to the same subject id and eye")
  nv   <- getnv()
  gpar <- getgpar() # get graphical parameters
  locs <- getlocini():ncol(vf)
  # left or right eye
  if(vf$eye[1] == "OS") gpar$tess$xlim <- gpar$tess$xlim[2:1]
  vfb <- vfmean(vfselect(vf, sel = "first", n = grp), by = "eye")
  vfl <- vfmean(vfselect(vf, sel = "last", n = grp), by = "eye")
  if(type == "s") {
    maxb <- nv$agem$model(vfb$age)
    maxl <- nv$agem$model(vfb$age)
    vflegoplotsens(gpar, vfb[,locs], vfl[,locs], maxb, maxl, ...)
  } else {
    if(type == "td") {
      devb  <- gettd(vfb)
      devbp <- gettdp(devb)
      devl  <- gettd(vfl)
      devlp <- gettdp(devl)
    } else if(type == "pd") {
      devb  <- getpd(gettd(vfb))
      devbp <- getpdp(devb)
      devl  <- getpd(gettd(vfl))
      devlp <- getpdp(devl)
    } else stop("wrong type of plot requested. Must be 's', 'td', or 'pd'")
    vflegoplotdev(gpar,
                  vfb[,locs], devb[,locs], devbp[,locs],
                  vfl[,locs], devl[,locs], devlp[,locs], ...)
  }
}

#' @rdname vfplots
#' @param vfb baseline visual field data
#' @param vfl last visual field data
#' @param maxb maximum value for the grayscale at baseline visual field data
#' @param maxl maximum value for the grayscale for last visual field data
#' @param crad radius of the circle in the legoplot
#' @return \code{vflegoplotsens} No return value
#' @export
vflegoplotsens <- function(gpar, vfb, vfl, maxb, maxl, crad = 2, digits = 1, ...) {
  bs <- getlocmap()$bs
  # baseline grayscale
  fcolb <- (vfb - gpar$tess$floor) / (maxb - gpar$tess$floor)
  fcolb[fcolb > 1] <- 1
  fcolb[fcolb < 0] <- 0
  fcolb[bs] <- 1   # blind spot
  fcolb <- rgb(fcolb, fcolb, fcolb)
  # final grayscale
  fcol <- (vfl - gpar$tess$floor) / (maxl - gpar$tess$floor)
  fcol[fcol > 1] <- 1
  fcol[fcol < 0] <- 0
  fcol[bs] <- 1
  # foreground text gray shades
  tcol <- rep(0.3, length(vfl))
  tcol[fcol < 0.5] <- 0.7
  # convert to hexadecimal color
  fcol <- rgb(fcol, fcol, fcol)
  tcol <- rgb(tcol, tcol, tcol)
  # values to present
  txt <- round(vfl - vfb, digits)
  # remove blind spot locations
  coord <- gpar$coord
  if(length(bs) > 0) {
    coord <- coord[-bs,]
    txt   <- txt[-bs]
    fcol  <- fcol[-bs]
    tcol  <- tcol[-bs]
  }
  # plot
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(mar = c(0, 0, 0, 0), ...)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE, axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles))
    polygon(gpar$tess$tiles[[i]], border = "lightgray", col = fcolb[i])
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")
  # draw circles
  for(i in 1:nrow(coord))
    draw.circle(coord$x[i], coord$y[i], radius = crad, col = fcol[i], lty = 0)
  text(coord$x, coord$y, txt, col = tcol, ...)
}

#' @rdname vfplots
#' @param devb baseline visual field (TD or PD) deviation values
#' @param devpb baseline visual field (TD or PD) deviation probability values
#' @param devl last visual field (TD or PD) deviation values
#' @param devpl last visual field (TD or PD) deviation probability values
#' @return \code{vflegoplotdev} No return value
#' @export
vflegoplotdev <- function(gpar, vfb, devb, devpb, vfl, devl, devpl, crad = 2, digits = 1, ...) {
  bs <- getlocmap()$bs
  # baseline colors
  colsb <- gpar$colmap$fun(vfb, devpb)
  # final colors
  colsl <- gpar$colmap$fun(vfl, devpl)
  # values to present
  txt <- round(devl - devb, digits)
  # text color
  tcol   <- rep(0.3, length(vfl))
  colrgb <- col2rgb(colsl) / 255
  tcol[(0.2126 * colrgb[1,] + 0.7152 * colrgb[2,] + 0.0722 * colrgb[3,]) < 0.4] <- 0.7
  tcol <- rgb(tcol, tcol, tcol)
  # remove blind spot locations
  coord <- gpar$coord
  if(length(bs) > 0) {
    coord <- coord[-bs,]
    txt   <- txt[-bs]
    colsl <- colsl[-bs]
    tcol  <- tcol[-bs]
  }
  # plot
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(mar = c(0, 0, 0, 0), ...)
  plot(gpar$coord$x, gpar$coord$y, typ = "n", ann = FALSE, axes = FALSE, asp = 1,
       xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
  # plot polygons
  for(i in 1:length(gpar$tess$tiles))
    polygon(gpar$tess$tiles[[i]], border = "lightgray", col = colsb[i])
  # add blind spot
  draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)
  # outer hull
  polygon(gpar$tess$hull, border = "lightgray")
  # draw circles
  for(i in 1:nrow(coord))
    draw.circle(coord$x[i], coord$y[i], radius = crad, col = colsl[i], lty = 0)
  text(coord$x, coord$y, txt, col = tcol, ...)
}

#' @rdname vfplots
#' @param thr threshold used for the median absolute deviation of residuals
#' from simple linear regression. If greater than the threshold, the
#' sparkline for that location is plotted in red and with a thicker line.
#' Default is `\code{2}` (dB)
#' @param width the width of each pointwise sparkline plot. Default is
#' `\code{4}` (degrees of visual angle)
#' @param height the height of each pointwise sparkline plot. Default is
#' `\code{2}` (degrees of visual angle)
#' @param add whether to generate a new plot (`\code{FALSE}`, as default)
#' or to add to an existing figure (`\code{TRUE}`)
#' @examples
#' # sparklines for the series of visual fields for the right eye of
#' # the subject in the dataset vfpwgSunyiu24d2
#' # with sensitivity values
#' vfplotsparklines(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "s")
#' # TD values
#' vfplotsparklines(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "td")
#' # PD values
#' vfplotsparklines(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "pd")
#' @return \code{vfplotsparklines} No return value
#' @export
vfplotsparklines <- function(vf, type = "td", thr = 2, width = 4,
                             height = 2, add = FALSE, ...) {
  if(nrow(unique(data.frame(vf$id, vf$eye))) != 1)
    stop("all visual fields must belong to the same subject id and eye")
  nv   <- getnv()
  gpar <- getgpar() # get graphical parameters
  locs <- getlocini():ncol(vf)
  # left or right eye
  x <- as.numeric(vf$date - vf$date[1]) / 365.25 # it should be difference in years from baseline date
  if(type == "s") {
    y <- vf[,locs]
  } else if(type == "td") {
    y <- gettd(vf)[,locs]
  } else if(type == "pd") {
    y <- getpd(gettd(vf))[,locs]
  } else stop("wrong type of plot requested. Must be 's', 'td', or 'pd'")
  # remove blind spot locations
  if(length(getlocmap()$bs) > 0) {
    gpar$coord <- gpar$coord[-getlocmap()$bs,]
    y <- y[,-getlocmap()$bs]
  }
  # left or right eye
  if(vf$eye[1] == "OS") gpar$tess$xlim <- gpar$tess$xlim[2:1]
  xlim <- c(0, max(x))
  ylim <- c(min(y), max(y))
  suppressWarnings(resmad <- sapply(as.list(y), function(y) mad(lm(y ~ x)$residuals)))
  cols <- rep("#4D4D4D", length(resmad))
  cols[resmad > thr] <- "#FF0000"
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  if(!add) {
    par(mar = c(0, 0, 0, 0), ...)
    plot(gpar$coord$x, gpar$coord$y, typ = "n",
         ann = FALSE, axes = FALSE, asp = 1,
         xlim = gpar$tess$xlim, ylim = gpar$tess$ylim, ...)
    # plot polygons
    lapply(gpar$tess$tiles, polygon, border = "lightgray", col = "#FFFFFF")
    # add blind spot
    draw.ellipse(15, -1.5, 2.75, 3.75, col = "lightgray", border = NA)
    # outer hull
    polygon(gpar$tess$hull, border = "lightgray")
  }
  # plot the spark lines:  for left eyes, handling figure positions is incompatible
  # with swapping the min and max x limits in the plot. We need a patch here
  if(gpar$tess$xlim[1] < gpar$tess$xlim[2]) {
    figs <- cbind(grconvertX(gpar$coord$x - width  / 2, to = "ndc"),
                  grconvertX(gpar$coord$x + width  / 2, to = "ndc"),
                  grconvertY(gpar$coord$y, to = "ndc"),
                  grconvertY(gpar$coord$y + height, to = "ndc"))
  } else {
    figs <- cbind(grconvertX(gpar$coord$x + width  / 2, to = "ndc"),
                  grconvertX(gpar$coord$x - width  / 2, to = "ndc"),
                  grconvertY(gpar$coord$y, to = "ndc"),
                  grconvertY(gpar$coord$y + height, to = "ndc"))
  }
  for(i in 1:nrow(figs)) {
    par(fig = figs[i,], new = TRUE)
    plot(x, y[,i], type = "l", xlim = xlim, ylim = ylim, axes = FALSE, col = cols[i], ...)
  }
}
