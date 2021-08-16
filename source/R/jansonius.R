#' @rdname jansonius
#' @title The Jansonius map for average path of nerve fiber bundles
#' @description Generates a function that renders the average path of
#' a nerve fiber bundle that exits through the optic nerve head (ONH)
#' with a particular angle
#' @details
#' \itemize{
#'   \item\code{cart2jpolar} converts the cartesian coordinates to the polar
#'     coordinates in the distorted space used in the Jansonius map
#'   \item\code{jpolar2cart} converts back from the Jansonius polar
#'     coordinates to cartesian coordinates
#'   \item\code{bundlePath} returns a function describing the expected fiber
#'     path given an angle of incidence on the ONH
#'   \item\code{loc2psi} returns the angle of incidence of the average bundle path
#'     that passes through specific locations of the visual field
#'   \item\code{psi2oct} returns the angle of OCT circular scans corresponding
#'     to average bundle paths with specific angle of incidence at the ONH
#'   \item\code{vf2gc} calculates ganglion-cell soma locations
#' }
#' @references
#' N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill,
#' P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze,
#' P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold,
#' and U. Schiefer. \emph{A mathematical description of nerve
#' fiber bundle trajectories and their variability in the human
#' retina}. Vision Research, 49(17):2157-2163, 2009
#' 
#' N. M. Jansonius, J. Nevalainen, B. Selig, L. M. Zangwill,
#' P. A. Sample, W. M. Budde, J. B. Jonas, W. A. Lagreze,
#' P. J. Airaksinen, R. Vonthein, L. A. Levin, J. Paetzold,
#' and U. Schiefer. \emph{Erratum to "A mathematical description
#' of nerve fiber bundle trajectories and their variability in
#' the human retina"}. Vision Research, 50:1501, 2010
#' 
#' N. M. Jansonius, J. Schiefer, J. Nevalainen, J. Paetzold,
#' and U. Schiefer. \emph{A mathematical model for describing
#' the retinal nerve fiber bundle trajectories in the human eye:
#' Average course, variability, and influence of refraction, optic
#' disc size and optic disc position}. Experimental Eye Research,
#' 105:70-78, 2012
#' 
#' N. Drasdo, C. L. Millican, C. R. Katholi, and C. A. Curcio. \emph{The
#' length of Henle fibers in the human retina and a model of ganglion
#' receptive field density in the visual field}. Vision Research,
#' 47:2901–2911, 2007
#' 
#' D. C. Hood, A. S. Raza, D. M. C. G. V., J. G. Odel, V. C. Greenstein,
#' J. M. Liebmann, and R. Ritch. \emph{Initial arcuate defects within the
#' central 10  degrees in glaucoma}. Investigative Ophthalmology and Visual
#' Science, 52(2):940-946, 2011
#' 
#' A. S. Raza, J. Cho, D. M. C. G. V., H. Wang, X. Zhang, R. H. Kardon,
#' J. M. Liebmann, R. Ritch, and D. C. Hood. \emph{Retinal ganglion cell
#' layer thickness and local visual field sensitivity in glaucoma}.
#' Archives of Ophthalmology, 129(12):1529-1536, 2011
#' 
#' G. Montesano, G. Ometto, R. E. Hogg, L. M. Rossetti, D. F. Garway-Heath,
#' and D. P. Crabb. \emph{Revisiting the Drasdo Model: Implications for 
#' Structure-Function Analysis of the Macular Region}. 
#' Translational Vision Science and Technology,
#' 9(10):15, 2020
#' @examples
#' # get ganglion-cell soma locations from visual field locations
#' vf2gc(locmaps$p10d2$coord)
#' # convert to polar of the distorted space used by Jansonius map and back
#' coord <- data.frame(x = c(3, 0, -3), y = c(0, 0, 0))
#' (rpsi <- cart2jpolar(coord))
#' jpolar2cart(rpsi)
#' 
#' # get an average bundle path from a specific angle of incidence in the ONH
#' # The object returned is a function that returns polar angles of the
#' # distorted space of the Jansonius map for distances from the ONH center
#' pathFun <- bundlePath(-125)
#' jpolar2cart(data.frame(10:20, pathFun(10:20)))
#' 
#' # get angle of incidence in the ONH from locations of the visual field
#' loc2psi(coord)
#' 
#' # get the OCT circular scan angles from the angle of incidence in the ONH
#' # for the 10-2 map of locations, ...
#' psi2oct(loc2psi(locmaps$p10d2$coord))
#' # the previous operation was actually fundamentally wrong! We need to
#' # obtain first the 
#' psi2oct(loc2psi(vf2gc(locmaps$p10d2$coord)))
#' @param coord coordinates of locations in the visual field
#' @return \code{cart2jpolar}: returns the Jansonius modified polar coordinates
#' @export
cart2jpolar <- function(coord) {
  x  <- coord[,1]
  y  <- coord[,2]
  xp <- x - 15
  yp <- y
  yp[x > 0] <- y[x > 0] - 2 * (x[x > 0] / 15)^2
  r <- sqrt(xp^2 + yp^2)
  psi <- atan2(yp, xp)
  return(data.frame(r = r, psi = 180 / pi * psi))
}

#' @rdname jansonius
#' @param rpsi visual field locations in polar coordinates of the
#'   distorted space of the Jansonius map
#' @return \code{jpolar2cart}: returns Cartesian coordinates
#' @export
jpolar2cart <- function(rpsi) {
  r   <- rpsi[,1]
  psi <- rpsi[,2]
  psi <- (pi / 180 * psi) %% (2 * pi)
  x <- sqrt(r^2 / (1 + (tan(psi))^2))
  y <- abs(x * tan(psi))
  x[psi > (pi / 2) & psi < (3 * pi / 2)] <- -x[psi > (pi / 2) & psi < (3 * pi / 2)]
  y[psi > pi] <- -y[psi > pi]
  x <- x + 15
  y[x > 0] <- y[x > 0] + 2 * (x[x > 0] / 15)^2
  return(data.frame(x = round(x, 6), y = round(y, 6)))
}

#' @rdname jansonius
#' @param psi0 angle of incidence at the ONH
#' @param r0 radius of the ONH. Its default value is \code{4}.
#'   Changing it changes the calculated average bundle paths.
#' @return \code{bundlePath}: returns a function describing a retinal ganglion cell bundle path
#' @export
bundlePath <- function(psi0, r0 = 4) {
  if(!(is.atomic(psi0) && length(psi0) == 1L) ||
     !is.numeric(psi0) || abs(psi0) > 180)
    stop("input argument psi0 must be a numeric scalar from -180 to 180")
  if(psi0 == -180) psi0 <- 180
  if(psi0 >= 0 & psi0 < 60) {
    b <- 0.00083 * psi0^2 + 0.020 * psi0 - 2.65
    c <- 1.9 + 1.4 * tanh((psi0 - 121) / 14)
  }
  if(psi0 >= 60 & psi0 <= 180) {
    b <- exp(-1.9 + 3.9 * tanh(-(psi0 - 121) / 14))
    c <- 1.9 + 1.4 * tanh((psi0 - 121) / 14)
  }
  if(psi0 > -180 & psi0 <= -60) {
    b <- -exp(0.7 + 1.5 * tanh(-(-psi0 - 90) / 25))
    c <- 1.0 + 0.5 * tanh((-psi0 - 90) / 25)
  }
  if(psi0 > -60 & psi0 < 0) {
    b <- 0.00083 * psi0^2 + 0.020 * psi0 - 2.65
    c <- 1.0 + 0.5 * tanh((-psi0 - 90) / 25)
  }
  fbfun <- function(r) {
    return(psi0 + b * (r - r0)^c)
  }
  return(fbfun)
}

#' @rdname jansonius
#' @param diam diamater in degrees of visual angle of the OCT circular
#'   scan centered at the center of the ONH
#' @return \code{loc2psi}: returns the angle of incidence on the ONH
#' @export
loc2psi <- function(coord, r0 = 4) {
  coord[,2] <- -coord[,2]
  psi0 <- sapply(split(cart2jpolar(coord), seq(nrow(coord))), function(rpsi) {
    fbpathinv <- function(psi0) {
    fb <- bundlePath(psi0, r0)
    psiest <- fb(rpsi$r)
    return((psiest - rpsi$psi)^2)}
    optimize(fbpathinv, interval = c(-179.99, 180))$minimum})
  return(as.numeric(psi0))
}

#' @rdname jansonius
#' @return \code{psi2oct}: returns the corresponding angle in the OCT circular scan
#' @export
psi2oct <- function(psi0, diam = 12) {
  radius <- diam / 2
  r <- seq(4, 10, 0.001)
  return(180 / pi * sapply(psi0, function(psi0) {
    fb <- bundlePath(psi0)
    coord <- jpolar2cart(data.frame(r = r, phi = fb(r)))
    coord <- coord - c(15, 2)
    coord <- coord[which.min(abs(coord$x^2 + coord$y^2 - radius^2)),]
    oct <- atan2(coord$y, -coord$x)
    oct[oct < 0] <- 2 * pi + oct[oct < 0]
    return(oct)}))
}

#' @rdname jansonius
#' @param angle fovea-disc angle in degrees
#' @return \code{vf2gc}: returns the ganglion cell soma corresponding to the photoreceptors
#'                       of a visual field location
#' @export
vf2gc <- function(coord, angle = 0) {
  Drasdo_LUT <- visualFields::drasdolut$Drasdo_LUT
  Degs <- visualFields::drasdolut$Degs
  MM <- visualFields::drasdolut$MM
  
  ############################################################################################################
  #Input to the function must be in visual field coordinates for a right eye:
  # - Negative Y is superior on the retina; 
  # - Negative X is temporal on the retina
  
  # Fovea is assumed to be at {0, 0}
  # The fovea-disc angle is 0 deg in the original Curcio and Allen map. If changed, it needs to be
  # in retinal coordinates for a right eye, meaning that a positive angle indicates a positive elevation
  # of the optic nerve over the fovea. This angle must be in degrees.
  
  #Assumes an axial length of 24 mm. Displacements in degrees are the same for any axial length
  #under the assumption of spherical global expansion and no changes to the anterior chamber of the schematic eye
  ############################################################################################################
  
  x <- coord[,1]
  y <- coord[,2]
  
  #Inverts y axis to convert to retinal coordinates (necessary for next steps)
  y <- -y
  #Converts Fovea-Disc angle to radians for following calculations
  angle <- angle/180*pi
  
  #Rotate so that it matches Drasdo assumptions (Fovea-Optic Nerve angle is 0)
  R <- rbind(c(cos(angle), -sin(angle)), c(sin(angle), cos(angle)))
  XY <- R%*%rbind(x, y)
  x <- XY[1,]
  y <- XY[2,]
  
  #Transforms from degrees to mm (using schematic eye from Montesano et al. 2020, axial length = 24 mm)
  xmm <- sign(x)*approx(x = Degs, y = MM, abs(x))$y
  ymm <- sign(y)*approx(x = Degs, y = MM, abs(y))$y
  
  #Calculates displacement in mm using precalculated two dimensional LUTs
  xdispMM <- interp2(x = Drasdo_LUT$Xv, y = Drasdo_LUT$Yv, Z = Drasdo_LUT$xlut, xp = xmm, yp = ymm)
  ydispMM <- interp2(x = Drasdo_LUT$Xv, y = Drasdo_LUT$Yv, Z = Drasdo_LUT$ylut, xp = xmm, yp = ymm)
  
  #Unchanged outside diplacement zone (no NAs within 30 degrees)
  isNa <- which(is.na(xdispMM))
  xdispMM[isNa] <- xmm[isNa]
  ydispMM[isNa] <- ymm[isNa]
  
  #Coverts back to degrees using the schematic eye
  xdispDeg <- sign(xdispMM)*approx(x = MM, y = Degs, abs(xdispMM))$y
  ydispDeg <- sign(ydispMM)*approx(x = MM, y = Degs, abs(ydispMM))$y
  
  #Rerotate to match original input
  XY <- inv(R)%*%rbind(x, y)
  x <- XY[1,]
  y <- XY[2,]
  
  XY <- inv(R)%*%rbind(xdispDeg, ydispDeg)
  xdispDeg <- XY[1,]
  ydispDeg <- XY[2,]
  
  #Inverts y axis to convert back to visual field coordinates
  y <- -y
  ydispDeg <- -ydispDeg
  
  #Stores results for output
  coord <- data.frame(x = xdispDeg, y = ydispDeg)
  
  return(coord)
}
