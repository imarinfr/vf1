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
#' #' N. Drasdo, C. L. Millican, C. R. Katholi, and C. A. Curcio. \emph{The
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
#' @export
vf2gc <- function(coord) {
  ###########################################################################
  # below is what IMF arrived at
  ###########################################################################
  #ecc   <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7)
  #displ <- c(0.00000000, 0.07852941, 0.15705882, 0.23558824, 0.31411765, 0.39264706, 0.46583333, 0.53472222, 0.60333333, 0.66750000, 0.72333333, 0.79000000, 0.85416667, 0.92333333, 0.98625000, 1.04666667, 1.11000000, 1.16666667, 1.21388889, 1.26833333, 1.33583333, 1.39000000, 1.44000000, 1.48888889, 1.53333333, 1.57583333, 1.62333333, 1.67333333, 1.71600000, 1.75916667, 1.80166667, 1.84555556, 1.89222222, 1.93250000, 1.96500000, 1.99833333, 2.02969697, 2.05083333, 2.07000000, 2.08583333, 2.10000000, 2.11190476, 2.12000000, 2.12333333, 2.12000000, 2.11600000, 2.10333333, 2.08296296, 2.06047619, 2.03388889, 1.99666667, 1.94833333, 1.90642857, 1.85833333, 1.80666667, 1.75958333, 1.70958333, 1.65916667, 1.60888889, 1.56333333, 1.51833333, 1.47444444, 1.42500000, 1.38333333, 1.34333333, 1.30750000, 1.27333333, 1.24083333, 1.21416667, 1.18666667, 1.16500000, 1.14222222, 1.12285714, 1.10428571, 1.08703704, 1.06833333, 1.05250000, 1.03666667, 1.02666667, 1.01533333, 1.00333333, 0.99333333, 0.98055556, 0.97083333, 0.96333333, 0.95666667, 0.95000000, 0.94153846, 0.93666667, 0.92833333, 0.92333333, 0.92000000, 0.91666667, 0.91333333, 0.91000000, 0.91000000, 0.90666667, 0.90333333, 0.89666667, 0.89000000, 0.88333333, 0.87666667, 0.86500000, 0.85250000, 0.83666667, 0.82190476, 0.80666667, 0.78666667, 0.76939394, 0.75000000, 0.73250000, 0.71333333, 0.69333333, 0.67416667, 0.65222222, 0.63111111, 0.61190476, 0.59000000, 0.57000000, 0.55037037, 0.52857143, 0.50714286, 0.48333333, 0.46250000, 0.44583333, 0.42666667, 0.40727273, 0.38787879, 0.36848485, 0.34909091, 0.32969697, 0.31030303, 0.29090909, 0.27151515, 0.25212121, 0.23272727, 0.21333333, 0.19393939, 0.17454545, 0.15515152, 0.13575758, 0.11636364, 0.09696970, 0.07757576, 0.05818182, 0.03878788, 0.01939394, 0.00000000)
  ###########################################################################
  # below is what AT arrived at (manipulated by IMF for steps of 0.1 degrees)
  ###########################################################################
  ecc   <- c(0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2.0,
                 2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,  3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4.0,
                 4.1,  4.2,  4.3,  4.4,  4.5,  4.6,  4.7,  4.8,  4.9,  5.0,  5.1,  5.2,  5.3,  5.4,  5.5,  5.6,  5.7,  5.8,  5.9,  6.0,
                 6.1,  6.2,  6.3,  6.4,  6.5,  6.6,  6.7,  6.8,  6.9,  7.0,  7.1,  7.2,  7.3,  7.4,  7.5,  7.6,  7.7,  7.8,  7.9,  8.0,
                 8.1,  8.2,  8.3,  8.4,  8.5,  8.6,  8.7,  8.8,  8.9,  9.0,  9.1,  9.2,  9.3,  9.4,  9.5,  9.6,  9.7,  9.8,  9.9, 10.0,
                10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0,
                12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0,
                14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0,
                16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0,
                18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0,
                20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7)
  displ <- c(0, 0.442228682, 0.654866680, 0.813841218, 1.023691357, 1.230218115, 1.393827498, 1.524305697, 1.657796512, 1.801673490,
                1.936625202, 2.045040405, 2.151033197, 2.273349180, 2.394764622, 2.478420950, 2.552366923, 2.632627210, 2.718063589,
                2.805373056, 2.874372871, 2.919787185, 2.963937875, 3.005817878, 3.047145830, 3.088718658, 3.111870869, 3.128545044,
                3.149989589, 3.161397491, 3.168765031, 3.174186334, 3.175177600, 3.175726209, 3.174781269, 3.169121629, 3.159110664,
                3.147766429, 3.134890806, 3.118395254, 3.099389724, 3.079521492, 3.059851108, 3.039782323, 3.014900890, 2.985654085,
                2.956410295, 2.929108071, 2.897298913, 2.857283115, 2.817945867, 2.784503024, 2.752566867, 2.715462980, 2.681927757,
                2.647076423, 2.606348459, 2.573205530, 2.540634663, 2.504214893, 2.467552496, 2.435851194, 2.407590327, 2.378353997,
                2.347825526, 2.315430535, 2.281852049, 2.248860853, 2.217813863, 2.188741954, 2.159359652, 2.126327991, 2.094083302,
                2.065156685, 2.038932990, 2.016372228, 1.994265982, 1.969521197, 1.945333993, 1.920584422, 1.894810674, 1.869990210,
                1.848228647, 1.824579708, 1.796628904, 1.776220012, 1.758086432, 1.738404997, 1.719438227, 1.699101384, 1.680364015,
                1.666417045, 1.651485798, 1.637162540, 1.623777193, 1.608915732, 1.594572880, 1.580402003, 1.565281197, 1.552129428,
                1.541537470, 1.531941372, 1.521742507, 1.511825323, 1.502137966, 1.492549755, 1.485157805, 1.478010386, 1.469928266,
                1.462111456, 1.454565238, 1.447514381, 1.440792109, 1.432802047, 1.426345113, 1.420707712, 1.412065551, 1.405283323,
                1.400699387, 1.394067432, 1.385987937, 1.381063058, 1.378329714, 1.372071637, 1.366568791, 1.362580380, 1.359310336,
                1.358949246, 1.357754316, 1.354350052, 1.351932213, 1.347480452, 1.341900096, 1.337936922, 1.334995455, 1.331130914,
                1.323619694, 1.311114060, 1.298294826, 1.288326299, 1.280810509, 1.270317742, 1.258155791, 1.246139463, 1.232280317,
                1.217776716, 1.203527941, 1.189466066, 1.174816812, 1.161636013, 1.149192340, 1.133022468, 1.116843087, 1.101576790,
                1.085251368, 1.066767800, 1.048863628, 1.032664258, 1.015141010, 0.997719851, 0.980921838, 0.964151010, 0.947662406,
                0.930475227, 0.912789729, 0.896562356, 0.878892415, 0.860469892, 0.844115616, 0.828064397, 0.810796538, 0.792973401,
                0.776483900, 0.759757283, 0.742592800, 0.726043444, 0.709451154, 0.691134122, 0.671915198, 0.655932488, 0.641106075,
                0.625729083, 0.609800030, 0.593357664, 0.576440731, 0.559087979, 0.541338154, 0.523230004, 0.504802275, 0.486093715,
                0.467143071, 0.447989090, 0.428670518, 0.409226103, 0.389694591, 0.370114731, 0.350525268, 0.330964951, 0.311472525,
                0.292086738, 0.272846337, 0.253790070, 0.234956682, 0.216384922, 0.198113536, 0.180181271, 0.162626874, 0.145489093,
                0.128806674, 0.112618364, 0.096962910, 0.081879061, 0.067405561, 0.053581159, 0.040444602, 0.028034636, 0.016390009, 0)
  # convert to polar coordinates
  r     <- sqrt(coord[,1]^2 + coord[,2]^2)
  theta <- atan2(coord[,2], coord[,1])
  # get soma GC eccentricity
  r <- r + approx(ecc, displ, r, rule = 2)$y
  #convert back to cartesian coordinates
  coord[,1] <- r * cos(theta)
  coord[,2] <- r * sin(theta)
  return(coord)
}