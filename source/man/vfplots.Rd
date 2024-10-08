% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vfplots.R
\name{vfgpar}
\alias{vfgpar}
\alias{vftess}
\alias{vfcolscheme}
\alias{vfprogcolscheme}
\alias{vfplot}
\alias{vfplotplr}
\alias{vflegoplot}
\alias{vfsparklines}
\title{Plots for visual fields data}
\usage{
vfgpar(
  coord,
  tess = vftess(coord),
  probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
  cols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
    colorRampPalette(c("#00FF00", "#008000"))(4)),
  floor = 0,
  ltprobs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 1),
  ltcols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
    "#008000"),
  gtprobs = c(0, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
  gtcols = c("#000000", "#FF0000", "#F7F0EB", colorRampPalette(c("#00FF00",
    "#008000"))(4)),
  neprobs = c(0, 0.0025, 0.005, 0.01, 0.25, 0.975, 0.99, 0.995, 0.9975, 1),
  necols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
    colorRampPalette(c("#FFFF00", "#FF0000"))(4)),
  bprobs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
  bcols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
    colorRampPalette(c("#00FF00", "#008000"))(4))
)

vftess(coord, floor = 0, delta = 3)

vfcolscheme(
  probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995, 1),
  cols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB",
    colorRampPalette(c("#00FF00", "#008000"))(4)),
  floor = 0
)

vfprogcolscheme(
  probs = c(0, 0.005, 0.01, 0.02, 0.05, 0.95, 1),
  cols = c("#000000", colorRampPalette(c("#FF0000", "#FFFF00"))(4), "#F7F0EB", "#008000")
)

vfplot(vf, td = NULL, tdp = NULL, pd = NULL, pdp = NULL, type = "td", ...)

vfplotplr(
  vf,
  alternative = "LT",
  xoffs = 0,
  yoffs = 0,
  addSpark = FALSE,
  thr = 2,
  width = 4,
  height = 2,
  ...
)

vflegoplot(
  vf,
  type = "td",
  grp = 3,
  addSpark = FALSE,
  thr = 2,
  width = 4,
  height = 2,
  ...
)

vfsparklines(vf, thr = 2, width = 4, height = 2, add = FALSE, ...)
}
\arguments{
\item{coord}{print x and y coordinates. Check section
\code{Structure of graphical parameters} for details}

\item{tess}{tesselation for the visual field maps. Check section
\code{Tesselation in visualFields} for details}

\item{probs}{probability scale to use for TD and PD values. It is a numeric vector
of probabilities with values in [0,1]. The values 0 and 1 must be included.
Although not technically necessary, it would be best if it is the same as for
the normative values used}

\item{cols}{corresponding colors for each of the probability levels}

\item{floor}{Flooring value, typically in dB. Default is 0}

\item{ltprobs, ltcols}{color map for progression with the alternative hypothesis
lower than (LT)}

\item{gtprobs, gtcols}{color map for progression with the alternative hypothesis
lower than (GT)}

\item{neprobs, necols}{color map for progression with the alternative hypothesis
not equal (NE)}

\item{bprobs, bcols}{color map for progression with blth alternative
hypotheses LT and GT (B for both)}

\item{delta}{Distance over which the boundary should be shifted. See for \code{polyclip}}

\item{vf}{the visual fields data to plot}

\item{td}{the total deviation values. If \code{NULL} (default) then use
visualFields normative values}

\item{tdp}{the total deviation probability values. If \code{NULL} (default)
then use visualFields normative values}

\item{pd}{the pattern deviation values. If \code{NULL} (default) then use
visualFields normative values}

\item{pdp}{the pattern deviation probability values. If \code{NULL} (default)
then use visualFields normative values}

\item{type}{the type of data to plot: sensitivities (`\code{s}`),
total deviation values (`\code{td}`), pattern deviation values (`\code{pd}`),
a hybrid plot that shows sensitivity grayscale with TD values and corresponding
probability levels (`\code{tds}`), or PD values and corresponding probability
levels (`\code{pds}`). Default is `\code{td}`.}

\item{...}{other graphical arguments. See \code{\link{plot}}}

\item{alternative}{alternative hypothesis used in progression analyses.
Allowed values are `\code{LT}` (as in "lower than", default),
`\code{GT}` (as in "greater than"), `\code{NE}` (as in "not equal"),
and `\code{both}` (both `\code{LT}` and `\code{GT}`)}

\item{xoffs, yoffs}{offset x and y where to print the slope values. That is,
the distance from the center of each Voronoy polygons in degrees of visual angle}

\item{addSpark}{whether to overlay a sparkline graph in each visual field location.
The parameters \code{thr}, \code{width}, and \code{height} are used only if
\code{addSpark} is \code{TRUE}. Default value is \code{FALSE}.}

\item{thr}{threshold used for the median absolute deviation of residuals
from simple linear regression. If greater than the threshold, the
sparkline for that location is plotted in red and with a thicker line.
Default is `\code{2}` (dB)}

\item{width}{the width of each pointwise sparkline plot. Default is
`\code{4}` (degrees of visual angle)}

\item{height}{the height of each pointwise sparkline plot. Default is
`\code{2}` (degrees of visual angle)}

\item{grp}{number of baseline (first) and last visual fields to group.
Default is `\code{3}`}

\item{add}{whether to generate a new plot (`\code{FALSE}`, as default)
or to add to an existing figure (`\code{TRUE}`)}
}
\value{
\code{vfgpar} returns a list with graphical parameters to be used for vfplots

\code{vftess} returns a list with the \code{xlim}, \code{ylim}, tessellation tiles and an outer hull
to be used for vfplots

\code{vfcolscheme} returns a list with a lookup table and a function that define the color scheme
to be used for vfplots

\code{vfprogcolscheme} returns the default \code{vfcolscheme} to be used for vfplots

\code{vfplot} No return value

\code{vfplotplr} No return value

\code{vflegoplot} No return value

\code{vfsparklines} No return value
}
\description{
Graphical tools for visualization and statistical analysis of
visual fields.
}
\details{
The following functions generate plots using visual fields data
\itemize{
  \item\code{vfgpar} generates simple graphical parameters
  \item\code{vftess} generates a structure to handle the visual field tessellation.
    Check section \code{Tesselation in visualFields} below for further details
  \item\code{vfcolscheme} generates the structures to handle the color scheme
    Check section \code{Color schemes in visualFields} below for further details
  \item\code{vfprogcolscheme} generates the structures to handle the color scheme
    for progression analysis. Check section \code{Color schemes in visualFields}
    below for further details
  \item\code{vfplot} plots a single test for visual field data
  \item\code{vfplotsens} plots a single test for visual field sensitivity data
    with a grayscale where darker means greater sensitivity loss
  \item\code{vfplotdev} plots a single test for visual field total or pattern
    deviation data with probability scales represented in color
  \item\code{vfplotplr} plots the results of pointwise linear regression for
    a series of visual fields for an eye from a subject
  \item\code{vflegoplot} the legoplot shows the differences between the average
    values of visual field tests taken as baseline and those at the end of
    follow up
  \item\code{vflegoplotsens} the legoplot for visual field sensitivity data with
    a grayscale where darker means greater sensitivity loss
  \item\code{vflegoplotdev} the legoplot for visual field total or pattern
    deviation data with probability scales represented in color
  \item\code{vfsparklines} the sparklines graph shows spark lines for the series
    of visual field sensitivities, or total or pattern deviation data for each
    location
}
}
\section{Structure of graphical parameters}{

Graphical parameters for visualFields must be a list containing
\itemize{
  \item\code{coord} print x and y coordinates. They could be different from the
    the real visual field location testing coordinates in complex visual field
    grids to help readability and improve visualization of statistical results
  \item\code{tess} tesselation for the visual field maps. Check section
    \code{Tesselation in visualFields}
  \item\code{colmap} color map representing the probability scale. Check section
    \code{Color schemes in visualFields}
}
A default graphical parameters can be generated with \code{generategpar}
}

\section{Tesselation in visualFields}{

A tesselation in visualFields must be defined with a list containing
\itemize{
  \item\code{xlim}, \item\code{ylim} 2-dimensional vectors containing the minimum
    and maximum x and y values
  \item\code{floor} the value to be assinged to any sensitivity value lower than
    \code{floor}
  \item\code{tiles} a list of as many tiles defining the tesselation as visual field
    test locations. Each element of the list is a table with x and y coordinates defining
    a polygon containing the corresponding test location. Each polygon is thus the tile
    for each visual field test location
  \item\code{hull} a table with x and y coordinates defining the outer hull of the
    tessellation
}
A default tessellation can be generated with \code{vftess}
}

\section{Color schemes in visualFields}{

A color scheme in visualFields must be defined with a list containing
\itemize{
  \item\code{map} a table mapping probabilities levels with colors defined
    in hexadecimal base
  \item\code{fun} a function that takes sensitivity values and deviation
    probability levels and returns the corresponding color code.
}
A default color scheme can be generated with \code{vfcolscheme}
}

\examples{
# generate a structure with default graphical parameters for the 30-2 map
vfgpar(locmaps$p30d2$coord)
# generate a structure with default tesselation for the 30-2 map
vftess(locmaps$p30d2$coord)
# default color scheme
vfcolscheme()
# default color scheme for progression
vfprogcolscheme()
# plot visual field values for the last field in the series for the first
# subject in the dataset vfpwgSunyiu24d2
# grayscale with sensitivity values
vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "s")
# TD values
vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "td")
# PD values
vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "pd")
# hybrid sensitivities and TD values
vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "tds")
# hybrid sensitivities and PD values
vfplot(vfselect(vffilter(vfpwgRetest24d2, id == 1), n = 1), type = "pds")
# plot results from pointwise linear regression for the series of
# visual fields for the right eye in the dataset vfpwgSunyiu24d2
# with sensitivity values
vfplotplr(vffilter(vfpwgSunyiu24d2, eye == "OD"))
# TD values
vfplotplr(gettd(vffilter(vfpwgSunyiu24d2, eye == "OD")))
# PD values
vfplotplr(getpd(gettd(vffilter(vfpwgSunyiu24d2, eye == "OD"))))
# legoplot for the series of visual fields for the right eye
# of the subject in the dataset vfpwgSunyiu24d2
# with sensitivity values
vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "s")
# TD values
vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "td")
# PD values
vflegoplot(vffilter(vfpwgSunyiu24d2, eye == "OD"), type = "pd")
# sparklines for the series of visual fields for the right eye of
# the subject in the dataset vfpwgSunyiu24d2
# with sensitivity values
vfsparklines(vffilter(vfpwgSunyiu24d2, eye == "OD"))
# TD values
vfsparklines(gettd(vffilter(vfpwgSunyiu24d2, eye == "OD")))
# PD values
vfsparklines(getpd(gettd(vffilter(vfpwgSunyiu24d2, eye == "OD"))))
}
