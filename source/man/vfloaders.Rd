% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vfloaders.R
\name{loadhfaxml}
\alias{loadhfaxml}
\alias{loadhfadicom}
\alias{loadoctopus}
\alias{loadhfaxmlbatch}
\alias{loadhfadicombatch}
\title{Loaders from perimeters}
\usage{
loadhfaxml(file, type = "pwg", repeated = mean)

loadhfadicom(file, type = "pwg", repeated = mean)

loadoctopus(file, type = "pwg", repeated = mean, dateFormat = "\%d.\%m.\%Y")

loadhfaxmlbatch(file, repeated = mean)

loadhfadicombatch(file, repeated = mean)
}
\arguments{
\item{file}{name of the csv file exported by the eyesuite software}

\item{type}{type of patient. It can be `\code{ctr}` (for control or healthy
subject-eye) or `\code{pwg}` (for patient with glaucoma) or other}

\item{repeated}{function to apply if there are repeated values in a particular location}

\item{dateFormat}{format to be used for date. Its default value is \%d.\%m.\%Y}
}
\value{
Visual field data
}
\description{
Functions to load from commercial perimeters
}
\details{
The XML loader for the Humphrery Field Analyser (HFA) by Carl Zeiss Meditec
is essentially a XML parser that reads in the XML generated with the scientific
export license. The DICOMM loader is also a parser to read HFA data generated in a
DICOMM file. The loader for the Octopus perimeter by Haag-Streit is a csv reader
from files generated with the Eyesuite software. The parser also extracts information
on visual field pattern deviation values and normative values. The list that is returned
with the \code{loadoctopus} loader contains data frames which are structured with keys so that
redundancy is minimized (similar to a relational database). Detailed examples for
\code{loadoctopus}: \url{https://huchzi.github.io/vf1/loadOctopus.html}
}
