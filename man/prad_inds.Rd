\name{prad_inds}
\alias{prad_inds}
\docType{data}
\title{normal (0) and tumor (1) classes associated with PRAD expression data}
\description{
 normal (0) and tumor (1) classes associated with PRAD expression data.
}
\source{
\url{https://www.bioconductor.org/packages/release/bioc/html/RTCGA.html}
}
\references{
Kosinski M, Biecek P (2021). RTCGA: The Cancer Genome Atlas Data Integration. R package version 1.22.0, https://rtcga.github.io/RTCGA.}
\format{A character vector contains the class label}
\usage{data("PRAD_data")}
\examples{
data("PRAD_data")
table(prad_inds)
}
\keyword{datasets}
