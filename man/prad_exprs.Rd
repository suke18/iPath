\name{prad_exprs}
\alias{prad_exprs}
\docType{data}
\title{expression matrix for PRAD cancer patients in TCGA}
\source{
\url{https://www.bioconductor.org/packages/release/bioc/html/RTCGA.html}
}
\description{
prad_exprs is the RPKM expression matrix which belongs to "matrix" class. The data includes 102 samples about human preimplantation embryos and embryonic stem cells. It contains 19304 genes after removing genes with extreme high dropout rate.
}
\references{
Kosinski M, Biecek P (2021). RTCGA: The Cancer Genome Atlas Data Integration. R package version 1.22.0, https://rtcga.github.io/RTCGA.}
\format{An object of "matrix" class contains the mRNA expressions}
\usage{data("PRAD_data")}
\examples{
data("PRAD_data")
prad_exprs[1:10, 1:4]
}
\keyword{datasets}
