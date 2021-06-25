\name{prad_cli}
\alias{prad_cli}
\docType{data}
\title{simulated clinical data for PRAD cancer patients}
\source{
\url{https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga}
}
\description{
prad_cli is the clinical data containing three variables times, bcr_patient_barcode, and patient.vital_status.
}
\references{
Kosinski M, Biecek P (2021). RTCGA: The Cancer Genome Atlas Data Integration. R package version 1.22.0, https://rtcga.github.io/RTCGA.}
\format{An object of "matrix" class contains the clinical outcomes}
\usage{data("PRAD_data")}
\examples{
data("PRAD_data")
prad_cli[1:10,]
}
\keyword{datasets}
