\name{LXgeo}
\alias{LXgeo}
\title{LXgeo}
\usage{
LXgeo(GSE_id,Platforms)
}
\description{
LXgeo can be used to get the GEO gene expression data.

}
\examples{
if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())

GSE_id="GSE136271"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
Platforms="GPL11533"

LXgeo(GSE_id,Platforms)

}
