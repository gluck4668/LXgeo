\name{LXgeo::LXgeo}
\alias{LXgeo}
\title{LXgeo::LXgeo}
\usage{
LXgeo(GSE_id,Platforms)
}
\description{
LXgeo can be used to download GEO data with gene symbol.
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())

setwd("D:/Desktop/R_example/LXgeo_example")

GSE_id="GSE133969"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

Platforms="GPL11180" # Please check if there is the gene symbol.


LXgeo(GSE_id,Platforms)

}
