
if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())

setwd("D:\Desktop\R_example\LXgeo_example")

GSE_id="GSE133969"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

Platforms="GPL11180"


LXgeo(GSE_id,Platforms)


#--------------------------------------------
