
if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())


GSE_id="GSE136271"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

Platforms="GPL11533"  # Please check if there is the gene symbol.


LXgeo(GSE_id,Platforms)

