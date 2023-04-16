
if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())


GSE_id="GSE222576"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

Platforms="GPL22598"  

options('download.file.method.GEOquery'='libcurl')
options('GEOquery.inmemory.gpl'=FALSE)

LXgeo(GSE_id,Platforms)


#-----自定义组间比较-----------
group_information <- read.xlsx(paste0("analysis results/",GSE_id," group information.xlsx")) #查看分组信息

 group1= "Control at 10 weeks, biological rep"
  
 group2= "CCl4-induced liver fibrosis for 10 weeks, biological rep"   
  
 LXgeo_02(group1,group2)  



