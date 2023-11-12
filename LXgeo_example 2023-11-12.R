

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgeo")

library(LXgeo)

??LXgeo

rm(list=ls())


GSE_id="GSE136271"  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

Platforms="GPL11533"  


#---GEO 基因表达矩阵下载及统计---------
LXgeo(GSE_id,Platforms)

#-----自定义组间比较------------------
group_information <- read.xlsx(paste0("analysis results/",GSE_id," group information.xlsx")) #查看分组信息

group1= "Liver, cis-unsaturated fat, replicate 1"

group2= "Liver, saturated fat, replicate 1"   

LXgeo_02(group1,group2)  

#-------下载 Supplementary files------

download_supplementary_file (GSE_id)


