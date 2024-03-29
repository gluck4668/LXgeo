
LXgeo <- function(GSE_id,Platforms){

  #--------------------------------------------

  # options ( "repos "="https://mirrors.ustc.edu.cn/CRAN/")
  # options(Bioc_mirror="https://mirrors.ustc.edu.cn/bioc/")


  all_packages <- installed.packages()[,1]
  
  com_pack <- c("BiocManager","openxlsx","psych","dplyr","tidyverse","stringr",
                "pacman","umap","nortest","car","Cairo","BiocGenerics","tidyr",
                "ggplot2","plyr","psych","ggrepel",'RCurl',
                "patchwork","raster","png","janitor")
  
  bioc_pack <- c("GEOquery","limma","affy")
  
  not_com <- com_pack[!com_pack %in% all_packages]
  not_bioc <- bioc_pack[!bioc_pack %in% all_packages]
  
  if (length(not_com)>0) {
    com_fun <- function(x){install.packages(x)}
    sapply(not_com,com_fun,simplify = T)
  }
  
  if (length(not_bioc)>0){
    bioc_fun <- function(x){BiocManager::install(x)}
    sapply(not_bioc,bioc_fun,simplify = T)
  }
  
  lib_fun <- function(x){library(x,character.only = T)}
  sapply(c(com_pack,bioc_pack),lib_fun)

  rm(all_packages,com_pack,bioc_pack,not_com,not_bioc,lib_fun)
  
  #-------------------------------------------------------------------------
  if(!dir.exists("analysis results"))
    dir.create("analysis results")

  #-------- 下载 GSE数据矩阵---------------
  GSE_id <- trimws(GSE_id)
  Platforms <- trimws(Platforms)

  options (stringsAsFactors = F)#不让字符串转化为因子
  Sys.setenv("VROOM_CONNETION_SIZE"=19999999999)#为内存不足，运行代码
  options('download.file.method.GEOquery'='libcurl') # 下载设置
  options('GEOquery.inmemory.gpl'=FALSE)
  options(timeout=6000000)
  getOption("timeout")

  gse <- getGEO (GSE_id, destdir = "analysis results/",GSEMatrix =T,
                # getGPL = F,
                AnnotGPL= T
                )

#----查看platforms ID-----------------------------------
  gse_names <- attr(gse, "names")
  
  plat_fun <- function(x){gse[[x]]@annotation}
  
  sapply(gse_names,plat_fun)
  
  Plats <- sapply(gse_names,plat_fun) %>% as.character()
 
 
 if (length(gse) > 1)
    idx <- grep(Platforms, attr(gse, "names"),ignore.case = T) else
      idx <- 1

  gse_df <- exprs(gse[[idx]]) %>% data.frame() # 表达矩阵

#-----判断是否有gene表达矩阵----------------------------------------------     
if(nrow(gse_df)>0)
  print( paste(GSE_id, "provided a gene expression series matrix.")) else 
    stop(paste0("There is no gene expression data in ",GSE_id,"_series_matrix",". Please check it." ))

#-----------------------------------------------------------------------------

  gse_col <- colnames(gse_df) %>% trimws()
  group_df <- pData(gse[[idx]]) # 实验分组情况
  write.xlsx(group_df,paste0("analysis results/",GSE_id," group information.xlsx"))

  for(i in 1:length(gse_col)){
       if(gse_col[i]== group_df[i,2])
          gse_col[i]<- group_df[i,1]
        }

  gse_col

  colnames(gse_df) <- gse_col

  gse_df$ID <- rownames(gse_df)
  

  #-------下载实验平台GPL数据---------------------------------------------
  #-----载入上述获得的annot.gz数据（AnnotGPL= T）

  gpl_labels <- fvarLabels(gse[[idx]]) %>% tolower() # 查看gpl包含的项目

  if(TRUE %in% grepl("symbol",gpl_labels))
    print(paste(Platforms," provided the gene symbol.")) else
    {symbol_txt <- paste(Platforms,"did not provide the gene symbol.")
    stop(symbol_txt)}

  gpl_ann_file <- paste0("analysis results/",Platforms,".annot.gz")

  if(file.exists(gpl_ann_file)){
  gpl_ann_df <- getGEO(filename=gpl_ann_file)
  gpl_ann_df =gpl_ann_df@dataTable@table

  gpl_ann_col <- colnames(gpl_ann_df) %>% tolower()
  ann_id  <- grep("id",gpl_ann_col) %>% as.numeric()
  ann_symbol <- grep("symbol",gpl_ann_col) %>% as.numeric()

  gpl_ann_df <- gpl_ann_df[,c(ann_id,ann_symbol)]
  colnames(gpl_ann_df)[1] <- "ID"

  gpl_df <- gpl_ann_df} else

  #----当然，也可以直接下载soft.gz注释文件---------
  {gpl_soft <- getGEO(Platforms, destdir = "analysis results/")
  gpl_soft_file <- paste0("analysis results/",Platforms,".soft.gz")

  gpl_soft_df <- getGEO(filename=gpl_soft_file)
  gpl_soft_df =gpl_soft@dataTable@table

  gpl_soft_col <- colnames(gpl_soft_df) %>% tolower()
  soft_id  <- grep("id",gpl_soft_col) %>% as.numeric()
  soft_symbol <- grep("symbol",gpl_soft_col) %>% as.numeric()

  gpl_soft_df <- gpl_soft_df[,c(soft_id,soft_symbol)]
  colnames(gpl_soft_df)[1] <- "ID"

  gpl_df <- gpl_soft_df}
  
  
  #------------合并ges和gpl数据，以ID为依据------------
  gpl_df$ID <- as.character(gpl_df$ID)
  gse_df$ID <- as.character(gse_df$ID)

  GEO_df <- dplyr::inner_join(gpl_df,gse_df,"ID")

  if(TRUE %in% grepl("uni",tolower(colnames(GEO_df)))){
      GEO_uni <- grep("unigene",tolower(colnames(GEO_df))) %>% as.numeric()
      GEO_df <- GEO_df[,-c(GEO_uni)]}

  GEO_id <- grep("id",tolower(colnames(GEO_df))) %>% as.numeric()
  GEO_df <- GEO_df[,-c(GEO_id)]
  colnames(GEO_df)[1] <- "gene_symbol"
  GEO_df <- tidyr::separate(data = GEO_df,col = gene_symbol,into = "gene_symbol",sep = "/",remove = T )

  write.xlsx(GEO_df, paste0("analysis results/",GSE_id,"-",Platforms," gene expression matrix.xlsx"))


  #----------- limma 差异基因分析----------------------

  LXlimma <- function(gene_matrix){

    exp_m <- gene_matrix

    # 去掉NA
    exp_m <- na.omit(exp_m)

    # 去掉symbol重复的行，并取重复行的平均值
    exp_m <- limma::avereps(exp_m[,-1],ID=exp_m[,1]) %>% data.frame()

    # log2 transformation
    qx <- as.numeric(quantile(exp_m, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||  (qx[6]-qx[1] > 50 && qx[2] > 0)

    exp_m2 <- log2(exp_m)

    if (LogC){exp_m[which(exp_m <= 0)] <- NaN
    exp_m2 <- log2(exp_m)}

    # 组别数据
    group_list <- colnames(exp_m2)
    groups <- gsub("\\d+$", "", group_list) # 去掉字符串后面的数字


    #构建design 样本分组矩阵
    design <- model.matrix(~ 0 + groups)
    rownames(design) = group_list
    colnames(design) <- levels(factor(groups))

    # 明确 model vs normal
    #cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
    #contrast.matrix <- makeContrasts(cts, levels = design)
    #contrast.matrix

    # 进行差异分析：
    #fit <- lmFit(exp_m, design)
    #fit <- contrasts.fit(fit, contrast.matrix)
    #fit <- eBayes(fit)
    #allDiff <- topTable(fit, number = Inf)


    fit <- lmFit(exp_m2, design)  # fit linear model

    # set up contrasts of interest and recalculate model coefficients
    cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)

    # compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by="B",number=Inf)
    tT$gene_symbol <- rownames(tT)

    exp_m$gene_symbol <- rownames(exp_m)
    exp_m <- exp_m[,c(ncol(exp_m),1:ncol(exp_m)-1)]

    tT <- dplyr::inner_join(exp_m,tT, "gene_symbol")

    file_save <- paste0("analysis results/",GSE_id,"-",Platforms," statistical data (",cts,").xlsx")

    FC_n <- grep("logfc",tolower( colnames(tT) )) %>% as.numeric()

    colnames(tT)[FC_n] <- c("log2FC")

    write.xlsx(tT, file_save, rowNames=F)

    #-----------差异基因：火山图-----------------

       vol_df <-tT

       vol_df$type <- case_when(vol_df$log2FC>0 & vol_df$P.Value<0.05 ~"Up",
                                vol_df$log2FC<0 & vol_df$P.Value<0.05 ~"Down",
                                TRUE~ "Not sig"
                                )

      change_n <- table(vol_df$type)
      change_nm <-data.frame(change_n)

      Down <- c(paste("Down",c(change_nm[1,2])))#Down
      NC <- c(paste("Not sig",c(change_nm[2,2])))#Not sig
      UP <- c(paste("Up",c(change_nm[3,2])))#Up

      vol_df$changes <- case_when(vol_df$type=="Up" ~UP,
                                  vol_df$type=="Down" ~Down,
                                  TRUE ~ NC
                                 )

      mytheme<-theme_bw()+
        theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
              #panel.border = element_rect (linewidth = 0.8,color = "gray30"),
              axis.line = element_blank(),
              axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
              axis.ticks.length = unit(1.5,units = "mm"))+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))

      xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
        theme(axis.text.y = element_text(face="bold",color="black",size=12))+
        theme(legend.text=element_text(face="bold",color="black",size=12))

      p1 <- ggplot(vol_df,aes(x=log2FC,y=-log10(P.Value)))+
        geom_point(aes(color=changes))+
        scale_color_manual(values = c("#00bfff","#c0c0c0","#ff4500"))+
        geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
        geom_vline(xintercept = c(-1,1),linetype="dashed",color="#808080")+
        labs(x="log2FC",y="-log10(p value)",title = "Vocano diagram") +
        labs(color="Types")+  #修改legned标题
        mytheme+xytheme

      p1


      # group_names <- paste0("(",group1," VS ",group2,")")

      vol_file <- paste0("analysis results/",GSE_id,"-",Platforms," DEGs data (",cts,").xlsx")

      DEGs <- dplyr::filter(vol_df,type!="Not sig")


      DEGs_file <- DEGs[,-ncol(DEGs)]


      write.xlsx(DEGs_file,vol_file)

      p1_file <- paste0("analysis results/",GSE_id,"-",Platforms," Gene Volcano graphics (",cts,").png")

      ggsave(p1_file,p1, width=1200, height =1000, dpi=180,units = "px")


  }

  #-----------判断分组情况--------------------

  col_name <- colnames(GEO_df) # 查看组别
  col_list <- gsub("\\d$","",col_name[-1]) # 去除字符串后的数字
  group_type <- col_list[!duplicated(col_list)] # 字符串去重
  group_num <- length(group_type) %>% as.numeric() # 组数

  if(group_num==2){
     gene_matrix = GEO_df
     LXlimma(gene_matrix) }

  if(group_num>=3){

    gene_matrix <- GEO_df[,c(1,
                        grep(group_type[1],col_name),
                        grep(group_type[2],col_name))]
    LXlimma(gene_matrix)


    for(i in 1:(group_num-2)){
      gene_matrix <- GEO_df[,c(1,
                          grep(group_type[2],col_name),
                          grep(group_type[i+2],col_name))]

      LXlimma(gene_matrix)

      }

    }

 #---------------------------------------------------------

print("The result could be found in the folder of <analysis results>")
  
}


#----------------自定义比较组别----------------------------

LXgeo_02 <- function(group1,group2){  


LXlimma_02 <- function(gene_matrix){
  
  exp_m <- gene_matrix
  
  # 去掉NA
  exp_m <- na.omit(exp_m)
  
  # 去掉symbol重复的行，并取重复行的平均值
  exp_m <- limma::avereps(exp_m[,-1],ID=exp_m[,1]) %>% data.frame()
  
  # log2 transformation
  qx <- as.numeric(quantile(exp_m, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||  (qx[6]-qx[1] > 50 && qx[2] > 0)
  
  exp_m2 <- log2(exp_m)
  
  if (LogC){exp_m[which(exp_m <= 0)] <- NaN
  exp_m2 <- log2(exp_m)}
  
  # 组别数据
  group_list <- colnames(exp_m2)
  groups <- gsub("\\d+$", "", group_list) # 去掉字符串后面的数字
  
  
  #构建design 样本分组矩阵
  design <- model.matrix(~ 0 + groups)
  rownames(design) = group_list
  colnames(design) <- levels(factor(groups))
  
  # 明确 model vs normal
  #cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
  #contrast.matrix <- makeContrasts(cts, levels = design)
  #contrast.matrix
  
  # 进行差异分析：
  #fit <- lmFit(exp_m, design)
  #fit <- contrasts.fit(fit, contrast.matrix)
  #fit <- eBayes(fit)
  #allDiff <- topTable(fit, number = Inf)
  
  
  fit <- lmFit(exp_m2, design)  # fit linear model
  
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(levels(factor(groups))[1], levels(factor(groups))[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B",number=Inf)
  tT$gene_symbol <- rownames(tT)
  
  exp_m$gene_symbol <- rownames(exp_m)
  exp_m <- exp_m[,c(ncol(exp_m),1:ncol(exp_m)-1)]
  
  tT <- dplyr::inner_join(exp_m,tT, "gene_symbol")
  
  if(!dir.exists("analysis results 02"))
      dir.create("analysis results 02")
  
  file_save <- paste0("analysis results 02/",GSE_id,"-",Platforms," statistical data (",cts,").xlsx")
  
  FC_n <- grep("logfc",tolower( colnames(tT) )) %>% as.numeric()
  
  colnames(tT)[FC_n] <- c("log2FC")
  
  write.xlsx(tT, file_save, rowNames=F)
  
  #-----------差异基因：火山图-----------------
  
  vol_df <-tT
  
  vol_df$type <- case_when(vol_df$log2FC>0 & vol_df$P.Value<0.05 ~"Up",
                           vol_df$log2FC<0 & vol_df$P.Value<0.05 ~"Down",
                           TRUE~ "Not sig"
  )
  
  change_n <- table(vol_df$type)
  change_nm <-data.frame(change_n)
  
  Down <- c(paste("Down",c(change_nm[1,2])))#Down
  NC <- c(paste("Not sig",c(change_nm[2,2])))#Not sig
  UP <- c(paste("Up",c(change_nm[3,2])))#Up
  
  vol_df$changes <- case_when(vol_df$type=="Up" ~UP,
                              vol_df$type=="Down" ~Down,
                              TRUE ~ NC
  )
  
  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
          #panel.border = element_rect (linewidth = 0.8,color = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))
  
  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=12))+
    theme(legend.text=element_text(face="bold",color="black",size=12))
  
  p1 <- ggplot(vol_df,aes(x=log2FC,y=-log10(P.Value)))+
    geom_point(aes(color=changes))+
    scale_color_manual(values = c("#00bfff","#c0c0c0","#ff4500"))+
    geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#808080")+
    geom_vline(xintercept = c(-1,1),linetype="dashed",color="#808080")+
    labs(x="log2FC",y="-log10(p value)",title = "Vocano diagram") +
    labs(color="Types")+  #修改legned标题
    mytheme+xytheme
  
  p1
  
  
  # group_names <- paste0("(",group1," VS ",group2,")")
  
  vol_file <- paste0("analysis results 02/",GSE_id,"-",Platforms," DEGs data (",cts,")_Definition.xlsx")
  
  DEGs <- dplyr::filter(vol_df,type!="Not sig")
  
  
  DEGs_file <- DEGs[,-ncol(DEGs)]
  
  
  write.xlsx(DEGs_file,vol_file)
  
  p1_file <- paste0("analysis results 02/",GSE_id,"-",Platforms," Gene Volcano graphics (",cts,")_Definition.png")
  
  ggsave(p1_file,p1, width=1200, height =1000, dpi=180,units = "px")
  
  
}

#-----------判断分组情况--------------------

group_n1 <- gsub("\\d$","",group1) %>% tolower()
group_n2 <- gsub("\\d$","",group2) %>% tolower()

GEO_df <- read.xlsx(paste0("analysis results/",GSE_id,"-",Platforms," gene expression matrix.xlsx"))

colnames(GEO_df)[2:ncol(GEO_df)] <- group_information[,1]
col_name <- colnames(GEO_df)
col_name <- gsub("\\d$","",col_name) %>% tolower()

gene_matrix <- GEO_df[,c(1,
                           grep(group_n1,col_name),
                           grep(group_n2,col_name))]
LXlimma_02(gene_matrix)
  

#---------------------------------------------------------

print("The result could be found in the folder of <analysis results 02>")

}




#--------下载 supplementary file---------------------------------------

download_supplementary_file <- function(GSE_id){
  
dir_name <- paste0("Supplementary file (",GSE_id,")")

if (!dir.exists(dir_name))
  dir.create(dir_name)

website=paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",substr(GSE_id,1,6),"nnn/",GSE_id,"/suppl/")

html=getURL(website) %>% as.character()

#--截取以 GSE开头，<结尾的字符串----
files <- str_extract_all(html,"(?<=>GSE).*(?=<)") %>% unlist() %>% paste0("GSE",.)

#---批下载文件-----------------    
down_fun <- function(x)
{url=paste0(website,x)
file_path=paste0(dir_name,"/",x)
if (!file.exists(file_path))
 download.file(url,file_path,method = "libcurl",quiet = F ,mode = "wb") else
   print (paste0(x," was existed, and thus it was not neccessary to be downloaded again."))
}

sapply(files,down_fun)

print(paste0("The supplementary files can be found in the folder of ",dir_name))  
}

