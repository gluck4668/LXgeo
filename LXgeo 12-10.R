

LXgeo <- function(GSE_id,Platforms){

#--------------------------------------------
 options (stringsAsFactors = F)#不让字符串转化为因子
 Sys.setenv("VROOM_CONNETION_SIZE"=19999999999)#为内存不足，运行代码
 options('download.file.method.GEOquery'='libcurl')

 # options ( "repos "="https://mirrors.ustc.edu.cn/CRAN/")
 # options(Bioc_mirror="https://mirrors.ustc.edu.cn/bioc/")


 all_packages <- data.frame(installed.packages())

 pack <- data.frame(c("BiocManager","openxlsx","psych","dplyr","tidyverse","stringr",
                      "pacman"))

 bioc_pack <- data.frame(c("GEOquery","limma","affy"))

 pack$type <- pack[,1] %in% all_packages$Package

 for (i in 1:nrow(pack)){if (!requireNamespace(pack[i,1], quietly=TRUE))
                 install.packages(pack[i,1],update = F,ask = F)}
 rm(i)

 for (i in 1:nrow(bioc_pack)){if (!requireNamespace(bioc_pack[i,1], quietly=TRUE))
   BiocManager::install (bioc_pack[i,1],update = F,ask = F) }

 rm(i)

 packages <- c(pack[,1],bioc_pack[,1])

 for(i in packages){
   library(i, character.only = T)}

 rm(i)


#-------------------------------------------------

if(dir.exists("analysis results")==F)
  dir.create("analysis results")


#-------- GEO数据下载 ---------------

  # getGEO函数会加载GSE的matrix文件，默认会下载其注释探针信息，
  # 并对表达矩阵中的探针予以注释，但往往注释文件比较大，会出现parse保存的问题，
  # 所以一般建议把注释关掉了：getGPL=F，然后在后续分析步骤里进行手动注释。

 GSE_id <- trimws(GSE_id)
 Platforms <- trimws(Platforms)

 gse_file=paste0("analysis results/",GSE_id, ".Rdata")

 if(!file.exists(gse_file)){
    gse <- getGEO (GSE_id, destdir = "analysis results/", GSEMatrix =T,
                   getGPL = F, # 平台文件
                   AnnotGPL = F) # 注释文件

    save(gse,file=gse_file)
    }

    load(gse_file)


#-------- GE0数据整理 ---------------
 exp_matrix <- exprs(gse[[1]]) #获取表达矩阵

 if(nrow(exp_matrix)==0)
      {stop_text <- paste("There is no gene expression data in",GSE_id,"series matrix",". Please check it." )
      stop(stop_text)}

 group_infor <-pData(gse[[1]]) #获取分数等实验基本数据


#------- 实验基本信息----------------

 matrix_infor <- paste0(GSE_id,"_series_matrix.txt.gz")

 experiment_infor <-  gse[[matrix_infor]]@experimentData

 black <- c("                ")

 title <- paste("Title:",experiment_infor@title)

 abstract <- paste("Abstract:",experiment_infor@abstract)

 sample <- paste("Samples:",nrow(group_infor))

 group_name <- c("Group Name / ID:", paste(group_infor$title,"/",group_infor$geo_accession))

 # group_file <- c(sample,group_name)

 experiment_infor_file <- c(black,title,black,abstract,black,sample,black,group_name)

 write(x = experiment_infor_file,file = paste0("analysis results/",GSE_id,"_experiment information.txt"))

 sample_n <- as.character(nrow(group_infor))

 if(nrow(group_infor)<6)
     {sample_pr <- paste("The samples are lese than 6,",GSE_id,"may not be suitable")
      print(sample_pr)} else
          {sample_pr <- paste("The samples are",sample_n, ", which are adequate!")
           print(sample_pr)}


#------soft文件列名不统一，活学活用，有的表格里没有symbol列，也有的GPL平台没有提供注释

  gpl_soft = getGEO(Platforms, destdir = "analysis results/")
  gpl_id =gpl_soft@dataTable@table

  if(nrow(gpl_id)==0)
       {stop_text <- paste("There is no platform data in",Platforms,". Please check it." )
        stop(stop_text)}

  g_symbol <- colnames (gpl_id) %>% tolower() %>% data.frame()
  colnames(g_symbol) <- "Item"
  # g_symbol <- dplyr::filter(g_symbol,str_count(title)>=6)

  if("gene_symbol" %in% g_symbol[,1]==T | "gene symbol" %in% g_symbol[,1]==T | "symbol" %in% g_symbol[,1] ==T)
    {symbol_txt <- paste(GSE_id,"provided the gene symbol, which can be used.")
    print(symbol_txt)} else
         {symbol_txt <- paste(GSE_id,"did not provide the gene symbol.")
           stop(symbol_txt)}

  colnames(gpl_id) <- colnames (gpl_id) %>% tolower()

  if("gene_symbol" %in% g_symbol[,1]==T)
    gpl_df = gpl_id[, c("id","gene_symbol")] else
      {if("gene symbol" %in% g_symbol[,1]==T)
         gpl_df = gpl_id[, c("id","gene symbol")] else
            gpl_df = gpl_id[, c("id","symbol")]
       }

  colnames(gpl_df) <- c("ID","gene_symbol")

  gpl_df <- gpl_df[gpl_df$gene_symbol!="" & !str_detect(gpl_df$gene_symbol, "---"),]


#-------- 芯片id和gene symbol 进行慕匹配-------------------

  exp_df<-as.data.frame (exp_matrix)#数据框
  exp_df$ID<-rownames (exp_matrix)#行名为ID

  exp_symbol<-inner_join(exp_df,gpl_df, by="ID")

  exp_symbol<-na.omit(exp_symbol) #去掉空值

 # exp_symbol <- separate(exp_symbol,gene_symbol, "symbol",remove = F) # 去掉行中含有多个字符串，如 Dlg1///LOC100047603，只保留Dlg1

  exp_symbol <- separate(exp_symbol,gene_symbol, "gene_symbol",remove = T) # 参数remove=T 表示，去掉表中原有的gene_symbol列

  table(duplicated(exp_symbol$gene_symbol)) # 查看gene_symbol重复数量

  exp_result <-avereps(exp_symbol,ID=exp_symbol$gene_symbol)  %>% data.frame() # 对重复的symbol,取其平均值，同时也有去重功能

  rownames(exp_result) <- exp_result$gene_symbol

  n1 <- ncol(exp_result)-1
  n2 <- ncol(exp_result)

  exp_result <- exp_result[,-c(n1:n2)]
  colnames(exp_result) <- group_infor$title

  exp_df <- as.data.frame(lapply(exp_result,as.numeric)) #将数据框（data.frame）中字符型数据转化为数值型
  rownames(exp_df) <-  rownames(exp_result)


#------------- 保存表达数据文件 ----------------------------------

 write.xlsx(group_infor, paste0("analysis results/",GSE_id,"_groups information.xlsx"), rowNames=T,colNames=T)

 write.xlsx(exp_df, paste0("analysis results/",GSE_id,"_gene_expression_Set.xlsx"), rowNames=T,colNames=T)

 print("--------------------------------------------------------------------")

 print("The gene sets can be found in the folder of <analysis results>")

 print("--------------------------------------------------------------------")


 }



###-----------------  Thne end -----------------------------------##








