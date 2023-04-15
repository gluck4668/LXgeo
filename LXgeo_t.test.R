

LXgeo <- function(GSE_id,Platforms){

#--------------------------------------------
 options (stringsAsFactors = F)#不让字符串转化为因子
 Sys.setenv("VROOM_CONNETION_SIZE"=19999999999)#为内存不足，运行代码
 options('download.file.method.GEOquery'='libcurl')

 # options ( "repos "="https://mirrors.ustc.edu.cn/CRAN/")
 # options(Bioc_mirror="https://mirrors.ustc.edu.cn/bioc/")


 all_packages <- data.frame(installed.packages())

 pack <- data.frame(c("BiocManager","openxlsx","psych","dplyr","tidyverse","stringr",
                      "pacman","umap","nortest","car","Cairo"))

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
 
 options('download.file.method.GEOquery'='libcurl')
 
 if(!file.exists(gse_file)){
    gse <- getGEO (GSE_id, destdir = "analysis results/", 
                   GSEMatrix =TRUE, getGPL=FALSE)

    save(gse,file=gse_file)
    }

 load(gse_file)

#-------- GE0数据整理 ---------------
 
 if (length(gse) > 1) idx <- grep(Platforms, attr(gse, "names")) else idx <- 1
 
 exp_matrix <- exprs(gse[[idx]]) #获取表达矩阵

 if(nrow(exp_matrix)==0)
      {stop_text <- paste("There is no gene expression data in",GSE_id,"series matrix",". Please check it." )
      stop(stop_text)}

 group_infor <-pData(gse[[1]]) #获取分组等实验基本数据

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
  options('download.file.method.GEOquery'='libcurl')
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

#-------------------------------------------------------------  
  # log2 transform
  
  qx <- as.numeric(quantile(exp_df, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  
  LogC <- (qx[5] > 100) ||  (qx[6]-qx[1] > 50 && qx[2] > 0)
  
  if (LogC) { exp_df[which(exp_df <= 0)] <- NaN
              ex_log2 <- log2(exp_df) }
  
  # box-and-whisker plot
  par(mar=c(7,4,2,1))
  title <- paste (GSE_id, "/", Platforms, sep ="")
  
  dev.off() 
  
  Cairo::CairoPNG( 
    filename = paste0("analysis results/boxplot_",GSE_id,".png"), # 文件名称
    width = 9,           # 宽
    height = 7,          # 高
    units = "in",        # 单位
    dpi = 300)           # 分辨率
  
  boxplot(ex_log2, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
  
  dev.off() 
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  title <- paste (GSE_id, "/", Platforms, " value distribution", sep ="")
  
  Cairo::CairoPNG( 
    filename = paste0("analysis results/plotDensities_",GSE_id,".png"), # 文件名称
    width = 8,           # 宽
    height = 9,          # 高
    units = "in",        # 单位
    dpi = 300)           # 分辨率
  
  plotDensities(ex_log2, main=title, legend=F)
  
  dev.off()
  
  # mean-variance trend
  ex_log2 <- na.omit(ex_log2) # eliminate rows with NAs
  
  Cairo::CairoPNG( 
    filename = paste0("analysis results/plotSA_",GSE_id,".png"), # 文件名称
    width = 8,           # 宽
    height = 9,          # 高
    units = "in",        # 单位
    dpi = 300)           # 分辨率
  
  plotSA(lmFit(ex_log2), main=paste0("Mean variance trend, ",GSE_id))
  
  dev.off()
  
 #---t.test------------------------- 
  
  exp_df_test <- exp_df
  
  exp_df_test$id <- rownames(exp_df)
  
  nn <- ncol(exp_df_test) %>% as.numeric()
  
  exp_df_test <- exp_df_test[,c(nn,1:nn-1)]
  
  rownames(exp_df_test) <- c(1:nrow(exp_df_test))
  
  df <- limma::avereps(exp_df_test[,-1],exp_df_test[,1]) %>% data.frame() # 对重复的项，取其平均值，同时也有去重功能
  
  group <- colnames(df)
  group <- gsub("\\d+$","",group)  #去掉字符串后面的数字
  table <- data.frame(table(group))
  group <- distinct(data.frame(group))
  group <- dplyr::left_join(group,table,"group")
  
  df$mean1 <- NA #添加列名
  df$sd1 <- NA
  
  df$mean2 <- NA
  df$sd2 <- NA
  
  df$FoldChange <- NA
  df$pvalue <- NA
  
  colnames(df)[ncol(df)-5] <- paste0(group[1,1],"_mean") #改列名
  colnames(df)[ncol(df)-4] <- paste0(group[1,1],"_sd")
  
  colnames(df)[ncol(df)-3] <- paste0(group[2,1],"_mean")
  colnames(df)[ncol(df)-2] <- paste0(group[2,1],"_sd")
  
  df[ncol(df)-5] <- rowMeans(df[,c(1:group[1,2])]) # 计算row的平均值
  
  n1 <- ncol(df)-group[2,2]-6+1 # group2起始列
  n2 <- group[1,2]+group[2,2] %>% as.numeric() # group2起末列
  df[ncol(df)-3] <- rowMeans(df[,c(n1:n2)])
  
  df[ncol(df)-1] <- df[ncol(df)-3]/df[ncol(df)-5] # 计算FoldChange
  
  
  #计算sd和p值
  n1_sd <- ncol(df)-4 # 第一组sd所在的列
  n2_sd <- ncol(df)-2 # 第二组sd所在的列
  
  equal_variance <- df
  equal_variance$equal_variance <- NA
  equal_variance$normality_test <- NA
  
  for(i in 1:nrow(df)){
    df[i,n1_sd] <- sd(df[i,c(1:group[1,2])]) # 计算第一组sd
    df[i,n2_sd] <- sd(df[i,c(n1:n2)]) # 计算第二组sd
    
    test_group <- c(rep(group[1,1],group[1,2]),rep(group[2,1],group[2,2]))
    test_data <- df[i,c(1:n2)] %>% as.numeric()
    
    # 正态检验(SPSS 规定: 当样本含量3 ≤ n ≤ 5000时, 结果以Shapiro-Wilk为准, 当样本含量n > 5000结果以Kolmogorov-Smirnov为准)
    n <- table(test_data) %>% sum() %>% as.numeric()
    
    if(n<=5000){
      normality_test <- shapiro.test(test_data)  # shapiro.test(x)；x为数据，长度为3-5000。
      norm.test.p <- normality_test[["p.value"]]} else
      {normality_test <- lillie.test(test_data) #nortest包的lillie.test函数:(Kolmogorov-Smirnov)正态性检验
      norm.test.p <- normality_test[["p.value"]]}
    
    # 方差齐性分析
    #对于正态分布的样本，Bartlette检验极其灵敏，但是对于非正态分布的样本，检验非常不准确；
    #equal_var <- bartlett.test(test_data,test_group)
    #equal.var.p <- equal_var[["p.value"]]
    
    #Levene检验是一种更为稳健的检验方法，既可用于正态分布的样本，也可用于非正态分布的样本，
    #同时对比较的各组样本量可以相等或不等；是sPSs的黑默t认方差齐性检验方法。
    equal_var <- leveneTest(test_data~factor(test_group))
    equal.var.p <- equal_var[1,3]
    
    #前者是对原始数据的方差进行检验的，leveneTest是对方差模型的残差进行组间齐性检验.
    #一般认为是要求残差的方差齐，所以一般的统计软件都做的是leveneTest
    
    if(equal.var.p>0.05){
      t_test <- t.test(df[i,c(1:group[1,2])],df[i,c(n1:n2)],
                       paired = FALSE,
                       var.equal = T, # 方差齐 var.equal = T
                       conf.level = 0.95) } else
                       {t_test <- t.test(df[i,c(1:group[1,2])],df[i,c(n1:n2)],
                                         paired = FALSE,
                                         var.equal = F, # 方差不齐 var.equal = F
                                         conf.level = 0.95)}
    
    df[i,ncol(df)] <- t_test[["p.value"]]
    
    equal_variance[i,ncol(equal_variance)] <- norm.test.p
    equal_variance[i,ncol(equal_variance)-1] <- equal.var.p
    equal_variance[i,ncol(equal_variance)-2] <- t_test[["p.value"]]
    equal_variance[i,ncol(equal_variance)-4] <- sd(df[i,c(n1:n2)])
    equal_variance[i,ncol(equal_variance)-6] <- sd(df[i,c(1:group[1,2])])
    
  }
  
  
  df$p.adjust <- p.adjust(df[,ncol(df)],method = "BH")
  
  
#------------- 保存表达数据文件 ----------------------------------
 write.xlsx(group_infor, paste0("analysis results/",GSE_id,"_groups information.xlsx"), rowNames=T,colNames=T)
 write.xlsx(df, paste0("analysis results/",GSE_id,"_gene_expression_Set.xlsx"), rowNames=T,colNames=T)

 print("--------------------------------------------------------------------")

 print("The gene sets can be found in the folder of <analysis results>")

 print("--------------------------------------------------------------------")


 }



###-----------------  Thne end -----------------------------------##








