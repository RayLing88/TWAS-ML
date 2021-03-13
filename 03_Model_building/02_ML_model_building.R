#' @title Machine lerning models building to predict gene expression levels
#' @description  This function can construct ML models (RF,SVM,NN) to predict gene expression levels  

#' 
#' @author Lei Ling
#' @keywords machine learning

ReadGeneExp <- function(x){
  #' @title To read gene expression levels
  #' @description: Read a gene expression file
  #' @details Input an gene expression filename
  #' @param x A character
  #' @return A DataFrame which stores gene expression levels
  #' @export ReadGeneExp
  require(data.table)
  geneExp <- fread(x,header = T)
  geneExp <- data.frame(geneExp)
  rownames(geneExp) <- geneExp$ID
  return(geneExp)
}

ReadBinaryGenotype <- function(plink_prefix){
  #' @title Read binary genotype data in plink format
  #' @description: Read binary genotype data in plink format
  #' @details Input prefix of plink genotype data
  #' @param x A character
  #' @return A list which stores genotype data
  #' @export ReadBinaryGenotype
  require(plink2R)
  geno <- read_plink(plink_prefix)
  return(geno)
}
# ara_fusion <- fread('ara_FUSION_expr.txt',header = T)
# ara_fusion <- data.frame(ara_fusion)
# rownames(ara_fusion) <- ara_fusion$ID

CVgroup <- function(k,datasize,seed){
  #' description: Do k-fold cross validation
  #' @param k (numeric, integer) k-fold (such as k = 5)
  #' @param datasize (numeric,integer) the total numbers of sample
  #' @param seed (numeric,integer) the random seed
  #' @return A list storing k-fold samples split result
  #' @export CVgroup
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]
  temp <- sample(n,datasize)
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}

# CalR2 <- function(y,y_pre){
#    R2 <- sum((y_pre -mean(y))^2)/ sum( (y-mean(y))^2)
#    return(R2)
# }

MachineLearn <- function(i, phenotypes,genotype,bim,size = 1e5,k = 5,resDic){
  #' @title Build machine learning models for cis snp and gene expression levels
  #' @description: Build machine learning models for cis snp and gene expression levels
  #' @param i (numeric,integer) which represents the ith gene in gene expression matrix
  #' @param phenotypes (DataFrame) each rows represents a gene  and each column means a sample
  #' @param genotype (numeric,matrix). Genotypes are coded as {0,1,2}. Row is sample and column is SNP
  #' @param bim (DataFrame) SNP annotation file which can be from ReadBinaryGenotype
  #' @param size (numeric,integer) bp . cis window size such as 1e5 bp in ath.
  #' @param k (numeric, integer) k-fold (such as k = 5)
  #' @export MachineLearn
  ## phenotypes ia a matrix the fourth col is gene annotation
  curY <- t(phenotypes[i,5:ncol(phenotypes)])
  gene_start = as.vector(as.matrix(phenotypes[i,4]))
  cis_snp = bim[which(bim$V4 > gene_start -size & bim$V4 < gene_start + size),2]
  curX <- as.matrix(genotype[,cis_snp])
  tmp <- CVgroup(k = k,datasize = nrow(curX),seed = 123)
  
  if (!is.data.frame(curX)){
    curX = data.frame(curX)
  }
  result = matrix(0,nrow = k,ncol = 3)
  colnames(result) <- c('RF','SVM','NN')
  
  for(j in 1:k){
    X_train <- curX[-tmp[[j]],]
    X_test <- curX[tmp[[j]],]
    y_train <- curY[-tmp[[j]]]
    y_test <- curY[tmp[[j]]]
    require("randomForest")
    rf_model <- randomForest(X_train,y_train)
    y_pre <- predict(rf_model,X_test)
    
    result[j,1] <- cor(y_pre,y_test)^2
    #### SVM modeling
    require('e1071')
    model1 = svm(X_train,y_train,type = "nu-regression",gamma = 1/ncol(X_train))
    svm_y = predict(model1,X_test)
    result[j,2] <- cor(svm_y,y_test)^2
    #### NN Modeling
    require('AMORE')
    net <- newff(n.neurons=c(ncol(X_train),20,1),
                 learning.rate.global=1e-4,
                 momentum.global=0.001,error.criterium="LMS",
                 Stao=NA, hidden.layer="tansig", output.layer="purelin",
                 method="ADAPTgdwm")
    nn_model <- AMORE::train(net,X_train,y_train,error.criterium="LMS", report=TRUE,n.shows = 5,
                             n.threads = 10,show.step = 100)
    NN_y <-sim(nn_model$net,X_test)
    result[j,3] <- cor(NN_y,y_test)^2
  }
  result_data <- c(phenotypes[i,c(1,3,4)],colMeans(result) )
  if (is.list(result_data)){
    result_data <- do.call(c,result_data)
    result_data <- t(result_data)
  }
  write.table(result_data,file = paste0(resDic,rownames(phenotypes)[i], ".txt"),
              quote = F,col.names = F,row.names = F)
  
}

runMachineLearn <- function(phenoVec, phenotypes, genotype,bim,size = 1e5,k = 5, cpus,
                            resDic){
  #' @title Parallel building machine learning models
  #' @description: Parallel build machine learning models
  #' @param phenoVec (numeric,vector) such as 1:100
  #' @param phenotypes (DataFrame) each rows represents a gene  and each column means a sample
  #' @param genotype (numeric,matrix). Genotypes are coded as {0,1,2}. Row is sample and column is SNP
  #' @param bim (DataFrame) SNP annotation file which can be from ReadBinaryGenotype
  #' @param size (numeric,integer) bp . cis window size such as 1e5 bp in ath.
  #' @param k (numeric, integer) k-fold (such as k = 5)
  #' @param cpus (integer) the number of cpus used in parallel computing
  #' @param resDic (character) the output file path
  #' @export runMachineLearn
  require(snowfall)
  if(cpus > 1){
    
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("MachineLearn")
    sfExport("CVgroup")
    sfLibrary("randomForest", character.only = TRUE )
    sfLibrary('AMORE',character.only = TRUE)
    sfLibrary("e1071",character.only = TRUE)
    
    res <- sfSapply(x = phenoVec, MachineLearn, phenotypes = phenotypes,
                    genotype = genotype,bim = bim,
                    size = 1e5,k =5,resDic = resDic)
    sfStop()
  }else{
    for(i in 1:length(phenoVec)){
      res <- sfLapply(x = phenoVec[i],MachineLearn, phenotypes = phenotypes,
                      genotype = genotype,
                      bim  = bim,size = 1e5,k =5,resDic = resDic)
    }
  }
}


tt <- runMachineLearn(phenoVec = phenoVec, phenotypes = sub,bim = bim,genotype = genotype,
                         cpus = 1,resDic = resDic)


cat('Finishing ML task!!! congratulations!!!!!!!!!')