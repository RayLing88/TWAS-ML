#' @title Data Quality control of Transcriptome-wide Association Study (TWAS) 
#' @description  This function can examine and summary the quality of TWAS data. And can be used for imputation for genotype.  
#' @param markers  (numeric, matrix)row is sample well column is SNP information (feature).Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param phenotype  (numeric)the phenotype value of each individual.
#' @param impute  (logical)if TRUE, imputation, default F.
#' @param filter (logical)if TRUE, filter the markers data with the setting of MAF and missing value.
#' @param NABound (numeric, [0,1])max missing value percentage.
#' @param MAFBound (numerix, [0,1])min MAF percentage in each marker local sit.
#' @param imputeMethod  (character)the method of imputation, "mean", "median" or "KNN", default "mean".
#' @param round  (numeric)rounds the values in its first argument to the specified number of decimal places, default 4.
#' @param k,maxp,rowmax,colmax,rng.seed (numeric)KNN method parameters.
#' @seealso
#' \pkg{impute}
#' @return
#' A list of the data quality information.
#' 
#' @author Lei Ling, Qian Cheng
#' @keywords Quality Control
#' @export 

TWASDataQC <- function(markers, geneExp,phenotype, impute = F, filter = F, NABound = 0.8, MAFBound = 0.005, imputeMethod = "mean", round = 4 ,k = 10, maxp = 1500, rowmax = 0.5, colmax = 0.8, rng.seed = 362436069){
  
  if (! is.matrix(markers)){
    markerFormat <- warning(paste0("The formation of markers is wrong! Please check!","\n"))
  }else{
    markerFormat <- cat(paste0("The formation of markers is right!","\n"))
  }
  if (! is.numeric(markers)){
    markersClass <- warning(paste0("The data in matrix is not numeric! Please transform!","\n"))
  }else{
    markersClass <- cat(paste0("The data in matrix is right!","\n"))
  }
  if(! is.numeric(geneExp)){
    geneExpFormat <- warning(paste0("The geneExp data is not numeric!Please transform!","\n"))
  }
  if(! is.numeric(phenotype)){
    phenotypeFormat <- warning(paste0("The phenotype data is not numeric!Please transform!","\n"))
  }
  if (is.matrix(markers) & is.numeric(markers) & is.numeric(phenotype)) {
    MAF <- round(MAFSummary(markers),2)
    if (length(which(is.na(markers)) == TRUE) == 0){
      NACount <- 0 
      NAPercentTotal <- 0
      NACountRow <- 0
      NACountCol <- 0
      NACountRowPercent <- 0
      NACountColPercent <- 0
      resImpute <- markers
      NASummary <- cat(paste0("The marker matrix have no missing value!","\n"))
      NAIdx <- NA
      tableRes <- round(table(markers))
    }else{
      NACount <- round(length(which(is.na(markers)) == TRUE))
      NAPercentTotal <- round(NACount/(ncol(markers)*nrow(markers)),round)
      
      NACountRow <- round(apply(markers,1,function(x){length(which(is.na(x)) == TRUE)}))
      NACountCol <- round(apply(markers,2,function(x){length(which(is.na(x)) == TRUE)}))
      
      NACountRowPercent <- round(NACountRow/ncol(markers),2)
      NACountColPercent <- round(NACountCol/nrow(markers),2)
      
      NAIdx <- which(is.na(markers) == TRUE)
      
      tableRes <- round(table(markers))
      
      if(filter){
        NAFliIdx <- which(NACountColPercent >= NABound)
        MAFFliIdx <- which(MAF <= MAFBound)
        fliterIdx <- unique(c(NAFliIdx,MAFFliIdx))
        markers <- markers[,-fliterIdx]
      }
      
      if(impute){
        resImpute <- GSImputation(x = markers, imputeMethod = imputeMethod , k = k,rowmax = rowmax,colmax = colmax, maxp = maxp, rng.seed = rng.seed)
      }else{
        resImpute <- markers
      }
    }
    
    phenotypeNACount <- round(length(which(is.na(phenotype) == TRUE)))
    phenotypeNAIdx <- which(is.na(phenotype) == TRUE)
    
    NASum <- c(NACount,NAPercentTotal)
    names(NASum) <- c("missValueCount","missValuePercent") 
    
    markersTable <- c(tableRes,NASum)
    
    ## summarize
    #     summarize <- paste0(markerFormat,markersClass,"The data have ",NACount," missing value!","\n","Missing value percent:",NAPercentTotal,"\n",
    #                             "The phenotype data have ",phenotypeNACount," missing value!","\n",
    #                             "The data have ",length(tableRes)," element.","\n","The detail in markerReport!",
    #                             "\n")
    markerReport <- list(Global = markersTable, Imputation = resImpute,MAF = MAF, NACountRow = NACountRow,NACountRowPercent = NACountRowPercent,NACountCol = NACountCol,
                         NACountColPercent = NACountColPercent, NAIdx = NAIdx,penotypeNAIdx = phenotypeNAIdx)
    cat(paste0(markerFormat,markersClass,"The data have ",NACount," missing value!","\n","Missing value percent:",NAPercentTotal,"\n",
               "The phenotype data have ",phenotypeNACount," missing value!","\n",
               "The data have ",length(tableRes)," element.","\n","The detail in markerReport!",
               "\n"))
    markerReport
  }else{
    stop(paste0("The format of data is wrong! There may be the following above cases:","\n","The format of maker is not a matrix;","\n",
                "The class of markers or phenotype is not numeric."))
  }
}
Imputation <- function(x, imputeMethod = "mean", k = 10,rowmax = 0.5,colmax = 0.8,maxp = 1500, rng.seed = 362436069){
  require("impute")
  if(is.numeric(x)){
    if(imputeMethod == "KNN"){
      requireNamespace("impute", quietly = TRUE)
      x <- impute.knn(data = t(x), k = k, rowmax = rowmax, colmax = colmax, maxp = maxp, rng.seed = rng.seed)
      x <- round(t(x$data))
    }else{
      imputeMethod <- imputeMethod
    }
    if(imputeMethod == "mean"){
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }else if(imputeMethod=="median"){
      x[which(is.na(x))] <- median(x,na.rm=TRUE)
    }
  }else{
    if(imputeMethod =="mean"){
      stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
    }else if(imputeMethod=="median"){
      tt <- table(x)
      x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
    }else{
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }
  }
  return(x)
}
MAFSummary <- function(x){
  MAF <- apply(x,2,function(x){
    level <- table(x)
    if (length(level) == 2){
      index <- order(level,decreasing = T)
      MAF <- level[index[2]]/sum(level)
    } else if(length(level) == 3){
      index <- order(level[-2],decreasing = T)
      MAF <- (level[-2][index[2]]*2 + level[2])/(sum(level)*2)
    }
  })
  MAF
}

ReadGeneExp <- function(x){
  #' @title To read gene expression matrix
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
  #' @return A list which stores genotype data (bed,bim,fam)
  #' @export ReadBinaryGenotype
  require(plink2R)
  geno <- read_plink(plink_prefix)
  return(geno)
}

FilterLowExp <- function(ExpMat = ExpMat,cutoff = 1,nsize = 5){
  #' @title Filter genes having low expression levels
  #' @param ExpMat geneExpMat from ReadGeneExp function
  #' @param cutoff cutoff of TPM value
  #' @param nsize the number of samples whose TPM value > cutoff 
  #' @return A high level of gene expression matrix 
  ### row represents a  gene and column represent a sample
  ### and input from fread function, thus the first col is gene ID
  ExpMat <- data.frame(ExpMat)
  # transcript_id <- ExpMat$Transcript_ID
  # ExpMat <- ExpMat[,-1]
  ExpMat <- t(ExpMat) ### row is sample and col is gene
  tmp <- apply(ExpMat, 2, function(x){length(which(x > cutoff))})
  keep_index <- which(tmp > nsize)
  cat('Remove the number of genes is :',ncol(ExpMat) - length(keep_index),'\n')
  cat('Keep the number of genes: ',length(keep_index),'\n')
  HighExpMat <- ExpMat[,keep_index] 
  colnames(HighExpMat) = transcript_id[keep_index]
  return(list(HighExpMat = HighExpMat,keep_index = keep_index))
}

QuantileNormalization <- function(ExpMat){
  #' @title Filter genes having low expression levels
  #' @param ExpMat geneExpMat from ReadGeneExp function
  #' @return A gene expression matrix processed by quantile normalization 
  
  quantile_mat <- preprocessCore::normalize.quantiles(ExpMat,copy = F)
  col_names = colnames(ExpMat)
  #New_mat = cbind(New_mat[,1:4],quantile_mat)
  colnames(quantile_mat) = col_names
  return(quantile_mat)
}

