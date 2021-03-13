#' @title Simulation study of Transcriptome-wide Association Study (TWAS) 
#' @description  This function can perform simulation study of TWAS.  


library(plink2R)
library(optparse)
library(randomForest)
library(e1071)
library(AMORE)
library(impute)
library(glmnet)

option_list = list(
  make_option("--GFF", action="store", default=NA, type='character',
              help="GFF annotation file (minus txt) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--PATH_PLINK", action="store", default=NA, type='character',
              help="Path to PLINK software [required]"),
  make_option("--PATH_GEMMA", action="store", default=NA, type='character',
              help="Path to GEMMA software [required]"),
  make_option("--geno", action="store", default="plink", type='character',
              help="Path to genotype data of PLINK format  [%default]"),
  make_option("--window_size", action="store", default="50000", type='double',
              help="cis SNPs of each gene [%default]"),
  make_option("--hg2", action="store", default="0.5", type='double',
              help="Heritability of phenotype [%default]"),
  make_option("--he2", action="store", default="0.5", type='double',
              help="Heritability of gene expression [%default]"),
  make_option("--eff", action="store", default="0.2", type='double',
              help="effect size (SNP to gene expression) [%default]"),
  make_option("--Gmodel", action="store", default="Add", type='character',
              help="Genetic model used for generating simulated TWAS data (Add,Dom,Rec) [%default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]")
)

opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

####################################
geno <- read_plink(opt$geno)
genotype <- geno$bed
bim <- geno$bim
fam <- geno$fam
rm(geno)
bim <- bim[match(unique(bim$V2),bim[,2]),]
genotype <- genotype[,bim[,2]]
if (length(is.na(genotype)) > 0 ){
  genotype <- impute.knn(genotype,k = 10,rowmax = 0.5,
                         colmax = 0.8,maxp = 1500, rng.seed = 362436069)
  genotype <- genotype$data
}


## read gff file (chr,left,right,strand,type,geneID)
gff <- read.table(opt$GFF,header = F,as.is = T)

CisSNPCount <- function(chr ,gff,size){
  gff <- gff[which(gff[,1] == chr),]
  rownames(gff) <- gff$V6
  result <- NULL
  for (gene in gff$V6){
    
      P0 <- as.numeric(gff[gene,2]) - as.numeric(size)
      P1 <- as.numeric(gff[gene,3]) + as.numeric(size)
      if (P0 <0){P0 = 1}
    
    cmd <- paste0(opt$PATH_PLINK,' --bfile ',opt$geno, ' --make-bed --out ',opt$tmp,' --chr ',bim[1,1], ' --from-bp ',P0,' --to-bp ',P1,' --extract ',opt$geno,'.bim')
    system(cmd,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    #cat(P0,P1,'\n')
    if (file.exists(paste0(opt$tmp,'.bim'))){
      tmp_bim <- read.table(paste0(opt$tmp,'.bim'),header = F,as.is = T)
      cis_snp <- tmp_bim[,2]
      result <- rbind(result,c(gene,length(cis_snp)))
      system(paste0('rm -r ',opt$tmp,'*'))
    }else{
      next
    }
    
    
  }
  result <- data.frame(result)
  result[,2] <- as.numeric(result[,2])
  return(result)
 
}

tmp <- CisSNPCount(chr = bim[1,1],gff = gff,size = opt$window_size)
cat('Task of sampling genes associated with trait has been done!!!!!!-----------\n')
gene_name <- tmp[ sample(which(tmp$X2 > 300 & tmp$X2 < 600 ),50),1]
subgff <- gff[match(gene_name,gff$V6),]
rownames(subgff) <- subgff$V6

### Task of preparing genotype data has been done!!!!!!!


SimulationStudy <- function(subgff = subgff,genotype = genotype,
                              bim = bim,resDic,ngwas,K = 5,
                              hg2=0.5,he2 = 0.22,size,ncau = 0.01,
                              big_eff = 0.3,Gmodel = 'Add'){
  
  
  
  map <- bim[,c(2,1,4)]
  rownames(map) <- map[,1]
  
  if (ngwas == nrow(genotype)){
    subgenotype <- genotype
    result <- NULL
    Acc_result <- NULL
    P_result <- NULL
  }else{
    subgenotype <- genotype[sample(1:nrow(genotype),ngwas),]
    #subgenotype <- matrix(,nrow = ngwas,ncol = ncol(genotype))
    #colnames(subgenotype) <- colnames(genotype)
    result <- NULL
    P_result <- NULL
    for (j in 1:ncol(subgenotype)){
      
      index2 <- sample(1:ngwas,round(0.35*ngwas))
      subgenotype[index2,j] = 2
      subgenotype[setdiff(1:ngwas,index2),j] = 0
    }
  }
  
  
  #subgenotype <- genotype[sample_index,]
  
  
  for (i in 1:nrow(subgff)){
    lasso_count = 0;gwas_count = 0;smr_count = 0
    svm_count  = 0;rf_count = 0;nn_count = 0;fusion_count = 0
    
    curgff <- subgff[i,]
    gene <- curgff[1,6]
    P0 <- as.numeric(curgff[gene,2]) - as.numeric(size)
    P1 <- as.numeric(curgff[gene,3]) + as.numeric(size)
    if (P0 <0 ){P0 = 1}
    ### simulate gene expression
    #t1 <- Sys.time()
    
    
    cmd <- paste0(opt$PATH_PLINK,' --bfile ',opt$geno, ' --make-bed --out ',opt$tmp,' --chr ',bim[1,1], ' --from-bp ',P0,' --to-bp ',P1,' --extract ',opt$geno,'.bim')
    system(cmd,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    tmp_bim <- read.table(paste0(opt$tmp,'.bim'),header = F,as.is = T)
    cis_snp <- tmp_bim[,2]
    
    # if (length(cis_snp) < 100){
    #   next
    # }
    
    causal_SNP <- sample(cis_snp,1)
    X  = subgenotype[,causal_SNP]
    #beta = rnorm(length(causal_SNP),0,sd = sqrt(he2))
    X = as.matrix(X)
    
    if (Gmodel == 'Dom'){
      ef_vec <- vector(mode = 'numeric',length = nrow(X))
      idx <- which(X[,1] == names(which.max(table(X[,1]))) )
      ef_vec[idx] = big_eff
      X[X ==0] = -2
      E <- X * ef_vec + rnorm(nrow(X),0,sqrt(1-he2))
      
    }else if(Gmodel == 'Rec'){
      ef_vec <- vector(mode = 'numeric',length = nrow(X))
      idx <- which(X[,1] == names(which.min(table(X[,1]))) )
      ef_vec[idx] = big_eff
      E <- X * ef_vec + rnorm(nrow(X),0,sqrt(1-he2))
    }else{
      Gmodel = 'Add'
      E <- X*big_eff + rnorm(nrow(X),0,sqrt(1-he2))
    }
    #min_eff <- runif((length(causal_SNP) -1),0,0.1)
   
    #E = X %*% beta + rnorm(nrow(X),0,sd = sqrt(1-he2))
    y = E * rnorm(1,0,sd = sqrt(hg2)) + rnorm(nrow(subgenotype),0,sqrt(1-hg2))
    
    ####
    fam <- read.table(paste0(opt$tmp,'.fam'),header = F,as.is = T)
    fam[,6] <- y
    ### phenotype determine
    write.table(fam,file = paste0(opt$tmp,'.fam'),
                quote = F,col.names = F,row.names = F)
    ### run gwas
    gwas_cmd <- paste0(opt$PATH_GEMMA,' -bfile ',opt$tmp,' -lm 1 -c pca.eigenvec -o ',opt$tmp)
    system(gwas_cmd,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    gwas <- read.table(paste0('output/',opt$tmp,'.assoc.txt'),header = T,as.is = T)
    gwas_snp = gwas[which.min(gwas$p_wald),2]
    z_gwas = qnorm(gwas[match(gwas_snp,gwas$rs),'p_wald'] / 2 )
    
    ## run eqtl
    fam <- read.table(paste0(opt$tmp,'.fam'),header = F,as.is = T)
    fam[,6] <- E
    ### phenotype determine
    write.table(fam,file = paste0(opt$tmp,'.fam'),
                quote = F,col.names = F,row.names = F)
    ###
    cmd <- paste0(opt$PATH_GEMMA,' -bfile ',opt$tmp,' -lm 1 -c pca.eigenvec -o ',opt$tmp)
    system(cmd,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    eqtl <- read.table(paste0('output/',opt$tmp,'.assoc.txt'),header = T,as.is = T)
    
    ### run SMR
    #eqtl <- oracle(E,subgenotype[,cis_snp],confounder = confounder,K=5)
    eqtl_snp = eqtl[which.min(eqtl$p_wald),2]
    z_eqtl = qnorm(eqtl[match(eqtl_snp,eqtl$rs),'p_wald'] / 2)
    
    T_smr <- (z_eqtl^2 * z_gwas^2)/(z_eqtl^2 + z_gwas^2)
    #z_smr = (z_eqtl * z_gwas)/sqrt(z_eqtl^2 + z_gwas^2)
    p_smr = 1 - pchisq(T_smr,df=1)
    
   
    
    if(is.na(p_smr)){
      p_smr = NA
    }
    if (p_smr < 0.05){
      smr_count <- smr_count + 1
    }
    
    
    F_X = subgenotype[,cis_snp]
    F_X = apply(F_X, 2, as.numeric)
    index <- sample(1:nrow(genotype),round(0.7*nrow(genotype)))
    X_train = F_X[index,]
    X_test = F_X[-index,]
    E_train = E[index]
    E_test = E[-index]
    
    
    #### fusion 
    # fam <- read.table(paste0(opt$tmp,'.fam'),header = F,as.is = T)
    # fam[,6] <- E
    # write.table(fam,file = paste0(opt$tmp,'.fam'),
    #             quote = F,col.names = F,row.names = F)
    bslmm_cmd <- paste0(opt$PATH_GEMMA,' -miss 1 -maf 0 -r2 1 -bfile ',opt$tmp,' -bslmm 2 -o ',opt$tmp)
    system(bslmm_cmd,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
    weights <- read.table(paste0('output/',opt$tmp,'.param.txt'),header = T,as.is = T)
    weights <- weights$alpha + weights$beta*weights$gamma
    E_bslmm_pre <- X_test %*% weights + 1
    bslmm_r <- cor(E_bslmm_pre,E_test)
    
    lm_fusion <- lm(y[-index] ~ E_bslmm_pre)
    p_twas <- summary(lm_fusion)$coefficients[2,4]
   
    if (p_twas < 0.05){
      fusion_count  = fusion_count + 1
    }
    system(paste0('rm -r ',opt$tmp,'*'))
    system(paste0('rm -r output/',opt$tmp,'*'))
    #### svm
    
    model1 = svm(X_train,E_train,type = "nu-regression",gamma = 1/ncol(X_train))
    svm_y = predict(model1,F_X)
    svm_r = cor(svm_y,E)
    lm_svm <- lm(y ~ svm_y)
    p_svm <- summary(lm_svm)$coefficients[2,4]
    
    if (p_svm < 0.05){
      svm_count = svm_count + 1
    }
    ### net 
    net <- newff(n.neurons=c(ncol(X_train),50,50,1),learning.rate.global=1e-4,momentum.global=0.001,error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin",
                 method="ADAPTgdwm")
    nn_model <- train(net,X_train,as.numeric(E_train),error.criterium="LMS", report=TRUE,n.shows = 5,
                      n.threads = 10,show.step = 100)
    NN_y <-sim(nn_model$net,F_X)
    nn_r <- cor(NN_y,E)
    lm_nn <- lm(y ~ NN_y)
    p_nn <- summary(lm_nn)$coefficients[2,4]
     
    if (p_nn < 0.05){
      nn_count = nn_count + 1
    }
    
    rf_model = randomForest(x = X_train,y = E_train,ntree = 300)
    rf_y = predict(rf_model,F_X)
    rf_r <- cor(E,rf_y)
    lm_rf <- lm(y ~ rf_y)
    p_rf <- summary(lm_rf)$coefficients[2,4]
    
    if (p_rf < 0.05){
      rf_count = rf_count + 1
    }
    
    ### enet
    
    enet = cv.glmnet( x= X_train,y = E_train,alpha=0.5,nfold=5,intercept=T,standardize=F )
    enet_y = predict(enet,X_test)
    enet_y = as.vector(enet_y)
    enet_r <- cor(E[-index],enet_y)
    tmp1 = summary(lm(y[-index] ~ enet_y))$coefficients
    if (length(tmp1) != 8){
      cat("SKip enet prediction of gene:",gene,'\n')
      p_enet = 0.5
    }else{
      p_enet =  as.vector(tmp1)[8]
    }
    
    
    
    if (p_enet < 0.05){
      lasso_count = lasso_count + 1
    }
    tmp =matrix(c(Gmodel,hg2,he2,ncau,big_eff,ngwas,
                  p_twas,p_smr,p_svm,p_rf,p_nn,p_enet,
                  bslmm_r,svm_r,rf_r,nn_r,enet_r),nrow = 1,ncol = 17)
    #tmp2 <- matrix(c(p_smr,pmr_result$causal_pvalue,
                     #p_svm,p_rf,p_nn,p_enet,p_twas),nrow = 1,ncol = 7)
    #P_result <- rbind(P_result,tmp2)
    CurAcc <- c(bslmm_r,svm_r,nn_r,rf_r,enet_r)
    names(CurAcc) <- c('bslmm','svm','nn','rf','enet')
    Acc_result <- rbind(Acc_result,c(hg2,he2,ncau,ngwas,CurAcc))
    
    result <- rbind(result,tmp)
    #cat(tmp,'\n')
    cat(gene,"finishing",match(gene,subgff[,6]),"round" ,"\n")
    cat(names(CurAcc),'\n')
    cat(CurAcc,'\n')
    cat(c(p_twas,p_svm,p_nn,p_rf,p_enet),'\n')
    #write.table(result,file = paste0(resDic,gene),quote = F,
    # row.names = F,col.names = F)
  }
  return(result = result)
  
}

### perform simulation study


### pow1 based on causality genetic model

simulation_result <- SimulationStudy(subgff = subgff,genotype,bim,ngwas = nrow(genotype),K = 5,
                            Gmodel = opt$Gmodel,hg2 = opt$hg2,he2 = opt$he2,
                            big_eff = opt$eff,
                            size = opt$window_size)
  
colnames(simulation_result) <- c('Gmodel','hg2','he2','ncau','eff','ngwas',
                             'p_twas','p_smr','p_svm','p_rf','p_nn','p_enet',
                             'bslmm_r','svm_r','rf_r','nn_r','enet_r')
write.table(simulation_result,
            file = paste0(opt$out,'_',opt$Gmodel,'_',opt$hg2,'_',
                          opt$he2,'_',opt$eff,'_.txt'),
            quote = F,col.names = T,row.names = F)









