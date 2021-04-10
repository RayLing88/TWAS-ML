# <center>**TWAS-ML**: An integrated R toolkit for plant transcriptome-wide association study</center>
<img src="https://github.com/RayLing88/TWAS-ML/tree/main/image/combine.png" title="R logo" style="zoom:50%;" />

### Brief introduction
TWAS-ML is a machine learning-based integrated R toolkit that aims to perfom the plant transcriptome -wide association study (TWAS). This toolkit contains a comprehensive collection of functions required for TWAS data quality control, simulation study, gene expression predictive model building, association study, visualizing for TWAS result. TWAS-ML also takes advantage of machine learning technologies for gene expression levels prediction, with high prediction accuracy, using several ML algorithm such as random forest (RF), support vector machine (SVM), and neural network (NN). 

### <center>Table of Contents</center>
<!-- TOC -->
[TWAS-ML installation](#twas-ml-installation)

[TWAS Data QC](#twas-data-qc)

[TWAS simulation study](#twas-simulation-study)

- Genetic model (additive, dominant, recessive)
- Heritability of gene expression
- Heritablity of trait  
- SNP effect size

[Machine learning model building](#machine-learning-model-building)

- Random Forest
- Support vector machine
- Neural network

[Association](#Association)

[Useful tools](#useful-tools)

[Source codes availability](#source-codes-availability)

[How to access help](#how-to-access-help)

[References](#references)

[Citation](#citation)

<!-- /TOC -->

### TWAS-ML installation
- Step 1: [installation](./tutorial/TWAS_ML_installation.md)
- Step 2: [TWAS-ML installation and quickly start](./tutorial/TWAS-ML_installation_and_quickly_start.md)


### TWAS Data QC
```R
# Loading TWAS data for quality control
cd 01_TWAS_Data_QC/
source("01_TWAS_Data_QC.R")

TWASDataQC <- TWASDataQC(markers = markers,geneExp = geneExp,impute = F,
                        phenotype = phenotype)

```
- this function return high quality TWAS data

### TWAS simualtion study

- Use following command to perform simulation study 

  ```R
  cd 02_Simulation_study/
  Rscript 02_TWAS_Simulation_study.R \
  --GFF gff_file --out simualtion_result --tmp tmpfile\
  --PATH_PLINK yourpath_plink --PATH_GEMMA yourpath_gemma \
  --geno genotype --window_size 50000 --hg2 0.5 --he2 0.5 --eff 0.2 \
  --Gmodel Add --verbose 0 
  ```

- GFF : A simple annotation file containing six columns is like the following table

  | chr  | start | end   | strand | type | GeneID           |
  | ---- | ----- | ----- | ------ | ---- | ---------------- |
  | 1    | 44289 | 49837 | +      | PCG  | *Zm00001d027230* |
  | 1    | 50877 | 55716 | -      | PCG  | *Zm00001d027231* |

- PATH_PLINK: the path of PLINK software; PATH_GEMMA: the path of GEMMA software

- geno: genotype data with PLINK format (.bed, .bim, .fam)

- window_size: cis SNP (x-kb window for a gene)

- hg2, he2, eff, Gmodel were simualted parameters for generating simualtied gene expression data and phenotype data

- The output of this module is like the following table

- | Gmodel   | hg2  | he2​  | ncau | eff  | *p_bslmm* | *p_smr* | *p_rf* | *p_svm* | *p_nn* | *p_lasso* |      |      |      |      |      |      |
  | -------- | ---- | ---- | ---- | ---- | --------- | ------- | ------ | ------- | ------ | --------- | ---- | ---- | ---- | ---- | ---- | ---- |
  | Additive | 0.5  | 0.5  | 1    | 0.3  | 0.38      | 0.25    | 0.03   | 0.01    | 0.02   | 0.04      |      |      |      |      |      |      |
  | Additive | 0.5  | 0.5  | 1    | 0.3  | 0.42      | 0.35    | 0.01   | 0.02    | 0.01   | 0.09      |      |      |      |      |      |      |
  | ...      |      |      |      |      |           |         |        |         |        |           |      |      |      |      |      |      |
  | Additive | 0.5  | 0.5  | 1    | 0.3  | 0.08      | 0.32    | 0.01   | 0.03    | 0.06   | 0.09      |      |      |      |      |      |      |

  

  

### Machine learning model building

- keep high heritable genes for next-step ML model building
- this function is from FUSION software (Gusev et al., 2016)

```R
# 01 heritablity analysis of gene expression
cd 03_Model_building/
Rscript 01_Heritability_analysis.R \
--bfile $INP \
--tmp $TMP \
--out $OUT \
--models top1,blup,bslmm,lasso,enet

# 02 machine learning models (RF,SVM, NN) building
source(02_ML_model_building.R)
runMahineLearn(phenoVec = phenoVec,phenotypes = phenotypes,
              bim = bim,genotype = genotype,
              cpus = cpus,resDic = resDic)

```
### Association

- this module can perfom association study based on individual-level or summary statistics
- this module input file format is just like [PrediXcan software](https://github.com/hakyim/PrediXcan/tree/master/Software) (Gamazon et al., 2015)

```R
#  individual-level association
cd 04_Association/
Rscript 01_Association_individual.R \
--assoc --pheno phenotype_file\
--pred_exp predicted_expression_file\ 
--linear --filter filter_file filter_val --output_dir

# summary-level association
Rscript 02_Association_summary.R \
--sumstats myGWAS.sumstats \
--weights ./WEIGHTS/WEIGHTS.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/PAN503. \
--chr 1 \
--out TWAS.summary.dat

```

- the output of association module is like this table

  | GeneID           | chr  | locus start | locus end | Z score | S.E. | $P_{TWAS}$ |
  | ---------------- | ---- | ----------- | --------- | ------- | ---- | ---------- |
  | *Zm00001d027230* | 1    | 44289       | 49837     | -18.23  | 0.35 | 8.02E-60   |
  | *Zm00001d027231* | 1    | 50877       | 55716     | -9.40   | 0.97 | 1.06E-19   |
  | ...              |      |             |           |         |      |            |
  |                  |      |             |           |         |      |            |

  

### Useful tools

- this module can transform genotype data of Hapmap format to PLINK format

- this module can also visualize TWAS result

   ![](https://github.com/RayLing88/TWAS-ML/tree/main/image/manhatan.png 'manhatan plot')

- locus zoom plot (Zhu et al., 2016)

  ![](https://github.com/RayLing88/TWAS-ML/tree/main/image/locus_zoom.png 'locus zoom')

### Source codes availability

   The source codes of TWAS-ML are freely available at [TWAS-ML](<https://github.com/RayLing88/TWAS-ML>)
### How to access help
* If users encounter any bugs or issues, feel free to leave a message at Github [issues](<https://github.com/cma2015/PEA/issues>). We will try our best to deal with all issues as soon as possible.
* In addition, if any suggestions are available, feel free to contact: __*Lei Ling*__ <linglei@nwafu.edu.cn> or __*Chuang Ma*___ <chuangma2006@gmail.com>

### References
  * Gusev, A., Ko, A., Shi, H. *et al.* Integrative approaches for large-scale transcriptome-wide association studies. *Nat Genet* **48,** 245–252 (2016). https://doi.org/10.1038/ng.3506
  *  Gamazon, E., Wheeler, H., Shah, K. *et al.* A gene-based association method for mapping traits using reference transcriptome data. *Nat Genet* **47,** 1091–1098 (2015). https://doi.org/10.1038/ng.3367
  * Zhu, Z., Zhang, F., Hu, H. *et al.* Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nat Genet* **48,** 481–487 (2016). https://doi.org/10.1038/ng.3538
### Citation
Ling, L., Ma W., Zhang, X., Cheng, Q., Miao, Z., & Ma, C. **TWAS-ML: an integrated R toolkit for plant transcriptome-wide association study.**


