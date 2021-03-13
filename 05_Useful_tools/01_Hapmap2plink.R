
library(optparse)
option_list = list(
  make_option("--PATH_PLINK", action="store", default=NA, type='character',
              help="Path to PLINK software [required]"),
  make_option("--PATH_VCFTOOLS", action="store", default=NA, type='character',
              help="Path to vcftools software [required]"),
  make_option("--PATH_TASSEL", action="store", default=NA, type='character',
              help="Path to tassel software [required]"),
  make_option("--Hapmap", action="store", default="", type='character',
              help="Hapmap file  [%default]"),
  make_option("--VCF", action="store", default="", type='character',
              help="VCF file  [%default]"),
  make_option("--tmp", action="store", default="", type='character',
              help="cis SNPs of each gene [%default]"),
  make_option("--out", action="store", default="", type='character',
              help="output file [%default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]")
)

opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

cmd1 <- paste0(opt$PATH_TASSEL,' -Xms5g -Xmx10g -SortGenotypeFilePlugin -inputFile ',opt$Hapmap,' -outputFile ',opt$tmp,' -fileType Hapmap')
system(cmd1,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

cmd2 <- paste0(opt$PATH_TASSEL,' -Xmx10g -fork1 -h ',opt$tmp,' -export ',opt$VCF,' -exportType VCF -runfork1')
system(cmd2,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

cmd3 <- paste0(opt$PATH_VCFTOOLS,' --vcf ',opt$VCF,' --plink --out ',opt$out)
system(cmd3,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

cmd4 <- paste0(opt$PATH_PLINK,' --file ',opt$out,' --make-bed --out ',opt$out)
system(cmd4,ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
