#' @title Visualize result of Transcriptome-wide Association Study (TWAS) 
#' @description  This function can visualize result of TWAS.  

library(optparse)
library(CMplot)

option_list = list(
  make_option("--TWAS_Result", action="store", default=NA, type='character',
              help="GFF annotation file (minus txt) [required]"),
  make_option("--plot_type", action="store", default='m', type='character',
              help="Path to output files [required]"),
  make_option("--output_filetype", action="store", default='pdf', type='character',
              help="Path to output files [required]")
  
)

opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

twas_result <- read.table(opt$TWAS_Result,header = T,as.is = T)

CMplot(Pmap, col=c("#4197d8", "#f8c120", "#413496", "#495226",
                   "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
       bin.size=1e6, bin.range=NULL, pch=19, type="p", band=1, H=1.5, 
       ylim=NULL, cex.axis=1, lwd.axis=1.5, cex.lab=1.5, plot.type=opt$plot_type,
       multracks=FALSE, cex=c(0.5,1,1), r=0.3, outward=FALSE,
       ylab=expression(-log[10](italic(p))), ylab.pos=3, xticks.pos=1,
       mar = c(3,6,3,3), threshold = 0.05/nrow(twas_result), threshold.col="red", threshold.lwd=1, 
       threshold.lty=2, amplify= TRUE, signal.cex = 1.5, signal.pch = 19, 
       signal.col=NULL, signal.line=2, highlight=NULL, highlight.cex=1, 
       highlight.pch=19, highlight.type="p", highlight.col="red", 
       highlight.text=NULL, highlight.text.col="black", highlight.text.cex=1, 
       highlight.text.xadj=NULL, highlight.text.yadj=NULL, 
       highlight.text.font=3, chr.labels=NULL, chr.border=FALSE,
       chr.labels.angle=0, chr.den.col="black", cir.band=1, cir.chr=TRUE, 
       cir.chr.h=1.5, cir.legend=TRUE, cir.legend.cex=0.6, 
       cir.legend.col="black", LOG10=TRUE, box=FALSE, conf.int=TRUE, 
       conf.int.col=NULL, file.output=TRUE, file=opt$output_filetype, 
       dpi=300, height=NULL, width=NULL, memo="", main="", main.cex=1.5, 
       main.font=2, trait.legend.ncol=NULL, verbose=TRUE)