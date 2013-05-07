#Load libraries
source("http://bioconductor.org/biocLite.R")
library(DEXSeq)
library(Rsamtools)

#Load tables/data frames
setwd("/Users/stevensmith/Documents/School/Maryland/Classes/Comp_Bio_Systems/project/data")
load("exon_gene_mapping.Rda")
load("design_table.Rda")
load("exon_range.Rda")
load("counts.Rda")

#Define functions
plotDispEsts = function( cds, ymin, linecol="#ff000080",
                         xlab = "mean of normalized counts", ylab = "dispersion",
                         log = "xy", cex = 0.45, ... )
{
  px = rowMeans( counts( cds, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]
  
  py = fData(cds)$dispBeforeSharing[sel]
  if(missing(ymin))
    ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
  
  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
       log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
  xg = 10^seq( -.5, 5, length.out=100 )
  fun = function(x) { cds@dispFitCoefs[1] + cds@dispFitCoefs[2] / x }
  lines( xg, fun(xg), col=linecol, lwd=4)
}

plotMA = function(x, ylim,
                  col = ifelse(x$padj>=0.1, "gray32", "red3"),
                  linecol = "#ff000080",
                  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
                  log = "x", cex=0.45, ...)
{
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  
  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}

#Define gene and exon IDs and exon ranges
exon_ID<-exon_gene_mapping$EXONID
gene_ID<-exon_gene_mapping$GENEID

#Modify design table. Change depending on which predictors to use
current<-"class"
design_table<-design_table$class #the data from the design table must be obtained this way, i.e., a factor vector NOT a dataframe

#Define the exon count set/model
vaginal.model<-newExonCountSet(countData=counts,design=design_table,geneIDs=gene_ID,exonIDs=exon_ID,exonIntervals=exon_range)

#Normalize the exon counts
vaginal.model<-estimateSizeFactors(vaginal.model)
sizeFactors(vaginal.model)

#Run QC
#formula=count ~ sample + infection * exon
vaginal.model<-estimateDispersions(vaginal.model)
save(vaginal.model,file="vaginal_model_dispersions_byclass.Rda")
vaginal.model<-fitDispersionFunction(vaginal.model)

#plot Dispersion
#Define file name to be used in output files
fn.plot<-paste("Dispersion_",current,".jpg",sep="")
jpeg(filename=fn.plot,width=900,height=900, pointsize=14)
plotDispEsts(vaginal.model)
dev.off()


#Test for DE usage (change cores accordinly)
vaginal.model<-testForDEU(vaginal.model,nCores=2)
save(vaginal.model,file="vaginal_model_DEU_byclass.Rda")
vaginal.model<-estimatelog2FoldChanges(vaginal.model)
save(vaginal.model,file="vaginal_model_FC_byclass.Rda")
results<-DEUresultTable(vaginal.model)
save(results,file="vaginal_model_results_class.Rda")
results.rmna<-results[!is.na(results$padjust),]

#Write to output
fn<-paste("DEresults_",current,".txt",sep="")
write.table(results,file=fn,sep="\t",quote=F)


table(results.rmna$padjust<0.05)

results.rmna[results.rmna$padjust<0.1,]

plotMA(with(results.rmna,data.frame(baseMean=meanBase,log2FoldChange=`log2fold(Stable/Unstable)`,padj=padjust)),ylim=c(0,1),xlim=c(0,1),cex=0.8)
plotDEXSeq(vaginal.model,"6428")