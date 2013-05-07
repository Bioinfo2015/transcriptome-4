#Load model/results from using class
load("vaginal_model_results_class.Rda")
results.class<-results
load("vaginal_model_FC_byclass.Rda")
vaginal.model.class<-vaginal.model

sum(is.na(results.class$padjust))/nrow(results.class)

table(results.class$padjust<0.05)
results.class<-results.class[!is.na(results.class$padjust),]
top.exons.class<-results.class[results.class$padjust<0.05,]
write.table(top.exons.class,file="top_exons_from_class.txt",sep="\t",quote=F,row.names=F)
top.exons.class.nonzero<-top.exons.class$geneID[top.exons.class$padjust>0]

generate_geneview_plots(top.exons.class.nonzero,"class",vaginal.model.class)

#generate general stability
#Load model/results from using stability
load("vaginal_model_results_global_stability.Rda")
results.globalStability<-results
load("vaginal_model_FC.Rda")
vaginal.model.globalStability<-vaginal.model

sum(is.na(results.globalStability$padjust))/nrow(results.globalStability)

table(results.globalStability$padjust<0.05)
results.globalStability<-results.globalStability[!is.na(results.globalStability$padjust),]
top.exons.globalStability<-results.globalStability[results.globalStability$padjust<0.05,]
write.table(top.exons.globalStability,file="top_exons_from_globalStability.txt",sep="\t",quote=F,row.names=F)
top.exons.globalStability.nonzero<-top.exons.globalStability$geneID[top.exons.globalStability$padjust>0]

generate_geneview_plots(top.exons.globalStability.nonzero,"globalStability",vaginal.model.globalStability)

generate_geneview_plots<-function(top.exons.nonzero,title,vaginal.model){
gene.IDs<-unique(top.exons.nonzero)
for (gene in gene.IDs){
  fn.plot<-paste("GeneView_",gene,"_",title,".jpg",sep="")
  jpeg(filename=fn.plot,width=900,height=900, pointsize=14)
  plotDEXSeq(vaginal.model,gene,legend=TRUE,displayTranscripts=TRUE)
  dev.off()
}
}
head(results.globalStability)
plotDEXSeq(vaginal.model.globalStability,"51750",legend=TRUE,displayTranscripts=TRUE)
