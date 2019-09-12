library(stringr)
library(cluster)
library(fpc)
library(NMF)
library(ggfortify)
library(Rtsne)
library(ape)
library(pamr)
library(corrplot)

args = commandArgs(trailingOnly=TRUE)
table_file = args[1] #Data Points/Samples on rows and features on columns, must have rownames as first column
labels_file = args[2] #Same number of lines as data points
header_exist = args[3]
sep_str = args[4]
num_clust = as.integer(args[5])
out_folder = args[6]

labels = scan(labels_file,what="character")
labels_fact = as.factor(labels)
labels_fact_col = as.numeric(labels_fact)
lab_uniq = levels(labels_fact)
lab_uniq_col = as.numeric(as.factor(levels(labels_fact)))

table = c()
if(header_exist == "T")
{
table = read.table(table_file,sep = sep_str, row.names = 1, header=T)

} else if(header_exist == "F")
{
table = read.table(table_file,sep = sep_str, row.names = 1, header=F)
colnames(table) = paste("ft",1:ncol(table),sep="_")
} else {print("ERROR")}

#Get clusters#

fit1 = kmeans(table, num_clust,nstart = 25)
clust1 = fit1$cluster
count1 = table(labels,clust1)
fit2 = hclust(dist(table),method = "ward.D2")
clust2 = cutree(fit2, num_clust)
count2 = table(labels,clust2)
############################

#Plot phylo for hclust#
fit_temp = fit2
fit_temp$labels = labels
pdf(paste0(out_folder,"Unrooted_Dendogram.pdf"))
plot(as.phylo(fit_temp),type="unrooted", tip.color =labels_fact_col,label.offset = 3, cex = 0.4,align.tip.label=T)
legend(x = "bottomleft", inset=.02, title="Labels",legend = lab_uniq, fill=lab_uniq_col, horiz=TRUE, cex=0.8)
dev.off()
pdf(paste0(out_folder,"Fan_Dendogram.pdf"))
plot(as.phylo(fit_temp),type="fan", tip.color =labels_fact_col,label.offset = 3, cex = 0.4,align.tip.label=T)
legend(x = "bottomleft", inset=.02, title="Labels",legend = lab_uniq, fill=lab_uniq_col, horiz=TRUE, cex=0.8)
dev.off()
############################
#calculate
pur1 = round(purity(t(count1)),2)
pur2 = round(purity(t(count2)),2)
stat1 = round(cluster.stats(dist(table),clust1)$wb.ratio,2)
stat2 = round(cluster.stats(dist(table),clust2)$wb.ratio,2)
dunn1 = round(cluster.stats(dist(table),clust1)$dunn,2)
dunn2 = round(cluster.stats(dist(table),clust2)$dunn,2)
sil1 = round(cluster.stats(dist(table),clust1)$avg.silwidth,2)
sil2 = round(cluster.stats(dist(table),clust2)$avg.silwidth,2)
############################
#plot

pdf(paste0(out_folder,"Cluster_Histograms_Kmeans.pdf"))
barplot(count1, main=paste("Cluster groupings-Kmeans\n purity = ",pur1,"\nwithin/between = ",stat1,"\ndunn index = ",dunn1,"\nsilhouette width = ",sil1,sep=" "),
        xlab="Cluster", col=lab_uniq_col,cex.main=.7)
legend("bottom", inset=.02, title="Phenotype",
       rownames(count1), fill=lab_uniq_col, horiz=TRUE, cex=0.8,bg="white")
dev.off()
pdf(paste0(out_folder,"Cluster_Histograms_Hclust.pdf"))
barplot(count2, main=paste("Cluster groupings-Hclust\n purity = ",pur2,"\nwithin/between = ",stat2,"\ndunn index = ",dunn2,"\nsilhouette width = ",sil2,sep=" "),
        xlab="Cluster", col=lab_uniq_col,cex.main=.7)
legend("bottom", inset=.02, title="Phenotype",
       rownames(count2), fill=lab_uniq_col, horiz=TRUE, cex=0.8,bg="white")
        #legend.text = T,args.legend=list(x="right",bty = "o",horiz=F,xpd=T,inset=c(-.5,0)),cex.main=.7)
dev.off()
pdf(paste0(out_folder,"PCA.pdf"))
table_temp = table
colnames(table_temp) = gsub(pattern = "-",replacement = "__",x = colnames(table_temp))
autoplot(prcomp(table_temp),data = data.frame(table_temp,labels),colour="labels")
dev.off()
pdf(paste0(out_folder,"TSNE.pdf"))
tsne_out = Rtsne(as.matrix(table),perplexity = 10)
plot(tsne_out$Y,col=data.frame(table,labels)$labels,asp=1,pch=20)
legend(x = "bottomleft", inset=.02, title="Labels",legend = unique(data.frame(table,labels)$labels), fill=unique(data.frame(table,labels)$labels), horiz=TRUE, cex=0.8)
dev.off()
############################
table_pam = t(table)
data_for_pam_table = list(x=table_pam[complete.cases(table_pam),],y =labels_fact,genenames=rownames(table_pam[complete.cases(table_pam),]),geneid = rownames(table_pam[complete.cases(table_pam),]))
model = pamr.train(data_for_pam_table)
model.cv = pamr.cv(fit = model, data = data_for_pam_table)
pdf(paste0(out_folder,"PAM_CV.pdf"))
pamr.plotcv(model.cv)
dev.off()
thresh1 = round(quantile(model.cv$threshold,0),2)
thresh2 = round(quantile(model.cv$threshold,0.25),2)
thresh3 = round(quantile(model.cv$threshold,0.75),2)
thresh4 = round(quantile(model.cv$threshold,0.95),2)
thresh50 = round(quantile(model.cv$threshold,0.50),2)
pdf(paste0(out_folder,"PAM_Confusion_Matrices.pdf"))
par(mfrow=c(2,2))
corrplot(pamr.confusion(model.cv, threshold=thresh1,F)/rowSums(pamr.confusion(model.cv, threshold=thresh1,F)),method="pie",tl.col = "black",mar = c(2,2,2,2),title=paste0("Threshold = ",thresh1))
corrplot(pamr.confusion(model.cv, threshold=thresh2,F)/rowSums(pamr.confusion(model.cv, threshold=thresh2,F)),method="pie",tl.col = "black",mar = c(2,2,2,2),title=paste0("Threshold = ",thresh2))
corrplot(pamr.confusion(model.cv, threshold=thresh3,F)/rowSums(pamr.confusion(model.cv, threshold=thresh3,F)),method="pie",tl.col = "black",mar = c(2,2,2,2),title=paste0("Threshold = ",thresh3))
corrplot(pamr.confusion(model.cv, threshold=thresh4,F)/rowSums(pamr.confusion(model.cv, threshold=thresh4,F)),method="pie",tl.col = "black",mar = c(2,2,2,2),title=paste0("Threshold = ",thresh4))
dev.off()
pdf(paste0(out_folder,"PAM_CV_probs.pdf"))
usr <- par( "usr" )
pamr.plotcvprob(model.cv, data_for_pam_table, threshold=thresh1)
text( usr[ 1 ], usr[ 4 ], paste0("Threshold = ",thresh1), adj = c( 0, 1 ) )
pamr.plotcvprob(model.cv, data_for_pam_table, threshold=thresh2)
text( usr[ 1 ], usr[ 4 ], paste0("Threshold = ",thresh2), adj = c( 0, 1 ) )
pamr.plotcvprob(model.cv, data_for_pam_table, threshold=thresh3)
text( usr[ 1 ], usr[ 4 ], paste0("Threshold = ",thresh3), adj = c( 0, 1 ) )
pamr.plotcvprob(model.cv, data_for_pam_table, threshold=thresh4)
text( usr[ 1 ], usr[ 4 ], paste0("Threshold = ",thresh4), adj = c( 0, 1 ) )
dev.off()
pdf(paste0(out_folder,"PAM_Centroid_Threshold_",thresh50,".pdf"))
par(mfrow=c(1,1))
pamr.plotcen(model, data_for_pam_table, threshold=thresh50)
dev.off()
pdf(paste0(out_folder,"PAM_Example_Features_Threshold_",thresh3,".pdf"))
pamr.geneplot(model, data_for_pam_table, threshold=thresh3)
dev.off()


