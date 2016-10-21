# Analysis of Alena Shkumatava miRNASeq Data
Anton Enright  
`r format(Sys.Date())`  



# Experiment Setup

All data were pre-processed using *minion* to identify and check adapters, *reaper* to trim adapter sequences followed by *tally* to deduplicate reads while maintaining depth information. Subsequent to this all reads passed through the *mirmod* pipeline against all miRBase (Release 21) precursor sequences for Mouse and Zebrafish. Reads were summed across paired end sequences for the same read pair. Finally reads are loaded into R for final analysis.

This is the sample description file used for the analyses below.

**Name**|**FileName**|**Species**|**SampleName**|**QC**|**Genotype**   
--------|------------|-----------|--------------|------|------------
A543S1|A543S1.R1.fastq.gz|mouse|NREPc1WT0808|qc_all|wt  
A543S2|A543S2.R1.fastq.gz|mouse|NREPc8ho0808|qc_all|scrambled
A543S3|A543S3.R1.fastq.gz|mouse|NREPc1WT1605|qc_all|wt  
A543S4|A543S4.R1.fastq.gz|mouse|NREPc12WT1605|qc_all|wt 
A543S5|A543S5.R1.fastq.gz|mouse|NREPc8ho1605|qc_all|scrambled
A543S6|A543S6.R1.fastq.gz|mouse|NREPc29ho1605|qc_all|scrambled  
A543S7|A543S7.R1.fastq.gz|zebrafish|wt1zf|qc_all|wt
A543S8|A543S8.R1.fastq.gz|zebrafish|wt2zf|qc_all|wt
A543S9|A543S9.R1.fastq.gz|zebrafish|consdelcyrzf1|qc_all|del
A543S10|A543S10.R1.fastq.gz|zebrafish|consdelcyrzf2|qc_all|del  

## Preparation
We first load the R/BioConductor libraries that we need.

```r
library(RColorBrewer)
library(gplots)
library(DESeq2)
library(reshape2)
library(ggplot2)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
spectral <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")(100)
```

# Mouse Analysis

## Count Loading
We can now load all the count data

```r
setwd("/Users/aje/anton_r_notebook/alena_mirna_oct_2016/")
mircounts <- read.table("mouse_counts3.txt",header=TRUE,row.names=1)
mircounts=mircounts[-nrow(mircounts),]
```

As well as the pdata, which contains information on each sample. We will subset for one species at a time, now Mouse

```r
pdata <- read.table("pdata.txt",header=TRUE,row.names=1,sep="\t")
pdata=pdata[pdata$Species=="mouse",]

colnames(mircounts)=rownames(pdata)

conds=as.factor(as.character(pdata$Genotype))

# Filter S1 & S2 ??
#mircounts=mircounts[,c(3:6)]
#pdata=pdata[c(3:6),]
#conds=conds[3:6]
```

## Count Preparation & Normalisation
We are now ready to create a DESeq object from the counts table. 

```r
#Lets Load the Counts First
coldata = as.data.frame(t(t(conds)))
rownames(coldata)=colnames(mircounts)
colnames(coldata)='treatment'
dds <- DESeqDataSetFromMatrix(countData = mircounts, colData = coldata, design = ~ treatment)
```

We are  ready to normalise the data, but first we should look at the number of sequenced reads per sample.

```r
cond_colours = c("#E41A1C","#377EB8")[as.factor(conds)]
names(cond_colours)=conds

group_colours = brewer.pal(length(rownames(pdata)),"Accent")[as.factor(rownames(pdata))]
names(group_colours)=rownames(pdata)

quartz()
barplot(apply(mircounts,2,sum), las=2,col=cond_colours,main="Pre Normalised Counts",cex.names=0.4)
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
```

![](alena_files/figure-html/unnamed-chunk-5-1.png)<!-- -->



We will also estimate the negative binomial dispersion of the data.

```r
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```r
quartz()
plotDispEsts(dds)
```

![](alena_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

## Post Normalisation QC
Now we can normalise and plot the counts again.

```r
normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)


quartz()
barplot(apply(normcounts,2,sum), las=2,col=cond_colours,main="Post-Normalised Counts",cex.names=0.4)
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
```

![](alena_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

We will apply the Variance Stabilising Transformation (VST) it's better than log2 for counts.

```r
vsd <- varianceStabilizingTransformation(dds)
vstcounts <- assay(vsd)
vstcounts <- vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]
```

As an additional QC step we can calculate the sample-to-sample Pearson correlations and plot them in a heatmap.

```r
quartz()
heatmap.2(cor(rawcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (Raw Counts)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours, margins=c(9,7))
```

![](alena_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
quartz()
heatmap.2(cor(vstcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (VST)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours,margins=c(9,7))
```

![](alena_files/figure-html/unnamed-chunk-9-2.png)<!-- -->


We can also perform PCA.

```r
pca <- princomp(vstcounts)


quartz()
par(mfrow=c(2,1))
plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.4)
legend("topright",levels(conds),fill=cond_colours[levels(conds)],cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


PCA of the top 3 Principal Components.

```r
pca2=prcomp(t(vstcounts),center=TRUE)

quartz()
par(mfrow=c(1,3))
plot(pca2$x, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca2$x, as.vector(colnames(mircounts)), pos=3, cex=0.4)
plot(pca2$x[,1],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC1")
text(pca2$x[,1],pca2$x[,3], as.vector(colnames(mircounts)), pos=3, cex=0.4)
plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC2")
text(pca2$x[,2],pca2$x[,3], as.vector(colnames(mircounts)), pos=3, cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

## Initial Biological Analysis of the data
Here are the top10 microRNAs.

```r
top10=apply(mircounts,1,sum)[1:10]
top10[11]=sum(apply(mircounts,1,sum)[11:nrow(mircounts)])
names(top10)[11]="other"
pie(top10,col=brewer.pal(11,"Set3"),main="Top10 microRNAs")
```

![](alena_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

This is the expression of the top10 microRNAs sample to sample.

```r
heatmap.2(vstcounts[names(top10)[1:10],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

This is the expression of the miR-29 families of microRNAs sample to sample.

```r
heatmap.2(vstcounts[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
quartz()
barplot(t(vstcounts[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],]),beside=T,las=2,cex.names=0.5,col=cond_colours,main="miR-29 levels (VST)")
legend("topright",rownames(pdata),fill=cond_colours,cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Statistical Analysis
Run the statistical contrast on the count data

```r
p_threshold=0.05
lfc_threshold=0.75

cds <- nbinomWaldTest(dds)

res=results(cds,contrast=c("treatment","wt","scrambled"))
res <- res[order(res$padj),]
res
```

```
## log2 fold change (MAP): treatment wt vs scrambled 
## Wald test p-value: treatment wt vs scrambled 
## DataFrame with 1513 rows and 6 columns
##                   baseMean log2FoldChange     lfcSE       stat
##                  <numeric>      <numeric> <numeric>  <numeric>
## mmu-mir-196b-5p   138.1979       2.923347 0.3645313   8.019466
## mmu-mir-218-1-5p 6357.5799       1.586746 0.2740549   5.789882
## mmu-mir-150-5p    141.2804      -1.849097 0.3307271  -5.591004
## mmu-mir-139-5p    240.8153      -1.555849 0.2846775  -5.465304
## mmu-mir-10b-5p   4030.2773       2.055345 0.3886272   5.288731
## ...                    ...            ...       ...        ...
## mmu-mir-682-5p   0.6753636    -0.06668301 0.3730918 -0.1787308
## mmu-mir-669k-3p  0.4527021     0.05247907 0.3193682  0.1643216
## mmu-mir-208b-5p  0.2028985    -0.07803852 0.2393486 -0.3260454
## mmu-mir-7017-3p  0.4057970    -0.13881524 0.2429540 -0.5713644
## mmu-mir-7220-5p  0.1134521     0.03479149 0.2393486  0.1453591
##                        pvalue         padj
##                     <numeric>    <numeric>
## mmu-mir-196b-5p  1.062058e-15 6.722826e-13
## mmu-mir-218-1-5p 7.043605e-09 2.229301e-06
## mmu-mir-150-5p   2.257602e-08 4.763540e-06
## mmu-mir-139-5p   4.621129e-08 7.312936e-06
## mmu-mir-10b-5p   1.231676e-07 1.559302e-05
## ...                       ...          ...
## mmu-mir-682-5p      0.8581491           NA
## mmu-mir-669k-3p     0.8694780           NA
## mmu-mir-208b-5p     0.7443900           NA
## mmu-mir-7017-3p     0.5677527           NA
## mmu-mir-7220-5p     0.8844274           NA
```

```r
sig = rownames(res[(abs(res$log2FoldChange) > lfc_threshold) & (res$padj < p_threshold) & !is.na(res$padj),])
```


Volcanoplots of Significant Hits

```r
plot(res$log2FoldChange,-log(res$padj,10),ylab="-log10(Adjusted P)",xlab="Log2 FoldChange",main=paste("Volcano Plot","WT v Scr\nmir-29 in green\nsig. in red"),pch=19,cex=0.4)      
if (length(sig)>=1){
text(res[sig,]$log2FoldChange,-log(res[sig,]$padj,10),labels=rownames(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],]),pos=3,cex=0.6)
points(res[sig,"log2FoldChange"],-log(res[sig,"padj"],10),pch=19,cex=0.4,col="red")
points(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"log2FoldChange"],-log(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"padj"],10),pch=19,cex=0.6,col="green")

}
abline(h=-log10(p_threshold),lty=3)
abline(v=-lfc_threshold,lty=3)
abline(v=lfc_threshold,lty=3)   
```

![](alena_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
plot(apply(vstcounts[,grep("wt",conds)],1,median),apply(vstcounts[,grep("scrambled",conds)],1,median),col="darkblue",cex=0.4,pch=19)
if (length(sig)>1){
points(apply(vstcounts[sig,grep("wt",conds)],1,median),apply(vstcounts[sig,grep("scrambled",conds)],1,median),col="red")
points(apply(vstcounts[rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))],grep("wt",conds)],1,median),apply(vstcounts[rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))],grep("scrambled",conds)],1,median),col="green",cex=1.2)
text(apply(vstcounts[sig,grep("wt",conds)],1,median),apply(vstcounts[sig,grep("scrambled",conds)],1,median),labels=sig,cex=0.4,pos=1)
text(apply(vstcounts[rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))],grep("wt",conds)],1,median),apply(vstcounts[rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))],grep("scrambled",conds)],1,median),labels=rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))],cex=0.4,pos=1,col="darkgreen")
}
abline(a=0,b=1,lty=2,col="red")
```

![](alena_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

Heatmap of significant hits.


```r
heatmap.2(vstcounts[sig,],trace="none",ColSideColors = cond_colours,col=hmcol,margins=c(5,5),cexRow=0.5,cexCol=0.6,labCol=paste(rownames(pdata),pdata$SampleName,sep="\n"),main="Significant Hits Heatmap (VST)")
```

![](alena_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

## Result output to text file
Let's output the final results table with normalised expression values and stats listed

```r
write.table(cbind(as.matrix(counts(dds,normalized=T)[rownames(res),]),as.matrix(res)),"mouse_results.txt",quote=F,sep="\t")
```

## Length Analysis

Now we load in the length data from the mapping analysis separately to analyse.
We will analyse the top 10 expressed, top 10 differential miRs and the miR-29 family, (Excluding miRs with norm counts sum < 30)



```r
length_mouse=read.table("length_tables_mouse.txt",sep="\t",header=FALSE)
length_mouse$genotype=gsub("\\d+.lane.clean.processed.fa.gz","",gsub("mouse_","",length_mouse$V2))

mirlist=unique(c(rownames(counts(dds,normalized=T)[1:10,]),rownames(res[1:15,]),rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))]))

for (i in 1:length(mirlist)){
mir=mirlist[i]

if (median(normcounts[mirlist[i],]) >= 30){
length_table=as.matrix(length_mouse[length_mouse$V1==mir,4:34])/apply(as.matrix(length_mouse[length_mouse$V1==mir,4:34]),1,max)
rownames(length_table)=length_mouse[length_mouse$V1==mir,"genotype"]

length_table=length_table[order(rownames(length_table)),]

colours = c("#E41A1C","#377EB8")[as.factor(rownames(length_table))]
names(colours)=as.factor(rownames(length_table))

heatmap.2(length_table,col=spectral,trace="none",Rowv=F,Colv=F,dendrogram="none",labCol=paste(c(0:30),"nt"),main=paste("Length Table\n",mir),RowSideColors=colours)

barplot(length_table,beside=T,col=colours,names=paste(c(0:30),"nt"),las=2)

matplot(t(length_table),type="b",col=colours,pch=19,cex=0.4,lty=1,lwd=0.4,main=paste("Length Analysis:",mir),xlab=paste(c(0:30),"nt"),las=2,xaxt="n")
axis(1, at = 1:31, labels = paste(c(0:30),"nt"), cex.axis = 0.7,las=2)
legend("topright",levels(as.factor(rownames(length_table))),fill=colours[levels(as.factor(rownames(length_table)))])
}
}
```

![](alena_files/figure-html/unnamed-chunk-20-1.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-2.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-3.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-4.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-5.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-6.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-7.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-8.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-9.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-10.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-11.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-12.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-13.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-14.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-15.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-16.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-17.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-18.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-19.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-20.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-21.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-22.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-23.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-24.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-25.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-26.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-27.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-28.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-29.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-30.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-31.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-32.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-33.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-34.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-35.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-36.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-37.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-38.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-39.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-40.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-41.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-42.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-43.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-44.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-45.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-46.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-47.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-48.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-49.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-50.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-51.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-52.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-53.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-54.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-55.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-56.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-57.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-58.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-59.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-60.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-61.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-62.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-63.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-64.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-65.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-66.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-67.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-68.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-69.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-70.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-71.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-72.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-73.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-74.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-20-75.png)<!-- -->



# Zebrafish Analysis

## Count Loading

We can now load all the count data

```r
setwd("/Users/aje/anton_r_notebook/alena_mirna_oct_2016/")
mircounts <- read.table("zebrafish_counts3.txt",header=TRUE,row.names=1)
mircounts=mircounts[-nrow(mircounts),]
```

As well as the pdata, which contains information on each sample. We will subset for one species at a time, now Zebrafish.

```r
pdata <- read.table("pdata.txt",header=TRUE,row.names=1,sep="\t")
pdata=pdata[pdata$Species=="zebrafish",]

colnames(mircounts)=rownames(pdata)

conds=as.factor(as.character(pdata$Genotype))
```

## Count Preparation & Normalisation

We are now ready to create a DESeq object from the counts table. 

```r
#Lets Load the Counts First
coldata = as.data.frame(t(t(conds)))
rownames(coldata)=colnames(mircounts)
colnames(coldata)='treatment'
dds <- DESeqDataSetFromMatrix(countData = mircounts, colData = coldata, design = ~ treatment)
```

We are  ready to normalise the data, but first we should look at the number of sequenced reads per sample.

```r
cond_colours = c("#E41A1C","#377EB8")[as.factor(conds)]
names(cond_colours)=conds

group_colours = brewer.pal(length(rownames(pdata)),"Accent")[as.factor(rownames(pdata))]
names(group_colours)=rownames(pdata)

quartz()
barplot(apply(mircounts,2,sum), las=2,col=cond_colours,main="Pre Normalised Counts",cex.names=0.4)
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
```

![](alena_files/figure-html/unnamed-chunk-24-1.png)<!-- -->



We will also estimate the negative binomial dispersion of the data.

```r
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```r
quartz()
plotDispEsts(dds)
```

![](alena_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

## Post Normalisation QC
Now we can normalise and plot the counts again.

```r
normcounts <- counts(dds, normalized=TRUE)
rawcounts=counts(dds,normalized=FALSE)
log2counts=log2(normcounts+1)


quartz()
barplot(apply(normcounts,2,sum), las=2,col=cond_colours,main="Post-Normalised Counts",cex.names=0.4)
legend("topright",levels((conds)),cex=0.6,fill=cond_colours[levels(conds)])
```

![](alena_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

We will apply the Variance Stabilising Transformation (VST) it's better than log2 for counts.

```r
vsd <- varianceStabilizingTransformation(dds)
vstcounts <- assay(vsd)
vstcounts <- vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]
```

As an additional QC step we can calculate the sample-to-sample Pearson correlations and plot them in a heatmap.

```r
quartz()
heatmap.2(cor(rawcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (Raw Counts)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours, margins=c(9,7))
```

![](alena_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```r
quartz()
heatmap.2(cor(vstcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (VST)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours,margins=c(9,7))
```

![](alena_files/figure-html/unnamed-chunk-28-2.png)<!-- -->


We can also perform PCA.

```r
pca <- princomp(vstcounts)


quartz()
par(mfrow=c(2,1))
plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.4)
legend("topright",levels(conds),fill=cond_colours[levels(conds)],cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


PCA of the top 3 Principal Components.

```r
pca2=prcomp(t(vstcounts),center=TRUE)

quartz()
par(mfrow=c(1,3))
plot(pca2$x, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca2$x, as.vector(colnames(mircounts)), pos=3, cex=0.4)
plot(pca2$x[,1],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC1")
text(pca2$x[,1],pca2$x[,3], as.vector(colnames(mircounts)), pos=3, cex=0.4)
plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC2")
text(pca2$x[,2],pca2$x[,3], as.vector(colnames(mircounts)), pos=3, cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

```r
par(mfrow=c(1,1))
```

## Initial Biological Analysis of the data
Here are the top10 microRNAs.

```r
top10=apply(mircounts,1,sum)[1:10]
top10[11]=sum(apply(mircounts,1,sum)[11:nrow(mircounts)])
names(top10)[11]="other"
pie(top10,col=brewer.pal(11,"Set3"),main="Top10 microRNAs")
```

![](alena_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

This is the expression of the top10 microRNAs sample to sample.

```r
heatmap.2(vstcounts[names(top10)[1:10],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

This is the expression of the miR-29 families of microRNAs sample to sample.

```r
heatmap.2(vstcounts[rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_files/figure-html/unnamed-chunk-33-1.png)<!-- -->


```r
quartz()
barplot(t(vstcounts[rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))],]),beside=T,las=2,cex.names=0.5,col=cond_colours,main="miR-7 levels (VST)")
legend("topright",rownames(pdata),fill=cond_colours,cex=0.4)
```

![](alena_files/figure-html/unnamed-chunk-34-1.png)<!-- -->

## Statistical Analysis
Run the statistical contrast on the count data

```r
p_threshold=0.05
lfc_threshold=0.75

cds <- nbinomWaldTest(dds)

res=results(cds,contrast=c("treatment","wt","del"))
res <- res[order(res$padj),]
res
```

```
## log2 fold change (MAP): treatment wt vs del 
## Wald test p-value: treatment wt vs del 
## DataFrame with 545 rows and 6 columns
##                      baseMean log2FoldChange      lfcSE        stat
##                     <numeric>      <numeric>  <numeric>   <numeric>
## dre-let-7e-5p       1258.2221      1.1339745  0.1836661    6.174110
## dre-mir-2189-3p     4228.4625      1.1794140  0.2186925    5.393024
## dre-mir-430a-8-5p   2725.7388     -0.9606233  0.2093750   -4.588052
## dre-mir-738-5p       171.4140      0.9621898  0.2263565    4.250772
## dre-let-7c-2-3p      173.1958      0.7843486  0.2163074    3.626083
## ...                       ...            ...        ...         ...
## dre-mir-430b-19-3p 0.09529112   0.0005011098 0.04309577  0.01162782
## dre-mir-124-3-3p   0.48459106  -0.0186131929 0.05350431 -0.34788210
## dre-mir-193a-2-3p  0.09529112   0.0005011098 0.04309577  0.01162782
## dre-mir-124-4-3p   0.46846364  -0.0182240783 0.05293562 -0.34426873
## dre-mir-27e-5p     0.76247246  -0.0591662887 0.06270408 -0.94357963
##                          pvalue         padj
##                       <numeric>    <numeric>
## dre-let-7e-5p      6.653722e-10 3.626279e-07
## dre-mir-2189-3p    6.928177e-08 1.887928e-05
## dre-mir-430a-8-5p  4.474018e-06 8.127799e-04
## dre-mir-738-5p     2.130352e-05 2.902605e-03
## dre-let-7c-2-3p    2.877529e-04 3.136507e-02
## ...                         ...          ...
## dre-mir-430b-19-3p    0.9907226    0.9969942
## dre-mir-124-3-3p      0.7279287    0.9969942
## dre-mir-193a-2-3p     0.9907226    0.9969942
## dre-mir-124-4-3p      0.7306442    0.9969942
## dre-mir-27e-5p        0.3453845    0.9969942
```

```r
sig = rownames(res[(abs(res$log2FoldChange) > lfc_threshold) & (res$padj < p_threshold) & !is.na(res$padj),])
```


Volcanoplots of Significant Hits

```r
plot(res$log2FoldChange,-log(res$padj,10),ylab="-log10(Adjusted P)",xlab="Log2 FoldChange",main=paste("Volcano Plot","WT v Del\nmir-7 in green\nsig. in red"),pch=19,cex=0.4)      
if (length(sig)>=1){
text(res[sig,]$log2FoldChange,-log(res[sig,]$padj,10),labels=rownames(res[rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))],]),pos=3,cex=0.6)
points(res[sig,"log2FoldChange"],-log(res[sig,"padj"],10),pch=19,cex=0.4,col="red")
points(res[rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))],"log2FoldChange"],-log(res[rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))],"padj"],10),pch=19,cex=0.6,col="green")

}
abline(h=-log10(p_threshold),lty=3)
abline(v=-lfc_threshold,lty=3)
abline(v=lfc_threshold,lty=3)   
```

![](alena_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

```r
plot(apply(vstcounts[,grep("wt",conds)],1,median),apply(vstcounts[,grep("del",conds)],1,median),col="darkblue",cex=0.4,pch=19)
if (length(sig)>1){
points(apply(vstcounts[sig,grep("wt",conds)],1,median),apply(vstcounts[sig,grep("del",conds)],1,median),col="red")
points(apply(vstcounts[rownames(vstcounts)[grep("mir-7[a-z]",rownames(vstcounts))],grep("wt",conds)],1,median),apply(vstcounts[rownames(vstcounts)[grep("mir-7[a-z]",rownames(vstcounts))],grep("del",conds)],1,median),col="green",cex=1.2)
text(apply(vstcounts[sig,grep("wt",conds)],1,median),apply(vstcounts[sig,grep("del",conds)],1,median),labels=sig,cex=0.4,pos=1)
text(apply(vstcounts[rownames(vstcounts)[grep("mir-7[a-z]",rownames(vstcounts))],grep("wt",conds)],1,median),apply(vstcounts[rownames(vstcounts)[grep("mir-7[a-z]",rownames(vstcounts))],grep("del",conds)],1,median),labels=rownames(vstcounts)[grep("mir-7[a-z]",rownames(vstcounts))],cex=0.4,pos=1,col="darkgreen")
}
abline(a=0,b=1,lty=2,col="red")
```

![](alena_files/figure-html/unnamed-chunk-36-2.png)<!-- -->

Heatmap of significant hits

```r
heatmap.2(vstcounts[sig,],trace="none",ColSideColors = cond_colours,col=hmcol,margins=c(5,5),cexRow=0.5,cexCol=0.6,labCol=paste(rownames(pdata),pdata$SampleName,sep="\n"),main="Significant Hits Heatmap (VST)")
```

![](alena_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

## Result output to text file
Let's output the final results table with normalised expression values and stats listed

```r
write.table(cbind(as.matrix(counts(dds,normalized=T)[rownames(res),]),as.matrix(res)),"zebrafish_results.txt",quote=F,sep="\t")
```

## Length Analysis

Now we load in the length data from the mapping analysis separately to analyse.
We will analyse the top 10 expressed, top 10 differential miRs and the miR-7 family, (Excluding miRs with norm counts sum <50)


```r
length_zebrafish=read.table("length_tables_zebrafish.txt",sep="\t",header=FALSE)
length_zebrafish$genotype=gsub("\\d+.lane.clean.processed.fa.gz","",gsub("zfish_","",length_zebrafish$V2))

mirlist=unique(c(rownames(counts(dds,normalized=T)[1:10,]),rownames(res[1:15,]),rownames(mircounts)[grep("mir-7[a-z]",rownames(mircounts))]))

for (i in 1:length(mirlist)){
mir=mirlist[i]

if (median(normcounts[mirlist[i],]) >= 50){

length_table=as.matrix(length_zebrafish[length_zebrafish$V1==mir,4:34])/apply(as.matrix(length_zebrafish[length_zebrafish$V1==mir,4:34]),1,max)
rownames(length_table)=length_zebrafish[length_zebrafish$V1==mir,"genotype"]

length_table=length_table[order(rownames(length_table)),]

colours = c("#E41A1C","#377EB8")[as.factor(rownames(length_table))]
names(colours)=as.factor(rownames(length_table))

heatmap.2(length_table,col=spectral,trace="none",Rowv=F,Colv=F,dendrogram="none",labCol=paste(c(0:30),"nt"),main=paste("Length Table\n",mir),RowSideColors=colours)

barplot(length_table,beside=T,col=colours,names=paste(c(0:30),"nt"),las=2)

matplot(t(length_table),type="b",col=colours,pch=19,cex=0.4,lty=1,lwd=0.4,main=paste("Length Analysis:",mir),xlab=paste(c(0:30),"nt"),las=2,xaxt="n")
axis(1, at = 1:31, labels = paste(c(0:30),"nt"), cex.axis = 0.7,las=2)
legend("topright",levels(as.factor(rownames(length_table))),fill=colours[levels(as.factor(rownames(length_table)))])
}
}
```

![](alena_files/figure-html/unnamed-chunk-39-1.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-2.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-3.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-4.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-5.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-6.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-7.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-8.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-9.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-10.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-11.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-12.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-13.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-14.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-15.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-16.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-17.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-18.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-19.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-20.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-21.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-22.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-23.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-24.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-25.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-26.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-27.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-28.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-29.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-30.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-31.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-32.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-33.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-34.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-35.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-36.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-37.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-38.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-39.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-40.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-41.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-42.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-43.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-44.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-45.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-46.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-47.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-48.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-49.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-50.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-51.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-52.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-53.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-54.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-55.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-56.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-57.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-58.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-59.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-60.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-61.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-62.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-63.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-64.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-65.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-66.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-67.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-68.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-69.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-70.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-71.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-72.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-73.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-74.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-75.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-76.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-77.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-78.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-79.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-80.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-81.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-82.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-83.png)<!-- -->![](alena_files/figure-html/unnamed-chunk-39-84.png)<!-- -->

