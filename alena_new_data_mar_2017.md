# Analysis of Alena Shkumatava miRNASeq Data
Anton Enright  
`r format(Sys.Date())`  



# Experiment Setup

All data were pre-processed using *minion* to identify and check adapters, *reaper* to trim adapter sequences followed by *tally* to deduplicate reads while maintaining depth information. Subsequent to this all reads passed through the *mirmod* pipeline against all miRBase (Release 21) precursor sequences for Mouse and Zebrafish. Reads were summed across paired end sequences for the same read pair. Finally reads are loaded into R for final analysis.

This is the sample description file used for the analyses below.

**Name**|**File**|**Barcodes**|**3p_ad**
----|----|--------|-----
wt1|A638S1.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA
wt2|A638S2.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA
wt3|A638S3.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA
mt1|A638S4.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA
mt2|A638S5.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA
mt3|A638S6.R1.fastq.gz|no_barcode|AGATCGGAAGAGCACA

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
mircounts <- read.table("mouse_counts_mar_2017.txt",header=TRUE,row.names=1)
mircounts=mircounts[-nrow(mircounts),]
```

As well as the pdata, which contains information on each sample.

```r
pdata <- read.table("pdata_mar_2017.txt",header=TRUE,row.names=1)

#pdata=pdata[c(1,2,5,6),]
#mircounts=mircounts[,c(1,2,5,6)]

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

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-5-1.png)<!-- -->



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

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

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

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

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

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
quartz()
heatmap.2(cor(vstcounts),trace="none",col=hmcol,main="Sample to Sample Correlation (VST)",cexRow=0.5,cexCol=0.5,RowSideColors=cond_colours,margins=c(9,7))
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-9-2.png)<!-- -->


We can also perform PCA.

```r
pca <- princomp(vstcounts)


quartz()
plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.4)
legend("topright",levels(conds),fill=cond_colours[levels(conds)],cex=0.4)
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


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

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

## Initial Biological Analysis of the data
Here are the top10 microRNAs.

```r
top10=apply(mircounts,1,sum)[1:10]
top10[11]=sum(apply(mircounts,1,sum)[11:nrow(mircounts)])
names(top10)[11]="other"
pie(top10,col=brewer.pal(11,"Set3"),main="Top10 microRNAs")
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

This is the expression of the top10 microRNAs sample to sample.

```r
heatmap.2(vstcounts[names(top10)[1:10],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

This is the expression of the miR-29 families of microRNAs sample to sample.

```r
heatmap.2(vstcounts[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],],col=hmcol,trace="none",cexCol=0.4,cexRow=0.6,ColSideColors=cond_colours)
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
quartz()
barplot(t(vstcounts[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],]),beside=T,las=2,cex.names=0.5,col=cond_colours,main="miR-29 levels (VST)")
legend("topright",rownames(pdata),fill=cond_colours,cex=0.4)
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Statistical Analysis
Run the statistical contrast on the count data

```r
p_threshold=0.05
lfc_threshold=0.75

cds <- nbinomWaldTest(dds)

res=results(cds,contrast=c("treatment","wt","mut"))
res <- res[order(res$padj),]
res
```

```
## log2 fold change (MAP): treatment wt vs mut 
## Wald test p-value: treatment wt vs mut 
## DataFrame with 1471 rows and 6 columns
##                      baseMean log2FoldChange     lfcSE        stat
##                     <numeric>      <numeric> <numeric>   <numeric>
## mmu-mir-708-5p      370.21877      -3.880212 0.5208669   -7.449528
## mmu-mir-219-2-3p   2508.84257      -3.586149 0.6161118   -5.820615
## mmu-mir-204-5p      595.61606       3.605966 0.6224714    5.792982
## mmu-mir-219-2-5p     89.54399      -2.535594 0.4502091   -5.632036
## mmu-mir-10b-5p   632527.44731      -2.852560 0.5180408   -5.506439
## ...                       ...            ...       ...         ...
## mmu-mir-875-3p     0.12715984    -0.42370466  1.021558 -0.41476320
## mmu-mir-882-5p     0.22104968    -0.42370466  1.021558 -0.41476320
## mmu-mir-489-5p     0.08632002     0.05517826  1.021749  0.05400371
## mmu-mir-804-3p     0.21019972     0.26690316  1.024176  0.26060278
## mmu-mir-142-5p     0.12715984    -0.42370466  1.021558 -0.41476320
##                        pvalue         padj
##                     <numeric>    <numeric>
## mmu-mir-708-5p   9.367487e-14 8.814805e-11
## mmu-mir-219-2-3p 5.863146e-09 2.168921e-06
## mmu-mir-204-5p   6.914734e-09 2.168921e-06
## mmu-mir-219-2-5p 1.780944e-08 4.189672e-06
## mmu-mir-10b-5p   3.661657e-08 6.891239e-06
## ...                       ...          ...
## mmu-mir-875-3p      0.6783153           NA
## mmu-mir-882-5p      0.6783153           NA
## mmu-mir-489-5p      0.9569322           NA
## mmu-mir-804-3p      0.7943988           NA
## mmu-mir-142-5p      0.6783153           NA
```

```r
sig = rownames(res[(abs(res$log2FoldChange) > lfc_threshold) & (res$padj < p_threshold) & !is.na(res$padj),])
```


Volcanoplots of Significant Hits

```r
plot(res$log2FoldChange,-log(res$padj,10),ylab="-log10(Adjusted P)",xlab="Log2 FoldChange",main=paste("Volcano Plot","WT v Scr\nmir-29 in green\nsig. in red"),pch=19,cex=0.4)      
points(res[sig,"log2FoldChange"],-log(res[sig,"padj"],10),pch=19,cex=0.4,col="red")
text(res[sig[1:10],"log2FoldChange"],-log(res[sig[1:10],"padj"],10),pch=19,cex=0.4,pos=2,labels = rownames(res[sig[1:10],]))
points(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"log2FoldChange"],-log(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"padj"],10),pch=19,cex=0.6,col="green")
text(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"log2FoldChange"],-log(res[rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))],"padj"],10),pch=19,cex=0.4,pos=2,labels =rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))])
abline(h=-log10(p_threshold),lty=3)
abline(v=-lfc_threshold,lty=3)
abline(v=lfc_threshold,lty=3)   
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

## Scatter Plot

```r
wt_median = apply(vstcounts[,pdata$Genotype == "wt"],1,median)
mt_median = apply(vstcounts[,pdata$Genotype == "mut"],1,median)
plot(wt_median,mt_median,cex=0.4,pch=19,col="darkblue")
points(wt_median[grep("mir-29[a-z]",rownames(vstcounts))],mt_median[grep("mir-29[a-z]",rownames(vstcounts))],cex=0.4,pch=19,col="green")
points(wt_median[sig],mt_median[sig],cex=1,col="red")
text(wt_median[grep("mir-29[a-z]",rownames(vstcounts))],mt_median[grep("mir-29[a-z]",rownames(vstcounts))],cex=0.4,pos=3,labels=rownames(vstcounts)[grep("mir-29[a-z]",rownames(vstcounts))])
abline(a=0,b=1,lty=2,col="red")
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

Heatmap of significant hits.


```r
heatmap.2(vstcounts[sig,],trace="none",ColSideColors = cond_colours,col=hmcol,margins=c(5,5),cexRow=0.5,cexCol=0.6,labCol=paste(rownames(pdata),pdata$SampleName,sep="\n"),main="Significant Hits Heatmap (VST)")
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

## Result output to text file
Let's output the final results table with normalised expression values and stats listed

```r
write.table(cbind(as.matrix(counts(dds,normalized=T)[rownames(res),]),as.matrix(res)),"mouse_results.txt",quote=F,sep="\t")
```

## Length Analysis

Now we load in the length data from the mapping analysis separately to analyse.
We will analyse the top 5 expressed, top 5 differential miRs and the miR-29 family, (Excluding miRs with norm counts sum < 50)



```r
length_mouse=read.table("length_tables_mouse_mar_2017.txt",sep="\t",header=FALSE)
length_mouse$genotype=gsub("\\d+.lane.clean.uniquified.fa.gz","",gsub("mouse_","",length_mouse$V2))

mirlist=unique(c(rownames(counts(dds,normalized=T)[1:5,]),rownames(res[1:5,]),rownames(mircounts)[grep("mir-29[a-z]",rownames(mircounts))]))

for (i in 1:length(mirlist)){
mir=mirlist[i]

if (median(normcounts[mirlist[i],]) >= 50){
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


for (i in 16:18){
  for (j in 21:23){
    pvalue=NA
    pvalue=t.test(length_table[1:3,i:j],length_table[4:6,i:j])$p.value
    if (pvalue <= 0.05){
      sig="*"
    } else {
      sig="ns"
    }
    
    print(paste(mir," trimming-test ",i,"-",j," P-value:",pvalue," >>",sig,sep=""))
  }
}

for (i in 21:23){
  for (j in 26:28){
    pvalue=NA;
    pvalue=t.test(length_table[1:3,i:j],length_table[4:6,i:j])$p.value
    if (pvalue <= 0.05){
      sig="*"
    } else {
      sig="ns"
    }
    print(paste(mir," tailing-test ",i,"-",j," P-value:",pvalue," >>",sig,sep=""))
  }
}


}
}
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-1.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-2.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-3.png)<!-- -->

```
## [1] "mmu-mir-9-2-5p trimming-test 16-21 P-value:0.979147074753133 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 16-22 P-value:0.997584574003662 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 16-23 P-value:0.872215266188831 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 17-21 P-value:0.978627989650754 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 17-22 P-value:0.997543317648514 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 17-23 P-value:0.870460917312873 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 18-21 P-value:0.976687176024219 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 18-22 P-value:0.997273946801105 >>ns"
## [1] "mmu-mir-9-2-5p trimming-test 18-23 P-value:0.868021928599759 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 21-26 P-value:0.969498108549098 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 21-27 P-value:0.973160588838365 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 21-28 P-value:0.975019372480275 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 22-26 P-value:0.965984306687916 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 22-27 P-value:0.969632325471596 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 22-28 P-value:0.971431465613599 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 23-26 P-value:0.957795772136665 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 23-27 P-value:0.961556476420575 >>ns"
## [1] "mmu-mir-9-2-5p tailing-test 23-28 P-value:0.963322451065909 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-4.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-5.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-6.png)<!-- -->

```
## [1] "mmu-mir-148a-3p trimming-test 16-21 P-value:0.355061167786251 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 16-22 P-value:0.657956590637664 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 16-23 P-value:0.980677619822352 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 17-21 P-value:0.291448067194535 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 17-22 P-value:0.634011657381611 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 17-23 P-value:0.980094789720396 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 18-21 P-value:0.197331324938934 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 18-22 P-value:0.609321564452911 >>ns"
## [1] "mmu-mir-148a-3p trimming-test 18-23 P-value:0.98387447784671 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 21-26 P-value:0.9996680157295 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 21-27 P-value:0.999693250985347 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 21-28 P-value:0.999712230961797 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 22-26 P-value:0.999693299174668 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 22-27 P-value:0.99971287174773 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 22-28 P-value:0.999728711316913 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 23-26 P-value:0.999610327278387 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 23-27 P-value:0.999626098586824 >>ns"
## [1] "mmu-mir-148a-3p tailing-test 23-28 P-value:0.99964078082836 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-7.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-8.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-9.png)<!-- -->

```
## [1] "mmu-let-7f-2-5p trimming-test 16-21 P-value:0.943810792007168 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 16-22 P-value:0.91881834296344 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 16-23 P-value:0.933561746596804 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 17-21 P-value:0.940904058642806 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 17-22 P-value:0.913389013473437 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 17-23 P-value:0.931012007977119 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 18-21 P-value:0.935918373652784 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 18-22 P-value:0.903851511590817 >>ns"
## [1] "mmu-let-7f-2-5p trimming-test 18-23 P-value:0.927226996467204 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 21-26 P-value:0.995530815031397 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 21-27 P-value:0.995824276636423 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 21-28 P-value:0.996006052860702 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 22-26 P-value:0.996102031040233 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 22-27 P-value:0.996302361169858 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 22-28 P-value:0.996423992743325 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 23-26 P-value:0.997080446133547 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 23-27 P-value:0.997172797066433 >>ns"
## [1] "mmu-let-7f-2-5p tailing-test 23-28 P-value:0.997226453803235 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-10.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-11.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-12.png)<!-- -->

```
## [1] "mmu-let-7i-5p trimming-test 16-21 P-value:0.848743000061334 >>ns"
## [1] "mmu-let-7i-5p trimming-test 16-22 P-value:0.844411336808158 >>ns"
## [1] "mmu-let-7i-5p trimming-test 16-23 P-value:0.910115295354445 >>ns"
## [1] "mmu-let-7i-5p trimming-test 17-21 P-value:0.840133559493801 >>ns"
## [1] "mmu-let-7i-5p trimming-test 17-22 P-value:0.837499463693066 >>ns"
## [1] "mmu-let-7i-5p trimming-test 17-23 P-value:0.909125506674347 >>ns"
## [1] "mmu-let-7i-5p trimming-test 18-21 P-value:0.819208390830819 >>ns"
## [1] "mmu-let-7i-5p trimming-test 18-22 P-value:0.823666226938496 >>ns"
## [1] "mmu-let-7i-5p trimming-test 18-23 P-value:0.907473362525591 >>ns"
## [1] "mmu-let-7i-5p tailing-test 21-26 P-value:0.999631198394098 >>ns"
## [1] "mmu-let-7i-5p tailing-test 21-27 P-value:0.999656782931391 >>ns"
## [1] "mmu-let-7i-5p tailing-test 21-28 P-value:0.999681092171831 >>ns"
## [1] "mmu-let-7i-5p tailing-test 22-26 P-value:0.999800445996549 >>ns"
## [1] "mmu-let-7i-5p tailing-test 22-27 P-value:0.999811793155087 >>ns"
## [1] "mmu-let-7i-5p tailing-test 22-28 P-value:0.999828451552709 >>ns"
## [1] "mmu-let-7i-5p tailing-test 23-26 P-value:0.999870545885176 >>ns"
## [1] "mmu-let-7i-5p tailing-test 23-27 P-value:0.99987291778907 >>ns"
## [1] "mmu-let-7i-5p tailing-test 23-28 P-value:0.999861400997893 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-13.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-14.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-15.png)<!-- -->

```
## [1] "mmu-mir-21a-5p trimming-test 16-21 P-value:0.602235033654314 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 16-22 P-value:0.875857746593749 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 16-23 P-value:0.917650529732843 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 17-21 P-value:0.581016781598493 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 17-22 P-value:0.870092010270946 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 17-23 P-value:0.916435421364573 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 18-21 P-value:0.541905497352372 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 18-22 P-value:0.860933872359785 >>ns"
## [1] "mmu-mir-21a-5p trimming-test 18-23 P-value:0.914599969177603 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 21-26 P-value:0.978822137736616 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 21-27 P-value:0.981888601203621 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 21-28 P-value:0.983468238099258 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 22-26 P-value:0.977756369469599 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 22-27 P-value:0.980685341727589 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 22-28 P-value:0.98217182969097 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 23-26 P-value:0.97478917417996 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 23-27 P-value:0.977619016817308 >>ns"
## [1] "mmu-mir-21a-5p tailing-test 23-28 P-value:0.979022948296825 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-16.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-17.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-18.png)<!-- -->

```
## [1] "mmu-mir-708-5p trimming-test 16-21 P-value:0.112790349868752 >>ns"
## [1] "mmu-mir-708-5p trimming-test 16-22 P-value:0.154882207059156 >>ns"
## [1] "mmu-mir-708-5p trimming-test 16-23 P-value:0.192627887133128 >>ns"
## [1] "mmu-mir-708-5p trimming-test 17-21 P-value:0.10221552457954 >>ns"
## [1] "mmu-mir-708-5p trimming-test 17-22 P-value:0.142530288939413 >>ns"
## [1] "mmu-mir-708-5p trimming-test 17-23 P-value:0.179762189828129 >>ns"
## [1] "mmu-mir-708-5p trimming-test 18-21 P-value:0.123503063060658 >>ns"
## [1] "mmu-mir-708-5p trimming-test 18-22 P-value:0.157020025305967 >>ns"
## [1] "mmu-mir-708-5p trimming-test 18-23 P-value:0.18895193472291 >>ns"
## [1] "mmu-mir-708-5p tailing-test 21-26 P-value:0.73621757694301 >>ns"
## [1] "mmu-mir-708-5p tailing-test 21-27 P-value:0.775933469842892 >>ns"
## [1] "mmu-mir-708-5p tailing-test 21-28 P-value:0.791908454715575 >>ns"
## [1] "mmu-mir-708-5p tailing-test 22-26 P-value:0.743671739470898 >>ns"
## [1] "mmu-mir-708-5p tailing-test 22-27 P-value:0.779706751718029 >>ns"
## [1] "mmu-mir-708-5p tailing-test 22-28 P-value:0.792859599163142 >>ns"
## [1] "mmu-mir-708-5p tailing-test 23-26 P-value:0.735798875936854 >>ns"
## [1] "mmu-mir-708-5p tailing-test 23-27 P-value:0.768939569502606 >>ns"
## [1] "mmu-mir-708-5p tailing-test 23-28 P-value:0.779285957154257 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-19.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-20.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-21.png)<!-- -->

```
## [1] "mmu-mir-219-2-3p trimming-test 16-21 P-value:0.590951833562155 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 16-22 P-value:0.603594817225935 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 16-23 P-value:0.752650416928645 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 17-21 P-value:0.574521265412853 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 17-22 P-value:0.589921955379468 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 17-23 P-value:0.745296227257604 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 18-21 P-value:0.537310759363045 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 18-22 P-value:0.563960751819677 >>ns"
## [1] "mmu-mir-219-2-3p trimming-test 18-23 P-value:0.73213738638379 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 21-26 P-value:0.992979193491957 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 21-27 P-value:0.991920093661686 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 21-28 P-value:0.992195929751444 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 22-26 P-value:0.998037125577552 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 22-27 P-value:0.996319142976593 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 22-28 P-value:0.996306446257224 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 23-26 P-value:0.991809883478799 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 23-27 P-value:0.994565736090519 >>ns"
## [1] "mmu-mir-219-2-3p tailing-test 23-28 P-value:0.994923980723995 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-22.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-23.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-24.png)<!-- -->

```
## [1] "mmu-mir-204-5p trimming-test 16-21 P-value:0.114840601885707 >>ns"
## [1] "mmu-mir-204-5p trimming-test 16-22 P-value:0.681101245236908 >>ns"
## [1] "mmu-mir-204-5p trimming-test 16-23 P-value:0.973702418520574 >>ns"
## [1] "mmu-mir-204-5p trimming-test 17-21 P-value:0.11518200081588 >>ns"
## [1] "mmu-mir-204-5p trimming-test 17-22 P-value:0.680122738317857 >>ns"
## [1] "mmu-mir-204-5p trimming-test 17-23 P-value:0.973465649926833 >>ns"
## [1] "mmu-mir-204-5p trimming-test 18-21 P-value:0.115621732548279 >>ns"
## [1] "mmu-mir-204-5p trimming-test 18-22 P-value:0.678710515175847 >>ns"
## [1] "mmu-mir-204-5p trimming-test 18-23 P-value:0.973133543206085 >>ns"
## [1] "mmu-mir-204-5p tailing-test 21-26 P-value:0.906876675718502 >>ns"
## [1] "mmu-mir-204-5p tailing-test 21-27 P-value:0.900599401468799 >>ns"
## [1] "mmu-mir-204-5p tailing-test 21-28 P-value:0.907001586735046 >>ns"
## [1] "mmu-mir-204-5p tailing-test 22-26 P-value:0.895843068319167 >>ns"
## [1] "mmu-mir-204-5p tailing-test 22-27 P-value:0.887581842889889 >>ns"
## [1] "mmu-mir-204-5p tailing-test 22-28 P-value:0.893793844989741 >>ns"
## [1] "mmu-mir-204-5p tailing-test 23-26 P-value:0.869076561107417 >>ns"
## [1] "mmu-mir-204-5p tailing-test 23-27 P-value:0.856999520411682 >>ns"
## [1] "mmu-mir-204-5p tailing-test 23-28 P-value:0.863371773480411 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-25.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-26.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-27.png)<!-- -->

```
## [1] "mmu-mir-219-2-5p trimming-test 16-21 P-value:0.00552098381687736 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 16-22 P-value:0.00156603301836944 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 16-23 P-value:0.000526439592151446 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 17-21 P-value:0.00317577798368809 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 17-22 P-value:0.000810430195119382 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 17-23 P-value:0.000248434359046735 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 18-21 P-value:0.00321390752279475 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 18-22 P-value:0.000753320629120496 >>*"
## [1] "mmu-mir-219-2-5p trimming-test 18-23 P-value:0.000203140595087985 >>*"
## [1] "mmu-mir-219-2-5p tailing-test 21-26 P-value:0.111143947806822 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 21-27 P-value:0.151189034864951 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 21-28 P-value:0.181159649030453 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 22-26 P-value:0.160319565975624 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 22-27 P-value:0.200282143550049 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 22-28 P-value:0.228900923768157 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 23-26 P-value:0.232517632581453 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 23-27 P-value:0.266919650949807 >>ns"
## [1] "mmu-mir-219-2-5p tailing-test 23-28 P-value:0.291565749783691 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-28.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-29.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-30.png)<!-- -->

```
## [1] "mmu-mir-10b-5p trimming-test 16-21 P-value:0.708901030942516 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 16-22 P-value:0.630229704252306 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 16-23 P-value:0.826969720081803 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 17-21 P-value:0.700246115937416 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 17-22 P-value:0.618280126358844 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 17-23 P-value:0.823999287642605 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 18-21 P-value:0.685600423780837 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 18-22 P-value:0.599686574751033 >>ns"
## [1] "mmu-mir-10b-5p trimming-test 18-23 P-value:0.819463761597534 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 21-26 P-value:0.994571765055719 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 21-27 P-value:0.995044437221879 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 21-28 P-value:0.995337662535843 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 22-26 P-value:0.996017143300739 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 22-27 P-value:0.996278669940583 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 22-28 P-value:0.996450282801838 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 23-26 P-value:0.99890725021003 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 23-27 P-value:0.998867109940752 >>ns"
## [1] "mmu-mir-10b-5p tailing-test 23-28 P-value:0.99887364306616 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-31.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-32.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-33.png)<!-- -->

```
## [1] "mmu-mir-29a-3p trimming-test 16-21 P-value:0.650310237769317 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 16-22 P-value:0.699285860103387 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 16-23 P-value:0.994550124628535 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 17-21 P-value:0.602717223162233 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 17-22 P-value:0.670701562715975 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 17-23 P-value:0.994475666003533 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 18-21 P-value:0.564858603854099 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 18-22 P-value:0.648940130766574 >>ns"
## [1] "mmu-mir-29a-3p trimming-test 18-23 P-value:0.992853898798633 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 21-26 P-value:0.990598526082584 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 21-27 P-value:0.991235682000405 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 21-28 P-value:0.991775176505515 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 22-26 P-value:0.984638756851374 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 22-27 P-value:0.985439715142127 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 22-28 P-value:0.986108493253187 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 23-26 P-value:0.967028734972521 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 23-27 P-value:0.96820012773606 >>ns"
## [1] "mmu-mir-29a-3p tailing-test 23-28 P-value:0.969210112766676 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-34.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-35.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-36.png)<!-- -->

```
## [1] "mmu-mir-29b-1-3p trimming-test 16-21 P-value:0.0471457605204787 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 16-22 P-value:0.0224513888430102 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 16-23 P-value:0.0224980114905447 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 17-21 P-value:0.0367332820677074 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 17-22 P-value:0.0156237676775384 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 17-23 P-value:0.0162111754171668 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 18-21 P-value:0.0226174032737236 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 18-22 P-value:0.00791022416513233 >>*"
## [1] "mmu-mir-29b-1-3p trimming-test 18-23 P-value:0.00921310616314494 >>*"
## [1] "mmu-mir-29b-1-3p tailing-test 21-26 P-value:0.581612254402018 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 21-27 P-value:0.624616600416456 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 21-28 P-value:0.648505984204837 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 22-26 P-value:0.60286308695936 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 22-27 P-value:0.639487059671462 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 22-28 P-value:0.65918068146838 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 23-26 P-value:0.650394027784553 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 23-27 P-value:0.676870625034763 >>ns"
## [1] "mmu-mir-29b-1-3p tailing-test 23-28 P-value:0.690307198013953 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-37.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-38.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-39.png)<!-- -->

```
## [1] "mmu-mir-29c-5p trimming-test 16-21 P-value:0.28840716266688 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 16-22 P-value:0.596382441576874 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 16-23 P-value:0.62668614737863 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 17-21 P-value:0.26746324096747 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 17-22 P-value:0.591150726092031 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 17-23 P-value:0.620765254383842 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 18-21 P-value:0.23130240146071 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 18-22 P-value:0.58327260362377 >>ns"
## [1] "mmu-mir-29c-5p trimming-test 18-23 P-value:0.612251250780343 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 21-26 P-value:0.80054161096952 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 21-27 P-value:0.823695574434276 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 21-28 P-value:0.837995801144766 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 22-26 P-value:0.774182590209042 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 22-27 P-value:0.798583315455231 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 22-28 P-value:0.813416814152328 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 23-26 P-value:0.777970438100553 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 23-27 P-value:0.797674526342598 >>ns"
## [1] "mmu-mir-29c-5p tailing-test 23-28 P-value:0.810037538082939 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-40.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-41.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-42.png)<!-- -->

```
## [1] "mmu-mir-29c-3p trimming-test 16-21 P-value:0.587403722745256 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 16-22 P-value:0.784722457506194 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 16-23 P-value:0.894472623787016 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 17-21 P-value:0.553741539275495 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 17-22 P-value:0.76789158483482 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 17-23 P-value:0.892280078167803 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 18-21 P-value:0.584134649922235 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 18-22 P-value:0.799168700640831 >>ns"
## [1] "mmu-mir-29c-3p trimming-test 18-23 P-value:0.879449930673378 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 21-26 P-value:0.938130338956058 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 21-27 P-value:0.939824085963361 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 21-28 P-value:0.940448627377342 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 22-26 P-value:0.900023264502221 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 22-27 P-value:0.901939008763118 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 22-28 P-value:0.902492225942019 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 23-26 P-value:0.808011642254836 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 23-27 P-value:0.809078616050837 >>ns"
## [1] "mmu-mir-29c-3p tailing-test 23-28 P-value:0.808126889757535 >>ns"
```

![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-43.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-44.png)<!-- -->![](alena_new_data_mar_2017_files/figure-html/unnamed-chunk-21-45.png)<!-- -->

```
## [1] "mmu-mir-29b-2-5p trimming-test 16-21 P-value:0.739770903847523 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 16-22 P-value:0.816664295921607 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 16-23 P-value:0.579637315883919 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 17-21 P-value:0.737706410866572 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 17-22 P-value:0.815340029582888 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 17-23 P-value:0.575604941794864 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 18-21 P-value:0.734384925836335 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 18-22 P-value:0.813389968942076 >>ns"
## [1] "mmu-mir-29b-2-5p trimming-test 18-23 P-value:0.569940438324162 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 21-26 P-value:0.904059325822727 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 21-27 P-value:0.880875951315401 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 21-28 P-value:0.898497484626636 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 22-26 P-value:0.901108816631656 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 22-27 P-value:0.874951010107516 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 22-28 P-value:0.892430311627988 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 23-26 P-value:0.865840304113764 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 23-27 P-value:0.842411671373148 >>ns"
## [1] "mmu-mir-29b-2-5p tailing-test 23-28 P-value:0.862548008566446 >>ns"
```

