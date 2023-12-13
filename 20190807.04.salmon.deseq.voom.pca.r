## salmon results
setwd('~/cowork/derek/CICPT_1979_RNAseq_July2019')
salmon <- paste0(list.dirs('salmongc',recursive=F,full=T),'/quant.sf')[c(-5,-6,-8,-9,-16)]
names <- gsub('(salmongc/)|(-JULY)|(/quant.sf)','',salmon)
cells <- gsub('-[1-3]','',names)
cells <- gsub('SR-4835','SR4835',cells)
cells <- gsub('-','.',cells)
cells1 <- cells
cells1 <- gsub('l[1-2]','',cells1)
meta <- data.frame(cell1=cells,name=names,cell=cells1)

## tx2gene
gtf <- read.table('~/genomes/GRCh38/gencode.v28.primary_assembly.annotation.gtf',sep='\t',comment.char='#',stringsAsFactors=F)
gtftx <- gtf[gtf$V3=='transcript',]
ss <- strsplit(gtftx[,9],";")
tmp <- unlist(lapply(ss,function(x) x[grep("transcript_id",x)]))
tx <- gsub("\"","",gsub(" transcript_id ","",tmp))
tmp <- unlist(lapply(ss,function(x) x[grep("gene_id",x)]))
gene <- gsub("\"","",gsub("gene_id ","",tmp))
tmp2 <- unlist(lapply(ss,function(x) x[grep("gene_name",x)]))
genename<- gsub("\"","",gsub(" gene_name ","",tmp2))
tmp3 <- unlist(lapply(ss,function(x) x[grep("gene_type",x)]))
genetype<- gsub("\"","",gsub(" gene_type ","",tmp3))
tx2gene <- data.frame(tx=tx,gene=gene,stringsAsFactors=F)
## genes
gtfg <- gtf[gtf$V3=='gene',]
ss <- strsplit(gtfg[,9],";")
tmpgeneid <- lapply(ss,function(x) gsub("\"","",gsub("gene_id ","",x[grep("gene_id",x)])))
tmpgenename <- lapply(ss,function(x) gsub("\"","",gsub(" gene_name ","",x[grep("gene_name",x)])))
tmpgenetype <- lapply(ss,function(x) gsub("\"","",gsub(" gene_type ","",x[grep("gene_type",x)])))
genes <- data.frame(t(sapply(seq_along(tmpgeneid),function(i)
                             c(tmpgeneid[[i]],tmpgenename[[i]],tmpgenetype[[i]]))),
                    gtfg[,c(1,4,5,7)],stringsAsFactors=F)
colnames(genes) <- c('gene_id','gene_name','gene_type','chrom','start','end','strand')
save(genes,tx2gene,file='~/genomes/GRCh38/GRCh38.genes.tx.rda')

## deseq2
library(tximport)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
txi <- tximport(salmon, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi,colData = data.frame(meta),design = ~ cell)
ddsde <- DESeq(dds)
# nm <- assays(ddsde)[["avgTxLength"]]
# estimateSizeFactorsForMatrix(counts(ddsde) / nm)
vsd <- assay(rlog(ddsde, blind=FALSE))

## voom
library(limma)
library(edgeR)
cell <- meta$cell
dat0 <- assay(dds)
#dat0 = dat0[genes[match(rownames(dat0),genes[,1]),3]=='protein_coding',]
dge <- DGEList(counts=dat0)
dge <- calcNormFactors(dge)
vsd <- voom(dge,design=model.matrix(~ cell))$E

## pca & correlation
library(corpcor)
svd_tbl <- wt.scale(t(vsd),center=TRUE,scale=TRUE)
svd_run <- fast.svd(svd_tbl)
ds <- round(svd_run$d^2/sum(svd_run$d^2)*100,2)
par(mfrow=c(1,2),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
plot(svd_run$u[,1],svd_run$u[,2],col=as.integer(as.factor(meta[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC2: ',ds[2],'%'),main='')
text(svd_run$u[,1],svd_run$u[,2],meta[,2],cex=0.6,pos=2)
legend('topright',levels(as.factor(meta[,1])),col=as.integer(as.factor(levels(as.factor(meta[,1]))))+1,pch=20,bty='n',cex=0.6)
plot(svd_run$u[,1],svd_run$u[,3],col=as.integer(as.factor(meta[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC3: ',ds[3],'%'),main='')
text(svd_run$u[,1],svd_run$u[,3],meta[,2],cex=0.6,pos=2)
par(mfrow=c(4,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:4){
    for(j in 1:3){
        x <- (i-1)*3+j
        y <- (i-1)*3+j+1
        if(y > i*3) y <- y-3
        cat(x,y,'\n')
        plot(vsd[,x],vsd[,y],pch=20,col=rgb(0,0,0,alpha=0.1),xlab=meta$name[x],ylab=meta$name[y])
        #plot(log2(assay(dds)[,x]+1),log2(assay(dds)[,y]+1),pch=20,col=rgb(0,0,0,alpha=0.1),xlab=meta$name[x],ylab=meta$name[y])
        abline(a=0,b=1)
    }
}
par(mfrow=c(4,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:12){
    hist(vsd[,i],n=100,xlab=meta$name[i],main='')
}

## results
par(mfrow=c(2,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
dat <- data.frame(round(vsd,3),stringsAsFactors=F)
colnames(dat) <- meta$name
comp <- cbind(c('NRG1','SR4835','SR4835.NRG1',"SR4835.NRG1","SR4835.NRG1","SR4835"),c("DMSO","DMSO","DMSO","NRG1","SR4835","NRG1"))
for(i in seq_len(nrow(comp))){
    tmp <- lfcShrink(ddsde, contrast=c("cell",comp[i,1],comp[i,2]))
    DESeq2::plotMA(tmp, ylim=c(-5,5),main=paste0(comp[i,1],'-',comp[i,2]))
    tmp <- data.frame(tmp[,c(1,2,5,6)],stringsAsFactors=F)
    colnames(tmp) <- paste0(comp[i,1],'-',comp[i,2],c('.mean','.lfc','.pv','.fdr'))
    dat <- cbind(dat,tmp[match(rownames(dat),rownames(tmp)),])
}
out <- data.frame(genes[match(rownames(dat),genes$gene_id),],dat)
write.table(out,file='salmongc/CICPT_1979_RNAseq_July2019-srnrg3.txt',sep='\t',row.names=F,quote=F)

idx <- genes$gene_type[match(rownames(dat),genes$gene_id)]=='protein_coding' ## only protein coding genes
plot(log2(dat[idx,4]+1),log2(dat[idx,5]+1),pch=20,col=rgb(0,0,0,alpha=0.1))
