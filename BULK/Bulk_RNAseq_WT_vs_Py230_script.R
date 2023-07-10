#______________________
# SCRIPT RNA SEQ BULK
# MSC Ctrl vs Py230
#----------------------

####################### LIBRARIES
library("DESeq2")
library("gdata")
library("dplyr")
library("data.table")
library("ggplot2")
library("ggpubr")


####################### RAW COUNTS MATRIX
#load the raw counts 
dt2 <- read.csv("~/tablecounts_raw.csv", sep = ",")
#load the gene annotations 
Gene.annot <- read.csv("~/tableannot.csv", sep = ",")
Gene.annot = Gene.annot[,2:dim(Gene.annot)[2]]
#put the gene symbols
dt2 = merge(Gene.annot, dt2, by.x = "gene_id", by.y = "X", all = TRUE)
dt2 = dt2[,2:dim(dt2)[2]]
#merge the isoforms
dt2 = as.data.frame(dt2 %>% group_by(gene_name) %>% summarise_each(funs(mean)))
row.names(dt2)=dt2$gene_name
dt2 = dt2[,2:dim(dt2)[2]]
dt2 = round(dt2, digits = 0)
#get the samples metadata
sampleTable1 <- as.data.frame(fread("~/SampleSheet_Clean.csv", header = T))
#Ordering all sample in the same way
dt2 = dt2[,order(colnames(dt2))]
sampleTable1 = sampleTable1[order(sampleTable1$Sample_ID),]
#checks
colnames(dt2)
sampleTable1$Sample_ID
#Replace the column names
colnames(dt2) = sampleTable1$Sample_Name
#add the model info in the metadata
sampleTable1$Model = substr(sampleTable1$Sample_Name, 5, 9)
sampleTable1$Model[sampleTable1$Model=="Ctrl_"] = "WT"

####################### FILTERINGS AND DATASET CLEANNING
mito = which(rownames(dt2) %in% grep(pattern = "^mt-", x = rownames(x = dt2), value = TRUE))
Mirna = which(rownames(dt2) %in% grep(pattern = "mir", x = rownames(x = dt2), value = TRUE))
Ribo1 = which(rownames(dt2) %in% grep(pattern = "Rpl", x = rownames(x = dt2), value = TRUE))
Ribo2 = which(rownames(dt2) %in% grep(pattern = "Rps", x = rownames(x = dt2), value = TRUE))
dt2 = dt2[-c(mito, Mirna, Ribo1, Ribo2),]

####################### CREATE THE DESEQ2 OBJECT
ddsMat <- DESeqDataSetFromMatrix(countData = dt2,
                                 colData = sampleTable1,
                                 design = ~ Model)

dds = ddsMat
#look at the gene distribution
hist(log2(rowMeans(counts(dds))))
plot(density(log2(rowMeans(counts(dds)))))
abline(v=2, col="purple")
#remove the noisy gene expression
LowGenes = names(log2(rowMeans(counts(dds))))[which(log2(rowMeans(counts(dds))) < 2)]
dds <- dds[-which(rownames(counts(dds)) %in% LowGenes), ]

####################### PCA
rld <- rlog(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)

theme_blank = function(){ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) }
color = c("#A9A9A8", "#F39111")
TypVar = rld
ntop = 1500
pcaData <- DESeq2::plotPCA(TypVar, intgroup = c("Model"), ntop = ntop,returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Model, shape=Model)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw(base_line_size = 0,base_size = 15,base_rect_size = 2)+ 
  ggtitle(paste("Number of most variable genes:", ntop)) +
  theme(panel.grid.major= element_line(colour = "white"), panel.grid.minor= element_line(colour = "white"), 
        plot.title = element_text(color="grey30", size=14, face="bold")) + scale_color_manual(values = color,  breaks=c("WT","Py230"))


####################### DIFFERENTIAL GENE EXPRESSION
dds <- DESeq(dds)
#get the results 
res <- results(dds)
res = results(dds, contrast = c("Model","Py230","WT"), cooksCutoff = FALSE, independentFiltering = FALSE)
# Visualize the results
res$name = rownames(res)
ggmaplot(res, main = expression("Control" %->% "Py230"),
         fdr = 0.05, fc = 1.5, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res$name),
         legend = "top", top = 30,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

#get the matrix for the volcano plot
VolMat = as.data.frame(res)[, c(2,6)]




