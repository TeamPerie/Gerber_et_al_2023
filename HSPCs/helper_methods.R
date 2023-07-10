#returns a list of booleans stating TRUE if a given variable is NOT IN the list
"%ni%" <- Negate("%in%")

#this overlays cell density information onto an existing UMAP plot
createDensityPlot <- function(sobj){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "umap"))
  p <- ggplot(tmp.all, aes(x = UMAP_1, y = UMAP_2)) + geom_point(colour="#00000000") + stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=50) + scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) + theme_classic() + 
    theme(legend.position="bottom") + NoAxes()
  
  return(p)
}



geneSignatureScoring <- function(sobj,  gene.sets,gene.set.names,assay){
  
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(sobj@assays$RNA))
    sobj <- AddModuleScore(sobj , features = list(genes.of.interest),
                           name = gene.set.names[i],replace = TRUE,assay = assay)
  }
  
  return(sobj)
  
}



#calculate the background corrected mean expression score of a gene-set using Seurats 
#AddModuleScore method
createGeneSetScore <- function(sobj, gene.set.names, gene.set.fp = "genesets/genesets.csv" ){
  
  gene.sets <- read.csv(gene.set.fp)
  
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(sobj))
    sobj <- AddModuleScore(sobj, features = list(genes.of.interest),name = gene.set.names[i])
  }
  
  return(sobj)
  
}

#generate a volcano plot to visualise DEGs
volcanoPlot <- function(res,numberOfGenesToHighlight){
  res$gene <- rownames(res)
  
  res$sign <- 0
  res$sign[which(res$p_val_adj < 0.05 & res$avg_log2FC > 0.1)] <- 2
  res$sign[which(res$p_val_adj < 0.05 & res$avg_log2FC < -0.1)] <- 1
  res$p_val_adj[res$p_val_adj == 0] <- 1e-312
  
  p <- ggplot(data=res, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=as.factor(sign))) + geom_point( size=1) +
    scale_color_manual(name="", values=c("2" = rgb(67/255,138/255,201/255,1),"1"="#999998", "0"=rgb(220/255,220/255, 220/255,0.2))) +  
    theme(legend.position = "none") + xlim(-0.78,0.78) + 
    xlab("log2 fold change") + ylab("-log10 pvalue") + 
    geom_vline(xintercept=c(-0.1, 0.1), linetype=2) + 
    geom_hline(yintercept=-log10(0.05), linetype=2) 
  
  
  
  p <- p  + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  ) +geom_text_repel(data=res[res$gene %in% c(res$gene[1:numberOfGenesToHighlight],"Mpo","Ctsg"),]
                     , aes(label=gene), size = 4,force = 10)
  
  
  return(p)
  
}


#based on some pre-defined threshold determine cells that are most affected by the PyMT condition
setResponderStatus <- function(lsks,threshold){
  lsks@meta.data$PyMTstatus <- lsks@meta.data$condition
  levels(lsks@meta.data$PyMTstatus) <- c(levels(lsks@meta.data$condition), "PyMT_hi", "PyMT_lo")
  responder <- lsks@meta.data$condition == "PyMT" & lsks@meta.data$PyMT1 > quantile(lsks@meta.data$PyMT1,c(threshold))
  lsks@meta.data$PyMTstatus[responder == TRUE] <- "PyMT_hi"
  lsks@meta.data$PyMTstatus[lsks@meta.data$PyMTstatus == "PyMT"] <- "PyMT_lo"
  return(lsks)
}



#this creates randomised groupings of cells with equivalent numbers to PyMT and WT conditions
#this is needed for permutation testing
generateRandomGroups <- function(sobj,seed){
  
  set.seed(seed * 4342)
  random.index <- sample(1:ncol(sobj),size =  table(sobj@meta.data$condition)[2][[1]]) 
  
  #assign a binary variable to cells in the random.index group or not
  ident <- rep(0,1,ncol(sobj))
  names(ident) <- seq(1,ncol(sobj))
  ident[random.index] <- 1
  
  names(ident) <- colnames(sobj)
  sobj<- AddMetaData(object = sobj, metadata = ident, col.name = "random")
  return(sobj)
}



makeTernaryPlot <- function(sobj,color){
  
  normalize <- function(x){ return((x-min(x))/(max(x)-min(x)))}
  mpp2 <- normalize(sobj@meta.data$MPP2_Pietras1)
  mpp3 <- normalize(sobj@meta.data$MPP3_Pietras1)
  mpp4 <- normalize(sobj@meta.data$MPP4_Pietras1)
  
  df <- data.frame(cbind(mpp2,mpp3,mpp4))
  
  #now normalise the rows
  for(i in 1:nrow(df)){
    df[i,] <- df[i,]/sum(df[i,])
  }
  
  colnames(df) <- c("MPP2","MPP3","MPP4")
  
  triax.plot(df, col.symbols = color,show.grid = TRUE,pch = 1)

  
}

#this method performs a differential expression analysis and returns the numbers of genes that
#are differentially expressed, this is used in permutation testing
performDEA <- function(seed,sobj){
  
  sobj<- generateRandomGroups(sobj,seed)
  #create some random group assignments
  sobj@meta.data$ident <- as.factor(sobj@meta.data$random)
  names(sobj@meta.data$ident) <- colnames(sobj@assays$RNA@data)
  Idents(sobj) <- sobj@meta.data$ident
  
  markers <- FindAllMarkers(sobj, verbose = F,
                            test.use="LR",logfc.threshold = 0.1,
                            only.pos=F,return.thresh = 0.05)
  
  markers.filtered <- markers[markers$p_val_adj < 0.05,]
  
  return(nrow(markers.filtered))
  
}




#load in the cellranger outputs and process to make a seurat object. 
makeSeuratObject <- function(){
  
  pymt1 <- Read10X("rawdata/PyMT1")
  pymt2 <- Read10X("rawdata/PyMT2")
  pymt3 <- Read10X("rawdata/PyMT3")
  pymt4 <- Read10X("rawdata/PyMT4")
  
  #combine all of the data from each mice into one dataframe
  colnames(pymt1) <- paste(colnames(pymt1), "Pymt1", sep = "")
  colnames(pymt2) <- paste(colnames(pymt2), "Pymt2", sep = "")
  colnames(pymt3) <- paste(colnames(pymt3), "Pymt3", sep = "")
  colnames(pymt4) <- paste(colnames(pymt4), "Pymt4", sep = "")
  dataset <- cbind(pymt1,pymt2,pymt3,pymt4)
  
  #let's get rid of pseudogenes, and ribosomal genes, as these are not informative for the rest of the analysis. 
  genes.to.remove <- c(rownames(dataset)[grepl("Gm",rownames(dataset))],
                       rownames(dataset)[grepl("Rik",rownames(dataset))],
                       rownames(dataset)[grepl("Rp",rownames(dataset))])
  
  
  dataset <- as.matrix(dataset)
  dataset <- dataset[rownames(dataset) %ni% genes.to.remove,]
  
  
  
  #create a variable that says which mouse each cell comes from
  mouse    <- c(rep("pymt1",ncol(pymt1)), rep("pymt2",ncol(pymt2)),
                rep("pymt3",ncol(pymt3)),rep("pymt4",ncol(pymt4)))
  
  condition <- c(rep("WT",ncol(pymt1)   + ncol(pymt2)), 
                 rep("PyMT",ncol(pymt3) + ncol(pymt4)))
  
  cell_anns <- data.frame(mouse = mouse, condition = condition)
  rownames(cell_anns) <- make.unique(colnames(dataset))
  
  
  #now we have the dataset and the metadata we can make the seurat object
  lsks <- CreateSeuratObject(counts = dataset, min.cells =25, min.features =1000,meta.data = cell_anns, project = "PyMT")
  
  return(lsks)
  
  
}


# lets assign cell cycle phase using the cell cycle markers from Tirosh et al, 2015
# the method is in the scran package and was from John Marionis/Florian Buettner in 
# the Scialdone (2015) paper.  
runCellCycleAnnotation <- function(sobj){
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  ensembl  <-  mapIds(org.Mm.eg.db, keys=rownames(sobj@assays$RNA), keytype="SYMBOL", 
                      column="ENSEMBL")
  
  # assign cell cycle phase and store it in the metadata of the seurat object
  assignments <- cyclone(sobj@assays$RNA@data, mm.pairs, gene.names=ensembl)
  sobj@meta.data$phases <- assignments$phases
  sobj@meta.data$G1_score <- assignments$normalized.scores$G1
  sobj@meta.data$S_score <- assignments$normalized.scores$S
  sobj@meta.data$G2M_score <- assignments$normalized.scores$G2M
  table(sobj@meta.data$phases)
  return(sobj)
  
}