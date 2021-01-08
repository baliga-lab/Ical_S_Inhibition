
DESeq_analysis_v2 <- function(directory="/Volumes/omics4tb/sturkarslan/icalvum_assembly/results_c5/htseq-counts/",
                              metadata="/Volumes/omics4tb/sturkarslan/icalvum_assembly/icalvum_metadata_v3.txt"){
  
  library('DESeq2');library("RColorBrewer"); library("gplots");library('ggplot2');library("genefilter");library('pheatmap');library(gridExtra);library('apeglm');library(calibrate)
  # data directory
  directory <- directory
  # comparisons file
  comparisons <- read.delim("/Volumes/omics4tb/sturkarslan/icalvum_assembly/final_DET_comparisons.txt", sep="\t", header=F)
  # metadata
  meta <- read.delim(metadata, sep="\t", header=T)
  
  # name of the input
  input.names <- sapply(meta$Sample, function(i) paste(sprintf("%02d", i), "_htseqcounts.txt", sep = ""))
  # input files
  input.files <- sapply(input.names, function(j) list.files(directory, pattern=j))
  # samples
  input.samples <- meta$Sample
  # Conditions
  input.conditions <- meta$Condition
  # Growth Stage
  input.stage <- meta$Stage
  #Replicate
  input.replicate <- meta$Replicate
  
  # create analysis table
  input.table = data.frame(sampleName=input.samples, 
                           filename=input.files, 
                           replicate= input.replicate, 
                           condition=input.conditions, 
                           stage=input.stage)
  input.table$conditionfull <- apply(input.table, 1, function(i) as.factor(paste(i["condition"], i["stage"], sep = ".")))
  
  # collect reads from htseq counts with appropriate design
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = input.table, 
                                    directory = directory,
                                    design = ~ conditionfull)
  
  #### Data visualization plots
  dds.1 <- DESeq(dds)
  comparisons.1 <- resultsNames(dds.1)
  comparisons.1 <- comparisons.1[-1]
  res.1 <- results(dds.1, alpha=0.05)
  
  # normalization
  vsd.1 <- vst(dds.1, blind=FALSE)
  
  df <- as.data.frame(colData(dds.1)[,c("conditionfull","stage")])
  rownames(df) <- paste(vsd.1$condition,vsd.1$stage,vsd.1$replicate, sep = "_")
  sampleDists <- dist(t(assay(vsd.1)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.1$condition,vsd.1$stage,vsd.1$replicate, sep = "_")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, annotation_row=df)
  
  
  # PCA Plot
  pcaData <- plotPCA(vsd.1, intgroup=c("condition", "stage"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  q <- ggplot(pcaData, aes(PC1, PC2, color=stage, shape=condition)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme_gray() + coord_fixed()
  print(q)
  
  # heatmap
  select <- order(rowMeans(counts(dds.1,normalized=TRUE)),
                  decreasing=TRUE)[1:100]
  df.1 <- as.data.frame(colData(dds.1)[,c("condition","stage")])
  pheatmap(assay(vsd.1)[select,], cluster_rows=TRUE, show_rownames=FALSE,
           cluster_cols=TRUE, annotation_col=df.1)
  
  
  
  ## Loop through each commbination of conditions to do DE analysis and produce plots
  pdf(file = "/Volumes/omics4tb/sturkarslan/icalvum_assembly/DEG_Analysis_after_label_switch/icalvum_DE_analysis_with_stages_labelswitch.pdf")
  results <- data.frame()
  for(ref.condition in unique(input.table$conditionfull)){
    for(test.condition in unique(input.table$conditionfull)){
      # create a subset table for the specific comparison
      if(ref.condition == test.condition){
        next
      }
      subset.table <- input.table[which(input.table$conditionfull== ref.condition | input.table$conditionfull== test.condition),]
      # collect reads from htseq counts with appropriate design
      dds <- DESeqDataSetFromHTSeqCount(sampleTable = subset.table, 
                                        directory = directory,
                                        design = ~ conditionfull)
      
      # set the factor level for reference
      dds$conditionfull <- relevel(dds$conditionfull, ref = ref.condition)
      # run DESeq
      dds <- DESeq(dds)
      comparisons <- resultsNames(dds)
      # results names
      comparisons <- comparisons[-1]
      # normalization
      vsd <- vst(dds, blind=FALSE)
      
      # heatmap
      select <- order(rowMeans(counts(dds,normalized=TRUE)),
                      decreasing=TRUE)[1:200]
      df <- as.data.frame(colData(dds)[,c("condition","stage")])
      pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
               cluster_cols=FALSE, show_colnames = T, annotation_col=df, main = paste(comparisons))
      # PCA Plot
      pcaData <- plotPCA(vsd, intgroup=c("condition", "stage"), returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      q <- ggplot(pcaData, aes(PC1, PC2, color=stage, shape=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        theme_gray() + labs(title=paste(comparisons)) + coord_fixed()
      print(q)
      
      ## Results table
      res <- results(dds, name = comparisons, alpha = 0.05)
      #order based p-value
      resOrdered <- res[order(res$pvalue),]
      # LFC shrinkage for plotting
      #coef <- grep(comparisons, comparisons.3)
      resLFC <- lfcShrink(dds, coef= comparisons, res=res, type = "apeglm")
      scatterplot <- plotMA(resLFC, main=paste(comparisons))
      d <- plotCounts(dds, gene=which.min(res$padj), intgroup= c("condition", "stage"), 
                      returnData=TRUE)
      p <- ggplot(d, aes(x=condition, y=count, color=condition, shape=stage)) + 
        geom_point(position=position_jitter(w=0.1,h=0)) + theme_gray() +
        labs(title=row.names(as.data.frame(res)[which.min(res$padj),]))
      print(p)
      
      # volcano plot
      resPlot <- as.data.frame(resLFC)
      resPlot$Gene <- rownames(resPlot)
      maxim <- resPlot$log2FoldChange[order(abs(resPlot$log2FoldChange), decreasing = T)][1]
      with(resPlot, plot(log2FoldChange, -log10(padj), pch=20, main=paste(comparisons), sub=paste(dim(subset(resPlot, padj<.01 & abs(log2FoldChange)>2))[1], "genes with log2 fold.change > 2 & padj < .01"), xlim=c(-maxim,maxim), col="gray"))
      #minim <- res$log2FoldChange[order(res$log2FoldChange)][1]
      
      # Add colored points: red if padj<0.01, orange of log2FC>5, green if both)
      with(subset(resPlot, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="#2c7bb6"))
      with(subset(resPlot, abs(log2FoldChange) >2), points(log2FoldChange, -log10(padj), pch=20, col="#fdae61"))
      with(subset(resPlot, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="#d7191c"))
      
      # Label points with the textxy function from the calibrate plot
      with(subset(resPlot, padj<.01 & abs(log2FoldChange)>5), textxy(log2FoldChange, -log10(padj), labs=Gene, cex=.8))
      abline(v=c(-2, 2), col="gray", lty=2)
      
      tmp.3 <- as.data.frame(resOrdered)
      tmp.3 <- cbind(transcriptName=rownames(tmp.3), tmp.3)
      tmp.3 <-cbind(tmp.3, comparison=comparisons)
      results <- rbind(results, tmp.3)
      
    }
  }
  
  dev.off()
  
  ## Integrate results table with annotations
  library('progress')
  C5.annotations <- read.delim("/Volumes/omics4tb/share/C5/functionalAnnotation/C5.annotation.v2.txt", header=T, sep="\t")
  
  results$condition1 <- sapply(results$comparison, function(i) sub("conditionfull_", "", strsplit(as.character(i), split = "_vs_")[[1]][1]))
  results$condition2 <- sapply(results$comparison, function(i) sub("conditionfull_", "", strsplit(as.character(i), split = "_vs_")[[1]][2]))
  
  # pb <- progress_bar$new(
  #   format = " :what Remaining [:bar] :percent eta: :eta (:spin)",
  #   total = length(results$baseMean), clear = FALSE, width= 100)
  # 
  # results.annotated <- data.frame()
  # for(myrow in 1:length(results$baseMean)){
  #   pb$tick(tokens = list(what = paste("Sim # ", myrow, sep="")))
  #   Sys.sleep(1 / length(results$baseMean))
  #   
  #   name <- rownames(results[myrow,])
  #   annotation.line <- C5.annotations[which(C5.annotations$X.C5.transcript.name == name),]
  #   if(dim(annotation.line)[1] == 0){
  #     annotation.line <-  t(data.frame(myrow=sapply(annotation.line, function(i) paste("NA"))))
  #     results.annotated <- rbind(results.annotated, cbind(annotation.line, 
  #                                                         results[myrow,]))
  #   } else{
  #     results.annotated <- rbind(results.annotated, cbind(annotation.line, 
  #                                                         results[myrow,]))
  #   }
  # }
  
  #results$names <- rownames(results)
  results.annotated <- merge(results, C5.annotations, by.x="transcriptName", by.y="X.C5.transcript.name", all.x=T, all.y=F)
  results.annotated <- results.annotated[order(results.annotated$comparison),]
  results.filtered <- subset(results.annotated, padj<.01 & abs(log2FoldChange)>=2)
  results.filtered <- results.filtered[order(results.filtered$comparison),]
  
  write.table(results.annotated, file="/Volumes/omics4tb/sturkarslan/icalvum_assembly/DEG_Analysis_after_label_switch/icalvum_DE_analysis_with_stages_notfiltered.txt", sep="\t", row.names=F)
  write.table(results.filtered, file="/Volumes/omics4tb/sturkarslan/icalvum_assembly/DEG_Analysis_after_label_switch/icalvum_DE_analysis_with_stages_filtered.txt", sep="\t", row.names=F)
  
}





