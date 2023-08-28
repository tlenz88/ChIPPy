library(BiocParallel)
library(DiffBind)
library(tidyverse)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

####################################################################################################
######################################## Reading peaksets ##########################################
####################################################################################################

samples <- read.csv("chip_metadata.txt", sep="\t")
chip_df <- dba(sampleSheet=samples)

####################################################################################################
################################### Generate consensus peaksets ####################################
####################################################################################################

chip_consensus <- dba.peakset(chip_df, consensus=c(DBA_CONDITION), minOverlap=1)

####################################################################################################
################################### Generate consensus peaksets ####################################
####################################################################################################

chip_df <- dba.count(chip_df, bUseSummarizeOverlaps=TRUE, bParallel=TRUE)

####################################################################################################
############################# Establishing a model design and contrast #############################
####################################################################################################

chip_df <- dba.contrast(chip_df, reorderMeta=list(Condition="control"), minMembers=2)

####################################################################################################
############################### Performing the differential analysis ###############################
####################################################################################################

chip_df <- dba.analyze(chip_df, method=DBA_DESEQ2, bParallel=TRUE)

####################################################################################################
############################ Retrieving the differentially bound sites #############################
####################################################################################################

chip_df.report <- dba.report(chip_df, method=DBA_DESEQ2, th=0.05, bUsePval=FALSE, fold=1, 
                             bNormalized=TRUE, bFlip=FALSE, precision=0, bCalled=TRUE, bCounts=TRUE)
out <- as.data.frame(chip_df.report)
write.table(out, file="diff_peaks.txt", sep="\t", quote=F, row.names=F)

####################################################################################################
############################################# Plotting #############################################
####################################################################################################

pdf("figures/DiffBind_figures.pdf")

# Correlation Heatmap
dba.plotHeatmap(chip_df, correlations=TRUE, report=chip_df.report, colScheme="Reds")

# PCA plot
#dba.plotPCA(chip_df, contrast=1, th=0.05, report=chip_df.report, vColors=c("#f8766d", "#00bfc4"))

# Venn diagram
dba.plotVenn(chip_consensus, chip_consensus$masks$Consensus)

# MA plot
#debug(DiffBind:::pv.DBAplotMA)
#dba.plotMA(chip_df, th=0.05, fold=1, bNormalized=TRUE, dotSize=1, bSmooth=FALSE, bXY=FALSE)
#dba.plotMA(chip_df, th=0.05, fold=1, bNormalized=TRUE, dotSize=1, bSmooth=FALSE, bXY=TRUE)

# Boxplot
dba.plotBox(chip_df, notch=FALSE, vColors=c("#f8766d", "#00bfc4"))

# Volcano plot
# debug(DiffBind:::pv.DBAplotVolcano)
# When you get the line line that says "plot(p)", enter the following:
#p + geom_point(aes(col = Legend), size = 1) + 
#    xlim(-7.6, 7.6) + 
#    labs(title = "NF54CSAh vs NF54", caption = "FDR <= 0.05, FC >= 2") +
#    ylab(expression(-log[10]~italic(FDR))) + 
#    xlab(expression(log[2]~fold~change)) + 
#    scale_color_manual(values=c("grey60", "black"), labels = c("fail log2FC or FDR", "pass log2FC and FDR")) + 
#    geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "black", size = 0.4) + 
#    geom_vline(xintercept = log2(2), linetype = "longdash", color = "black", size = 0.4) + 
#    geom_vline(xintercept = -log2(2), linetype = "longdash", color = "black", size = 0.4) + 
#    theme(plot.title = element_text(size = 18, face = "bold", family = "serif", margin = margin(7, 0, 7, 0, "pt")),
#          axis.text.y = element_text(size = 10, family = "serif", margin = margin(0, 5, 0, 0, "pt")),
#          axis.text.x = element_text(size = 10, family = "serif", margin = margin(5, 0, 0, 0, "pt")),
#          axis.title.y = element_text(size = 10, family = "serif", margin = margin(0, 7, 0, 7, "pt")),
#          axis.title.x = element_text(size = 10, family = "serif", margin = margin(7, 0, 0, 0, "pt")),
#          axis.line = element_line(color = "black", size = 1),
#          axis.ticks = element_line(color = "black", linewidth = 1), 
#          axis.ticks.length = unit(4, "pt"),
#          panel.grid = element_line(color = "grey92", size = 1),
#          panel.background = element_blank(),
#          legend.position = "top",
#          legend.title = element_blank(),
#          legend.text = element_text(color = "black", size = 10, family = "serif"),
#          legend.background = element_blank(),
#          legend.key = element_blank(),
#          legend.margin = margin(14, 0, 18, 0, "pt"),
#          plot.caption = element_text(size = 10, family = "serif", margin = margin(14, 10, 7, 0, "pt")),
#          plot.margin = margin(7, 14, 7, 7, "pt"))
# then save as PDF size 7" x 7"
dba.plotVolcano(chip_df, th=0.05, fold=1)

# Heatmaps
hmap <- colorRampPalette(c("blue", "white", "red"))(n = 15)
dba.plotHeatmap(chip_df, contrast=1, correlations=FALSE, scale="row", 
                report=chip_df.report, colScheme=hmap)

dev.off()
