##############################
#                            #
#     MATQ-seq - 2023        #
#                            #
##############################


#Required libraries
library(edgeR)
library(rtracklayer)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(reshape2)
library(Cairo)
library(scran)
library(DESeq2)
library(PCAtools)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(gplots)
library(cowplot)
library(BiocParallel)



#--
#load in the counts
#--
setwd("exported_counts/")
counts_all <- readDGE(list.files(pattern = ".count$"), columns = c(1,2), check.names = FALSE,  sep="\t", header=FALSE)

#--
#load in GFF
#--
my_tags <- c("Name","locus_tag")
my_cols <- c("type","start","end")
my_filter <- list(type=c("CDS","sRNA","tRNA","rRNA","pseudogene","5UTR","3UTR"))
gff <- readGFF("salmonella_sl1344.gff", tags=my_tags, columns = my_cols, filter = my_filter)

#--
# Add genes and descriptions to count matrix
#--
counts_all_df <- as.data.frame(counts_all$counts)
counts_all_df$gene_id <- row.names(counts_all_df)
counts_all_df2 <- merge(x=counts_all_df,y=gff,by.x="gene_id",by.y="Name",all.x=T,sort=F)



#----
#
# Cell pheno data
#
#----

cell_pheno_data = readr::read_delim("cell_pheno_data.tsv")

#split into the four conditions
ID_0.1 = cell_pheno_data[cell_pheno_data$Condition == "OD 0.1",]$`Cell ID`
ID_0.3 = cell_pheno_data[cell_pheno_data$Condition == "OD 0.3",]$`Cell ID`
ID_1.0 = cell_pheno_data[cell_pheno_data$Condition == "OD 1.0",]$`Cell ID`
ID_2.0 = cell_pheno_data[cell_pheno_data$Condition == "OD 2.0",]$`Cell ID`




#------
#
# Identifying outliers
#
#------

#-----
#count number of detected genes to work out outliers
#-----

counts_det_features = counts_all_df2[which(!counts_all_df2$type %in% c("3UTR","5UTR")),] 
dim(counts_det_features)
det_genes <- colSums(counts_det_features[,c(2:385)] > 5)


#------
# Apply cutoffs to identify outliers
#------

det_genes_df = data.frame("cell" = names(det_genes),"det_genes" = unname(det_genes))

#0.1
det_genes_df__01 = det_genes_df[det_genes_df$cell %in% ID_0.1,]
limit_up_01 = mean(det_genes_df__01$det_genes) + 2*(sd(det_genes_df__01$det_genes))
limit_down_01 = mean(det_genes_df__01$det_genes) - 2*(sd(det_genes_df__01$det_genes))
outliers_01 = det_genes_df__01[which(det_genes_df__01$det_genes < limit_down_01 | det_genes_df__01$det_genes > limit_up_01),]
det_genes_df__01[which(det_genes_df__01$det_genes < limit_down_01 | det_genes_df__01$det_genes > limit_up_01),]

#0.3
det_genes_df__03 = det_genes_df[det_genes_df$cell %in% ID_0.3,]
limit_up_03 = mean(det_genes_df__03$det_genes) + 2*(sd(det_genes_df__03$det_genes))
limit_down_03 = mean(det_genes_df__03$det_genes) - 2*(sd(det_genes_df__03$det_genes))
outliers_03 = det_genes_df__03[which(det_genes_df__03$det_genes < limit_down_03 | det_genes_df__03$det_genes > limit_up_03),]
det_genes_df__03[which(det_genes_df__03$det_genes < limit_down_03 | det_genes_df__03$det_genes > limit_up_03),]

#1.0
det_genes_df__1 = det_genes_df[det_genes_df$cell %in% ID_1.0,]
limit_up_1 = mean(det_genes_df__1$det_genes) + 2*(sd(det_genes_df__1$det_genes))
limit_down_1 = mean(det_genes_df__1$det_genes) - 2*(sd(det_genes_df__1$det_genes))
outliers_1 = det_genes_df__1[which(det_genes_df__1$det_genes < limit_down_1 | det_genes_df__1$det_genes > limit_up_1),]
det_genes_df__1[which(det_genes_df__1$det_genes < limit_down_1 | det_genes_df__1$det_genes > limit_up_1),] 

#2.0
det_genes_df__2 = det_genes_df[det_genes_df$cell %in% ID_2.0,]
limit_up_2 = mean(det_genes_df__2$det_genes) + 2*(sd(det_genes_df__2$det_genes))
limit_down_2 = mean(det_genes_df__2$det_genes) - 2*(sd(det_genes_df__2$det_genes))
outliers_2 = det_genes_df__2[which(det_genes_df__2$det_genes < limit_down_2 | det_genes_df__2$det_genes > limit_up_2),]
det_genes_df__2[which(det_genes_df__2$det_genes < limit_down_2 | det_genes_df__2$det_genes > limit_up_2),]  

#All outliers
outlier_cells = c(outliers_01$cell, outliers_03$cell, outliers_1$cell, outliers_2$cell)
length(outlier_cells)

counts_all_df_final = counts_all_df2[,!colnames(counts_all_df2) %in% outlier_cells]
counts_all_df_final$length = counts_all_df_final$end -  counts_all_df_final$start


#Conditions with outliers removed
ID_0.1_v2 = ID_0.1[!ID_0.1 %in% outliers_01$cell]
ID_0.3_v2 = ID_0.3[!ID_0.3 %in% outliers_03$cell]
ID_1.0_v2 = ID_1.0[!ID_1.0 %in% outliers_1$cell]
ID_2.0_v2 = ID_2.0[!ID_2.0 %in% outliers_2$cell]




#----
#
# Save down separate counts per condition
#
#----


save_counts_01 = counts_all_df2[,c(ID_0.1_v2,"type","locus_tag","gene_id")]
save_counts_03 = counts_all_df2[,c(ID_0.3_v2,"type","locus_tag","gene_id")]
save_counts_1 = counts_all_df2[,c(ID_1.0_v2,"type","locus_tag","gene_id")]
save_counts_2 = counts_all_df2[,c(ID_2.0_v2,"type","locus_tag","gene_id")]





#-----
#
# Load in existing MATQ-seq data
#
#-----

#note
#Counts have been re-analysed with updated bioinformatics pipeline

setwd("Published_counts/")
counts_all_prev <- readDGE(list.files(pattern = ".count$"), columns = c(1,2), check.names = FALSE,  sep="\t", header=FALSE)

#load in the pheno data
pheno_prev = readr::read_delim("../prev_pheno_data.tsv")

#separate based upon condition and single cells
pheno_prev_sc_ana = pheno_prev[which(pheno_prev$No.cells == 1 & pheno_prev$Condition %in% "Anerobic_shock"),]
pheno_prev_sc_stat = pheno_prev[which(pheno_prev$No.cells == 1 & pheno_prev$Condition %in% "Stationary_phase"),]
pheno_prev_sc_nacl = pheno_prev[which(pheno_prev$No.cells == 1 & pheno_prev$Condition %in% "NACL_shock"),]

#make df
counts_prev_df = counts_all_prev$counts
counts_prev_df = as.data.frame(counts_prev_df)
counts_prev_df$gene = rownames(counts_prev_df)
dim(counts_prev_df)

#add in gff info
counts_prev_df2 = merge(x=counts_prev_df, y=gff , by.x="gene", by.y="Name", all.x=T, sort = F)
counts_prev_df3 = counts_prev_df2[!counts_prev_df2$type %in% c("5UTR","3UTR"),]

#re-analyse the number of detected genes and assigned reads
features_old_ana =  data.frame("ID" = names(colSums(counts_prev_df3[,pheno_prev_sc_ana$ID])),
                               "Colsums" = colSums(counts_prev_df3[,pheno_prev_sc_ana$ID]),
                               "Det_genes" = colSums(counts_prev_df3[,pheno_prev_sc_ana$ID] >5))

features_old_stat =  data.frame("ID" = names(colSums(counts_prev_df3[,pheno_prev_sc_stat$ID])),
                                "Colsums" = colSums(counts_prev_df3[,pheno_prev_sc_stat$ID]),
                                "Det_genes" = colSums(counts_prev_df3[,pheno_prev_sc_stat$ID] >5))

features_old_nacl =  data.frame("ID" = names(colSums(counts_prev_df3[,pheno_prev_sc_nacl$ID])),
                                "Colsums" = colSums(counts_prev_df3[,pheno_prev_sc_nacl$ID]),
                                "Det_genes" = colSums(counts_prev_df3[,pheno_prev_sc_nacl$ID] >5))

#add in sequenced reads
features_old_ana2 = merge(features_old_ana, y=pheno_prev, by.x="ID", by.y="ID", all.x=T, sort=F)
features_old_stat2 = merge(features_old_stat, y=pheno_prev, by.x="ID", by.y="ID", all.x=T, sort=F)
features_old_nacl2 = merge(features_old_nacl, y=pheno_prev, by.x="ID", by.y="ID", all.x=T, sort=F)



#outliers
limit_up_ana = mean(features_old_ana$Det_genes) + 2*(sd(features_old_ana$Det_genes))
limit_down_ana = mean(features_old_ana$Det_genes) - 2*(sd(features_old_ana$Det_genes))
outliers_ana = features_old_ana[which(features_old_ana$Det_genes < limit_down_ana | features_old_ana$Det_genes > limit_up_ana),]
features_old_ana[which(features_old_ana$Det_genes < limit_down_ana | features_old_ana$Det_genes > limit_up_ana),] #0

limit_up_stat = mean(features_old_stat$Det_genes) + 2*(sd(features_old_stat$Det_genes))
limit_down_stat = mean(features_old_stat$Det_genes) - 2*(sd(features_old_stat$Det_genes))
outliers_stat = features_old_stat[which(features_old_stat$Det_genes < limit_down_stat | features_old_stat$Det_genes > limit_up_stat),]
features_old_stat[which(features_old_stat$Det_genes < limit_down_stat | features_old_stat$Det_genes > limit_up_stat),] #0

limit_up_nacl = mean(features_old_nacl$Det_genes) + 2*(sd(features_old_nacl$Det_genes))
limit_down_nacl = mean(features_old_nacl$Det_genes) - 2*(sd(features_old_nacl$Det_genes))
outliers_nacl = features_old_nacl[which(features_old_nacl$Det_genes < limit_down_nacl | features_old_nacl$Det_genes > limit_up_nacl),]
features_old_nacl[which(features_old_nacl$Det_genes < limit_down_nacl | features_old_nacl$Det_genes > limit_up_nacl),] #1

#remove outlier
features_old_nacl2 = features_old_nacl2[!features_old_nacl2$ID %in% outliers_nacl$ID,]


#-------
#
# Comparing between methods (new vs old MATq-seq)
#
#-------



#------
# Detected genes vs sequencing depth
#------


new_seq_reads_v2 = data.frame("Condition" = c(features_0.1$Condition,
                                              features_0.3$Condition,
                                              features_1.0$Condition,
                                              features_2.0$Condition,
                                              features_old_ana2$Condition,
                                              features_old_nacl2$Condition,
                                              features_old_stat2$Condition),
                              "Sequenced_reads" = c(features_0.1$Sequenced_reads,
                                                    features_0.3$Sequenced_reads,
                                                    features_1.0$Sequenced_reads,
                                                    features_2.0$Sequenced_reads, 
                                                    features_old_ana2$fastq_reads,
                                                    features_old_nacl2$fastq_reads,
                                                    features_old_stat2$fastq_reads),
                              "Detected_genes" = c(features_0.1$Det_genes,
                                                   features_0.3$Det_genes,
                                                   features_1.0$Det_genes,
                                                   features_2.0$Det_genes, 
                                                   features_old_ana2$Det_genes,
                                                   features_old_nacl2$Det_genes,
                                                   features_old_stat2$Det_genes),
                              "Assigned_reads" = c(features_0.1$Colsums,
                                                   features_0.3$Colsums,
                                                   features_1.0$Colsums,
                                                   features_2.0$Colsums, 
                                                   features_old_ana2$Colsums,
                                                   features_old_nacl2$Colsums,
                                                   features_old_stat2$Colsums)
)

new_seq_reads_v2$Condition = factor(new_seq_reads_v2$Condition, levels = c("Anerobic_shock","NACL_shock","Stationary_phase","OD 0.1","OD 0.3", "OD 1.0","OD 2.0"))



#------
# Detected genes
#------

cairo_pdf(file="Detected_genes_old_new.pdf", width=12, height=6, bg = "white")
ggplot(new_seq_reads_v2, aes(x=Condition, y=Detected_genes, group=Condition, fill=Condition)) +
  geom_violin(width=0.9) +
  geom_boxplot(width=0.2, alpha=0.2) +
  theme_ipsum() + ylim(0,700) +
  theme(legend.position="none",
    plot.title = element_text(size=11), axis.title.x = element_text(size=9), axis.title.y = element_text(size=9)) +
  ggtitle("Detected genes per condition") +
  xlab("Growth") + ylab("Detected genes")
dev.off()



#------
#
# Looking at the proportion of Zeros per cell 
#
#------

#---
#Load in bulk RNA-seq data to work out which genes have zero counts
#---

#Load in counts

setwd("Bulk_counts/")

bulk_counts <- readDGE(list.files(pattern = ".count$"), columns = c(1,2), check.names = FALSE,  sep="\t", header=FALSE)

#add gene length
bulk_counts_df = as.data.frame(bulk_counts$counts)
bulk_counts_df$gene = row.names(bulk_counts_df)
bulk_counts_df$rowSums = rowSums(bulk_counts_df[,1:12])

bulk_genes_to_exclude = rownames(bulk_counts_df[bulk_counts_df$rowSums <1,])


#remove these genes from count matrix - current
counts_prop_zeros_df = counts_all_df_final[!counts_all_df_final$gene_id %in% bulk_genes_to_exclude,]
counts_prop_zeros_df = counts_prop_zeros_df[!counts_prop_zeros_df$type %in% c("5UTR","3UTR"),]

#remove these genes from count matrix - previous
counts_prop_zeros_prev_df = counts_prev_df3[!counts_prev_df3$gene %in% bulk_genes_to_exclude,]


#work out the number of zero counts per cell
prop_zeros_new = colSums(counts_prop_zeros_df[,2:371] < 1)
prop_zeros_prev = colSums(counts_prev_df3[,2:132] < 1)

prop_zeros_new_df = data.frame("ID" = names(prop_zeros_new), "zeros" = unname(prop_zeros_new))
prop_zeros_prev_df = data.frame("ID" = names(prop_zeros_prev), "zeros" = unname(prop_zeros_prev))

#merge condition info with df
prop_zeros_new_df2 = merge(x=prop_zeros_new_df, by.x="ID", y=new_seq_reads, by.y="ID", all.x=T)
prop_zeros_prev_df2 = merge(x=prop_zeros_prev_df, by.x="ID", y=pheno_prev, by.y="ID", all.x=T)

#filter and remove outliers
prop_zeros_new_df2 = prop_zeros_new_df2[!prop_zeros_new_df2$Outlier == "Y",]
prop_zeros_prev_df2 = prop_zeros_prev_df2[prop_zeros_prev_df2$No.cells == "1",]

#merge
prop_zeros_prev_df2$Sequenced_reads = prop_zeros_prev_df2$fastq_reads
prop_zeros_merged = rbind(prop_zeros_new_df2[,c(3,2,5)],prop_zeros_prev_df2[,c(4,2,6)])

#add factors for plotting
prop_zeros_merged$Condition = factor(prop_zeros_merged$Condition, levels = c("Anerobic_shock","NACL_shock","Stationary_phase","OD 0.1","OD 0.3", "OD 1.0","OD 2.0"))

#add in the percent
prop_zeros_merged$percent_zero = prop_zeros_merged$zeros / dim(counts_prop_zeros_df)[1]



#--
# Plots for zero proportions of genes
#--


#proportions
cairo_pdf(file="zero_prop_old_and_new.pdf", width=10, height=6, bg = "white")
ggplot(prop_zeros_merged, aes(x=Sequenced_reads, y=percent_zero,color=Condition)) + ylim(c(0.78,1.0)) +
  geom_point(size=2) + theme_ipsum() + xlab("Sequenced reads") + ylab("Proportion of zeros")
dev.off()


#--
# Detected genes per condition
#--

percent_det_vct = c(table(rowSums(counts_prop_zeros_df[,2:371] >5) >= 1)[2] / dim(counts_prop_zeros_df)[1],
                    table(rowSums(counts_prop_zeros_df[,ID_0.1_v2] >5) >= 1)[2] / dim(counts_prop_zeros_df)[1],
                    table(rowSums(counts_prop_zeros_df[,ID_0.3_v2] >5) >= 1)[2] / dim(counts_prop_zeros_df)[1],
                    table(rowSums(counts_prop_zeros_df[,ID_1.0_v2] >5) >= 1)[2] / dim(counts_prop_zeros_df)[1],
                    table(rowSums(counts_prop_zeros_df[,ID_2.0_v2] >5) >= 1)[2] / dim(counts_prop_zeros_df)[1],
                    
                    #prev
                    table(rowSums(counts_prop_zeros_prev_df[,c(features_old_ana2$ID,
                                                               features_old_nacl2$ID, 
                                                               features_old_stat2$ID)] > 5) >= 1 )[2] / dim(counts_prop_zeros_prev_df)[1],
                    table(rowSums(counts_prop_zeros_prev_df[,features_old_ana2$ID] > 5) >= 1 )[2] / dim(counts_prop_zeros_prev_df)[1],
                    table(rowSums(counts_prop_zeros_prev_df[,features_old_nacl2$ID] > 5) >= 1 )[2] / dim(counts_prop_zeros_prev_df)[1],
                    table(rowSums(counts_prop_zeros_prev_df[,features_old_stat2$ID] > 5) >= 1 )[2] / dim(counts_prop_zeros_prev_df)[1]
                    
  
)



percent_det_df = data.frame("Condition" = c("Combined_new","OD 0.1","OD 0.3","OD 1.0","OD 2.0", "Combined_old","ana","nacl","lsp"),
                            "Percent" = percent_det_vct
)
percent_det_df$Condition = factor(percent_det_df$Condition, levels = c("lsp","nacl","ana", "Combined_old","OD 2.0","OD 1.0","OD 0.3","OD 0.1","Combined_new"))


cairo_pdf(file="percent_genome_coverage_old_new.pdf", width=6, height=8, bg = "white")
ggplot(percent_det_df, aes(x=Condition, y=Percent))+ geom_bar(stat="identity") + coord_flip() + theme_ipsum() + xlab("") + ylim(c(0,1))
dev.off()



#-----
#
#Determine number of counted features
#
#-----

biotype_summary_new = counts_all_df2[,c(2:386)] %>% group_by(type) %>% summarise_each(funs(sum))
#remove bottom two rows containing all zeros
biotype_summary_new = biotype_summary_new[-c(7,8),]


biotype_summary_prev = counts_prev_df2[,2:133] %>% group_by(type) %>% summarise_each(funs(sum))
#remove bottom two rows containing all zeros
biotype_summary_prev = biotype_summary_prev[-c(7,8),]


percentage_mapped_reads = melt(data.frame("Gene_type" = biotype_summary_new$type,
                                          "OD 0.1" = prop.table(rowMeans(biotype_summary_new[,c(ID_0.1_v2)])),
                                          "OD 0.3" = prop.table(rowMeans(biotype_summary_new[,c(ID_0.3_v2)])),
                                          "OD 1.0" = prop.table(rowMeans(biotype_summary_new[,c(ID_1.0_v2)])),
                                          "OD 2.0" = prop.table(rowMeans(biotype_summary_new[,c(ID_2.0_v2)])),
                                          "Ana" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_ana$ID)])),
                                          "Nacl" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_nacl$ID)])),
                                          "Stat" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_stat$ID)]))
))



percentage_mapped_reads_rRNA = data.frame(row.names = biotype_summary_new$type,
           "OD 0.1" = prop.table(rowMeans(biotype_summary_new[,c(ID_0.1_v2)])),
           "OD 0.3" = prop.table(rowMeans(biotype_summary_new[,c(ID_0.3_v2)])),
           "OD 1.0" = prop.table(rowMeans(biotype_summary_new[,c(ID_1.0_v2)])),
           "OD 2.0" = prop.table(rowMeans(biotype_summary_new[,c(ID_2.0_v2)])),
           "Ana" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_ana$ID)])),
           "Nacl" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_nacl$ID)])),
           "Stat" = prop.table(rowMeans(biotype_summary_prev[,c(pheno_prev_sc_stat$ID)]))
)

percentage_mapped_reads_rRNA = rbind(percentage_mapped_reads_rRNA[1,] + 
                                      percentage_mapped_reads_rRNA[2,] + 
                                      percentage_mapped_reads_rRNA[3,] + 
                                      percentage_mapped_reads_rRNA[5,] + 
                                      percentage_mapped_reads_rRNA[6,],
                                    percentage_mapped_reads_rRNA[4,])

percentage_mapped_reads_rRNA$Biotype = c("Other","rRNA")
percentage_mapped_reads_rRNA = melt(percentage_mapped_reads_rRNA)

percentage_mapped_reads_rRNA$variable = factor(percentage_mapped_reads_rRNA$variable, levels = c("Ana","Nacl","Stat","OD.0.1","OD.0.3","OD.1.0","OD.2.0"))




#------
#
# DESeq2 - Scran - PCA plot
#
#------


counts_all_df_final_v2 = counts_all_df_final[which(!counts_all_df_final$type %in% c("3UTR","5UTR","rRNA","tRNA")),] 

#remove additional outlier from PCA plot
additional_outliers = c("III_D1_S37_R1_001")
counts_all_df_final_v3 = counts_all_df_final_v2[!colnames(counts_all_df_final_v2) %in% additional_outliers]

#prepare the metadata
cell_pheno_data = cell_pheno_data[!cell_pheno_data$`Cell ID` %in% outlier_cells,]
cell_pheno_data = cell_pheno_data[!cell_pheno_data$`Cell ID` %in% additional_outliers,]

#convert to matrix and add in gene names as rownames
counts_all_mtx = as.matrix(counts_all_df_final_v3[,c(2:370)])
#add in the gene names
row.names(counts_all_mtx) = counts_all_df_final_v3$gene_id
counts_all_mtx


#update the pheno data to group conditions
de_pheno = cell_pheno_data
de_pheno = de_pheno[match(colnames(counts_all_mtx),de_pheno$`Cell ID`),]

de_pheno$DE = "0"
de_pheno[de_pheno$Condition == "OD 0.1",]$DE <- "0.1" 
de_pheno[de_pheno$Condition == "OD 0.3",]$DE <- "0.3" 
de_pheno[de_pheno$Condition == "OD 1.0",]$DE <- "1.0" 
de_pheno[de_pheno$Condition == "OD 2.0",]$DE <- "2.0" 

#filter
counts_all_mtx = counts_all_mtx[rowSums(counts_all_mtx > 0) > 5,]


#---
#Deseq2
#---
dds = DESeqDataSetFromMatrix(countData = counts_all_mtx, colData = de_pheno, design = ~DE)
dds_de <- DESeq(dds, test="LRT", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf, reduced=~1, parallel = F)



#---
#add in the scran factors
#---
sce_deseq2 = SingleCellExperiment(assays=list(counts=counts(dds_de)), colData=de_pheno)
clusters_deseq2 <- scran::quickCluster(sce_deseq2, min.size=25)
sce_deseq2 <- scran::computeSumFactors(sce_deseq2, cluster = clusters_deseq2)

#Check to make sure in the in the same order
table(names(sizeFactors(dds_de)) == sce_deseq2$`Cell ID`) 
#add size factors
sizeFactors(dds_de) = sce_deseq2$sizeFactor



#---
# Variance stabilised data and PCA plot
#---


vst <- assay(vst(dds_de), blind=F)

cell_pheno_data_df = as.data.frame(de_pheno)
rownames(cell_pheno_data_df) = cell_pheno_data_df$`Cell ID`

table(rownames(cell_pheno_data_df) == colnames(vst))

pca_vst <- pca(vst, metadata = cell_pheno_data_df, removeVar = 0.1)
pca_cols =  c('OD 0.1' = '#BCBCBB', 'OD 0.3' = '#F3E226', 'OD 1.0' = '#1298CE', 'OD 2.0' = '#F26944')


cairo_pdf(file="PCA_conditions.pdf", width=9, height=6, bg = "white")
biplot(pca_vst, colby = "Condition", showLoadings = F, lab = NULL, legendPosition = 'right', title = "", titleLabSize = 20,
       colkey = pca_cols)
dev.off()


#---
# Explore PCA loadings 
#---

cairo_pdf(file="PCA_conditions_with_loadings.pdf", width=9, height=6, bg = "white")
biplot(pca_vst, colby = "Condition", showLoadings = T, lab = NULL, #ntopLoadings = 10,
       legendPosition = 'right', title = "PCA with loadings", titleLabSize = 20, 
       colkey = pca_cols, fillBoxedLoadings = alpha("white", 1), colLegendTitle = "Condition")
dev.off()




#---
# Making PCA plots with single genes
#---

#add genes
cell_pheno_data_df_with_expression = cell_pheno_data_df

#update to add in gene
vst_with_expression = data.frame(vst)
vst_with_expression$Genes = rownames(vst)


#add expression
expression_values = data.frame(t(vst_with_expression[vst_with_expression$Genes %in%  c("fliC","flaG","aceE","sipC","sipB"),]))
expression_values$ID = row.names(expression_values)

cell_pheno_data_df_with_expression = merge(x=cell_pheno_data_df_with_expression, by.x="Cell ID", y=expression_values, by.y="ID", all.x=T, sort=F)
rownames(cell_pheno_data_df_with_expression) = cell_pheno_data_df_with_expression$`Cell ID`

#make numeric
cell_pheno_data_df_with_expression$sipC = as.numeric(cell_pheno_data_df_with_expression$sipC)
cell_pheno_data_df_with_expression$sipB = as.numeric(cell_pheno_data_df_with_expression$sipB)
cell_pheno_data_df_with_expression$fliC = as.numeric(cell_pheno_data_df_with_expression$fliC)
cell_pheno_data_df_with_expression$flaG = as.numeric(cell_pheno_data_df_with_expression$flaG)
cell_pheno_data_df_with_expression$aceE = as.numeric(cell_pheno_data_df_with_expression$aceE)

#prepare pca object
pca_with_gene_exp <- pca(vst, metadata = cell_pheno_data_df_with_expression, removeVar = 0.1)


#Extract out OD0.1 and 0.3
separate_PCA_0.1_0.3 = vst[,colnames(vst) %in% c(ID_0.1_v2,ID_0.3_v2)]
#to make the other colours grey
separate_PCA_0.1_0.3_pca <- pca(separate_PCA_0.1_0.3, metadata = cell_pheno_data_df_with_expression[cell_pheno_data_df_with_expression$`Cell ID` %in% c(ID_0.1_v2,ID_0.3_v2),], removeVar = 0.1)

separate_PCA_1_2 = vst[,colnames(vst) %in% c(ID_1.0_v2,ID_2.0_v2)]
separate_PCA_1_2_pca <- pca(separate_PCA_1_2, metadata = cell_pheno_data_df_with_expression[cell_pheno_data_df_with_expression$`Cell ID` %in% c(ID_1.0_v2,ID_2.0_v2),], removeVar = 0.1)



cairo_pdf(file="PCA_gene_flaG.pdf", width=7, height=6, bg = "white")
biplot(pca_with_gene_exp, colby = "flaG", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "PCA - flaG expression", titleLabSize = 20) + scale_colour_gradient(low = 'khaki1', high = 'red2', limits = c(4,20))
dev.off()

cairo_pdf(file="PCA_gene_aceE.pdf", width=7, height=6, bg = "white")
biplot(pca_with_gene_exp, colby = "aceE", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "PCA - aceE expression", titleLabSize = 20) + scale_colour_gradient(low = 'khaki1', high = 'red2', limits = c(4,20))
dev.off()

cairo_pdf(file="PCA_gene_sipB.pdf", width=7, height=6, bg = "white")
biplot(pca_with_gene_exp, colby = "sipB", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "PCA - sipB expression", titleLabSize = 20) + scale_colour_gradient(low = 'khaki1', high = 'red2', limits = c(4,20))
dev.off()

cairo_pdf(file="PCA_gene_sipC.pdf", width=7, height=6, bg = "white")
biplot(pca_with_gene_exp, colby = "sipC", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "PCA - sipC expression", titleLabSize = 20) + scale_colour_gradient(low = 'khaki1', high = 'red2', limits = c(4,20))
dev.off()

cairo_pdf(file="PCA_0.1_0.3.pdf", width=7, height=6, bg = "white")
biplot(separate_PCA_0.1_0.3_pca, colby = "Condition", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "OD 0.1 and 0.3", titleLabSize = 20)
dev.off()

cairo_pdf(file="PCA_0.1_0.3_fliC.pdf", width=7, height=6, bg = "white")
biplot(separate_PCA_0.1_0.3_pca, colby = "fliC", showLoadings = F, lab = NULL, legendPosition = 'right', colLegendTitle = "log2(norm.)",
       title = "0.1 0.3 fliC", titleLabSize = 20) + scale_colour_gradient(low = 'khaki1', high = 'red2', limits = c(4,20))
dev.off()







#------
#
# #scran HVG
#
#------


scran_cells = log2(counts(dds_de, normalized=T)+1)
scran_cells = as.data.frame(scran_cells)

scran_cells__01 = scran_cells[,colnames(scran_cells) %in% c(ID_0.1_v2)]
scran_cells__03 = scran_cells[,colnames(scran_cells) %in% c(ID_0.3_v2)]
scran_cells__1 = scran_cells[,colnames(scran_cells) %in% c(ID_1.0_v2)]
scran_cells__2 = scran_cells[,colnames(scran_cells) %in% c(ID_2.0_v2)]


sce_HVG__01 =  scran::modelGeneVar(x = scran_cells__01)
sce_HVG__03 =  scran::modelGeneVar(x = scran_cells__03)
sce_HVG__1 =  scran::modelGeneVar(x = scran_cells__1)
sce_HVG__2 =  scran::modelGeneVar(x = scran_cells__2)


#numbers from modelGeneVar
top.hvgs__01 <- getTopHVGs(sce_HVG__01, prop=0.01) 
top.hvgs__03 <- getTopHVGs(sce_HVG__03, prop=0.01) 
top.hvgs__1 <- getTopHVGs(sce_HVG__1, prop=0.01) 
top.hvgs__2 <- getTopHVGs(sce_HVG__2, prop=0.01) 





#four heatmaps with all cells, but separated into timepoints

#to ensure the cells are in the correct order for the heatmaps, I re-merge each of the separate times into a df
scran_cells_ordered = cbind(scran_cells__01,scran_cells__03,scran_cells__1,scran_cells__2)
cell_pheno_data$Condition = factor(cell_pheno_data$Condition, levels = c("OD 0.1", "OD 0.3", "OD 1.0", "OD 2.0"))

#need to re-order the phenodata to match the matrix
table(colnames(scran_cells_ordered) == cell_pheno_data[match(colnames(scran_cells_ordered), cell_pheno_data$`Cell ID`),]$`Cell ID`)

hvg_phenotypes = cell_pheno_data[match(colnames(scran_cells_ordered), cell_pheno_data$`Cell ID`),]
hm_anno = HeatmapAnnotation(Condition = hvg_phenotypes$Condition)

col_fun = colorRamp2(c(0, 2.5,5,7.5,10,12.5,15,17.5,20), RColorBrewer::brewer.pal(9, "YlOrRd"))
lgd = Legend(col_fun = col_fun, title = "HM")



cairo_pdf(file="HM_all_cells__hvg_0.1.pdf", width=6, height=3, bg = "white")
ComplexHeatmap::Heatmap(as.matrix( scran_cells_ordered[rownames(scran_cells_ordered) %in% top.hvgs__01,]   ), cluster_rows = T, top_annotation = hm_anno,
                        cluster_columns = F, show_column_names = FALSE, col = col_fun, name="log2(norm.)")
dev.off()

cairo_pdf(file="HM_all_cells__hvg_0.3.pdf", width=6, height=3, bg = "white")
ComplexHeatmap::Heatmap(as.matrix( scran_cells_ordered[rownames(scran_cells_ordered) %in% top.hvgs__03,]   ), cluster_rows = T, top_annotation = hm_anno,
                        cluster_columns = F, show_column_names = FALSE, col = col_fun, name="log2(norm.)")
dev.off()

cairo_pdf(file="HM_all_cells__hvg_1.0.pdf", width=6, height=3, bg = "white")
ComplexHeatmap::Heatmap(as.matrix( scran_cells_ordered[rownames(scran_cells_ordered) %in% top.hvgs__1,]   ), cluster_rows = T, top_annotation = hm_anno,
                        cluster_columns = F, show_column_names = FALSE, col = col_fun, name="log2(norm.)")
dev.off()

cairo_pdf(file="HM_all_cells__hvg_2.0.pdf", width=6, height=3, bg = "white")
ComplexHeatmap::Heatmap(as.matrix( scran_cells_ordered[rownames(scran_cells_ordered) %in% top.hvgs__2,]   ), cluster_rows = T, top_annotation = hm_anno,
                        cluster_columns = F, show_column_names = FALSE, col = col_fun, name="log2(norm.)")
dev.off()




#----
#
#TPM normalisation - single cells
#
#----

#from here: https://www.biostars.org/p/335187/
tpm2 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x,na.rm = T)))
}



#exclude 5' and 3'
counts_all_df_final_tpm = counts_all_df_final_tpm[!counts_all_df_final_tpm$type %in% c("3UTR","5UTR"),]
rownames(counts_all_df_final_tpm) = counts_all_df_final_tpm$gene_id

#convert to tpm
counts_all_df_final_tpm2 = tpm2(counts = counts_all_df_final_tpm[,c(2:371)],len = counts_all_df_final_tpm$length)

#check all are 1m
summary(colSums(counts_all_df_final_tpm2))

#add the gene ID back in
counts_all_df_final_tpm2 = as.data.frame(counts_all_df_final_tpm2)
counts_all_df_final_tpm2$ID = row.names(counts_all_df_final_tpm2)

#add the gene info back in
counts_all_df_final_tpm2 = merge(x=counts_all_df_final_tpm2, y=counts_all_df_final_tpm[,c(1,372:376)], by.x="ID",by.y="gene_id", all.x=T)

#each timepoint
tpm_0.1 = counts_all_df_final_tpm2[,c("locus_tag","ID","type",ID_0.1_v2)]
tpm_0.3 = counts_all_df_final_tpm2[,c("locus_tag","ID","type",ID_0.3_v2)]
tpm_1.0 = counts_all_df_final_tpm2[,c("locus_tag","ID","type",ID_1.0_v2)]
tpm_2.0 = counts_all_df_final_tpm2[,c("locus_tag","ID","type",ID_2.0_v2)]







#------
#
# sRNA plots for detected genes and unique sRNA per condition
#
#------

sRNA_counts__0.1 = save_counts_01[save_counts_01$type %in% "sRNA",c(1:95,98)]
sRNA_counts__0.3 = save_counts_03[save_counts_03$type %in% "sRNA",c(1:92,95)]
sRNA_counts__1.0 = save_counts_1[save_counts_1$type %in% "sRNA",c(1:92,95)]
sRNA_counts__2.0 = save_counts_2[save_counts_2$type %in% "sRNA",c(1:91,94)]

numbers_of_sRNA_df = data.frame("Condition" = c("OD0.1","OD0.3","OD1.0","OD2.0"),
                                "Total_det_genes" = c(sum(colSums(sRNA_counts__0.1[,1:95] > 5)),
                                                      sum(colSums(sRNA_counts__0.3[,1:92] > 5)),
                                                      sum(colSums(sRNA_counts__1.0[,1:92] > 5)),
                                                      sum(colSums(sRNA_counts__2.0[,1:91] > 5))
                                ),
                                "Unique_sRNA" = c(sum(rowSums(sRNA_counts__0.1[,1:95] > 5) !=0),
                                                  sum(rowSums(sRNA_counts__0.3[,1:92] > 5) !=0),
                                                  sum(rowSums(sRNA_counts__1.0[,1:92] > 5) !=0),
                                                  sum(rowSums(sRNA_counts__2.0[,1:91] > 5) !=0)
                                )
)


cairo_pdf(file="srna_det_genes.pdf", width=6, height=4, bg = "white")
ggplot(numbers_of_sRNA_df, aes(Condition, Total_det_genes)) + geom_bar(stat="identity") + xlab("") + ylab("Total detected sRNA") + theme_ipsum()
dev.off()

cairo_pdf(file="srna_unique_genes.pdf", width=6, height=4, bg = "white")
ggplot(numbers_of_sRNA_df, aes(Condition, Unique_sRNA)) + geom_bar(stat="identity") + xlab("") + ylab("Unique sRNA") + theme_ipsum()
dev.off()





#------
#
# sRNA heatmap
#
#------

sRNA_TPM = counts_all_df_final_tpm2[counts_all_df_final_tpm2$type %in% "sRNA",]
sRNA_TPM$sum = rowSums(sRNA_TPM[,c(2:371)])
sRNA_TPM = sRNA_TPM[order(sRNA_TPM$sum, decreasing = T),]
#all sRNA ordered by TPM sum
sRNA_TPM_all_ordered = sRNA_TPM[,c(ID_0.1_v2,ID_0.3_v2,ID_1.0_v2,ID_2.0_v2)]


#filtering sRNA
sRNA_TPM_v2 = sRNA_TPM[which(sRNA_TPM$sum > 100),c(1:372)]

#put in order of timepoints (for the heatmap)
sRNA_TPM_v3 = sRNA_TPM_v2[,c(ID_0.1_v2,ID_0.3_v2,ID_1.0_v2,ID_2.0_v2)]

cell_pheno_data_sRNA = cell_pheno_data[match(colnames(sRNA_TPM_v3), cell_pheno_data$`Cell ID`),]

#make HM annotation
ha_srna = HeatmapAnnotation(Condition = cell_pheno_data_sRNA$Condition)

pdf(file="sRNA_heatmap.pdf", width=6, height=10, bg = "white")
ComplexHeatmap::Heatmap(as.matrix(log2(sRNA_TPM_v3 +1)),row_names_gp = gpar(fontsize = 8), 
                        col = RColorBrewer::brewer.pal(9, "YlOrRd"),top_annotation = ha_srna, cluster_columns = F,cluster_rows = T,
                        heatmap_legend_param = list(title = "log2(TPM)"), column_names_gp = gpar(fontsize = 0))

dev.off()











#----
#
#Calculate Pseudo bulk single cells
#
#----


#This isnt the TPM values, its still counts (even though the variable says TPM)
sc_grouped_counts = data.frame(ID_0.1_counts = rowSums(counts_all_df_final_tpm[,c(ID_0.1_v2)]),
                               ID_0.3_counts = rowSums(counts_all_df_final_tpm[,c(ID_0.3_v2)]),
                               ID_1.0_counts = rowSums(counts_all_df_final_tpm[,c(ID_1.0_v2)]),
                               ID_2.0_counts = rowSums(counts_all_df_final_tpm[,c(ID_2.0_v2)])
)

sc_grouped_counts$ID = rownames(sc_grouped_counts)

#add gene length
sc_grouped_counts = merge(x=sc_grouped_counts, y=counts_all_df_final_tpm[,c(1,376)], by.x="ID",by.y="gene_id",all.x=T)
row.names(sc_grouped_counts) = sc_grouped_counts$ID

#convert to TPM
sc_grouped_counts_tpm = tpm2(counts = sc_grouped_counts[,c(2:5)],len = sc_grouped_counts$length)
#check all are 1m
summary(colSums(sc_grouped_counts_tpm))

sc_grouped_counts_tpm = as.data.frame(sc_grouped_counts_tpm)
sc_grouped_counts_tpm$ID = row.names(sc_grouped_counts_tpm)



#---
#
# Convert bulk data to TPM
#
#---

bulk_counts_df = merge(x=bulk_counts_df, y=gff, by.x="gene",by.y="Name", all.x=T)
bulk_counts_df$length = bulk_counts_df$end - bulk_counts_df$start
rownames(bulk_counts_df) = bulk_counts_df$gene

bulk_tpms = tpm2(counts = bulk_counts_df[,c(2:13)],len = bulk_counts_df$length)
#chack all 1m
colSums(bulk_tpms)
bulk_tpms = as.data.frame(bulk_tpms)
bulk_tpms$ID = row.names(bulk_tpms)



#------
#
# Merge pseudo & bulk
#
#------



sc_pseudo_bulk = data.frame(ID_0.1_counts = rowSums(save_counts_01[,c(ID_0.1_v2)]),
                                 ID_0.3_counts = rowSums(save_counts_03[,c(ID_0.3_v2)]),
                                 ID_1.0_counts = rowSums(save_counts_1[,c(ID_1.0_v2)]),
                                 ID_2.0_counts = rowSums(save_counts_2[,c(ID_2.0_v2)]),
                                 ID = save_counts_01$gene_id)
#remove dup rows
sc_pseudo_bulk = sc_pseudo_bulk[!duplicated(sc_pseudo_bulk),]

#add gene length
gff$length = gff$end - gff$start
sc_pseudo_bulk = merge(x=sc_pseudo_bulk, y=gff[,c(4,6)], by.x="ID",by.y="Name",all.x=T)
#add genes to rownames
row.names(sc_pseudo_bulk) = sc_pseudo_bulk$ID
#convert to TPM
sc_pseudo_bulk_tpm = tpm2(counts = sc_pseudo_bulk[,c(2:5)],len = sc_pseudo_bulk$length)
#remove NAs
sc_pseudo_bulk_tpm = sc_pseudo_bulk_tpm[complete.cases(sc_pseudo_bulk_tpm),]
#check
summary(colSums(sc_pseudo_bulk_tpm))
#convert to df
sc_pseudo_bulk_tpm = as.data.frame(sc_pseudo_bulk_tpm)
sc_pseudo_bulk_tpm$ID = row.names(sc_pseudo_bulk_tpm)





#---
# Comparing pseudo-bulk vs bulk
#---

merged_bulk_sc = merge(x=bulk_tpms, by.x="ID", y=sc_pseudo_bulk_tpm, by.y="ID", all.x=T)

#correlations
cor(log2(merged_bulk_sc$OD_01_I+1), log2(merged_bulk_sc$ID_0.1_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_01_II+1), log2(merged_bulk_sc$ID_0.1_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_01_III+1), log2(merged_bulk_sc$ID_0.1_counts+1), method="spearman")

cor(log2(merged_bulk_sc$OD_03_I+1), log2(merged_bulk_sc$ID_0.3_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_03_II+1), log2(merged_bulk_sc$ID_0.3_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_03_III+1), log2(merged_bulk_sc$ID_0.3_counts+1), method="spearman")

cor(log2(merged_bulk_sc$OD_1_I+1), log2(merged_bulk_sc$ID_1.0_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_1_II+1), log2(merged_bulk_sc$ID_1.0_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_1_III+1), log2(merged_bulk_sc$ID_1.0_counts+1), method="spearman")

cor(log2(merged_bulk_sc$OD_2_I+1), log2(merged_bulk_sc$ID_2.0_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_2_II+1), log2(merged_bulk_sc$ID_2.0_counts+1), method="spearman")
cor(log2(merged_bulk_sc$OD_2_III+1), log2(merged_bulk_sc$ID_2.0_counts+1), method="spearman")


#plots
cor_p1 = ggplot(merged_bulk_sc, aes(x=log2(ID_0.1_counts+0.5), y=log2(OD_01_I+0.5), alpha=0.5)) + 
  geom_point(size=3, shape =16, colour="#C2C2C2") + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD0.1")+ geom_abline(intercept =0 , slope = 1) + theme(legend.position = "none")

cor_p2 = ggplot(merged_bulk_sc, aes(x=log2(ID_0.3_counts+0.5), y=log2(OD_03_I+0.5), alpha=0.5)) + 
  geom_point(size=3, shape =16, colour="#FAE10F") + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD0.3")+ geom_abline(intercept =0 , slope = 1)+ theme(legend.position = "none")

cor_p3 = ggplot(merged_bulk_sc, aes(x=log2(ID_1.0_counts+0.5), y=log2(OD_1_I+0.5), alpha=0.5)) + 
  geom_point(size=3, shape =16, colour="#469BCD") + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD1.0")+ geom_abline(intercept =0 , slope = 1)+ theme(legend.position = "none")

cor_p4 = ggplot(merged_bulk_sc, aes(x=log2(ID_2.0_counts+0.5), y=log2(OD_2_I+0.5), alpha=0.5)) + 
  geom_point(size=3, shape =16, colour="#F56437") + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD2.0")+ geom_abline(intercept =0 , slope = 1)+ theme(legend.position = "none")


cairo_pdf(file="correlations_sc_bulk.pdf", width=21, height=5, bg = "white")
plot_grid(cor_p1, cor_p2, cor_p3, cor_p4, nrow = 1)
dev.off()


#Contour

cor_p1 = ggplot(merged_bulk_sc, aes(x=log2(ID_0.1_counts+0.5), y=log2(OD_01_I+0.5))) + 
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis", limits=c(0, 23), breaks=seq(0,20,by=5)) + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD0.1")+ geom_abline(intercept =0 , slope = 1) 

cor_p2 = ggplot(merged_bulk_sc, aes(x=log2(ID_0.3_counts+0.5), y=log2(OD_03_I+0.5))) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", limits=c(0, 23), breaks=seq(0,20,by=5)) + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD0.3")+ geom_abline(intercept =0 , slope = 1) 

cor_p3 = ggplot(merged_bulk_sc, aes(x=log2(ID_1.0_counts+0.5), y=log2(OD_1_I+0.5))) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", limits=c(0, 23), breaks=seq(0,20,by=5)) + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD1.0")+ geom_abline(intercept =0 , slope = 1) 

cor_p4 = ggplot(merged_bulk_sc, aes(x=log2(ID_2.0_counts+0.5), y=log2(OD_2_I+0.5))) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", limits=c(0, 23), breaks=seq(0,20,by=5)) + theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("OD2.0")+ geom_abline(intercept =0 , slope = 1) 


melted_contour_plots = data.frame("ID" = rep(merged_bulk_sc$ID, 4),
                                  "Xl" = c(merged_bulk_sc$OD_01_I,
                                           merged_bulk_sc$OD_03_I,
                                           merged_bulk_sc$OD_1_I,
                                           merged_bulk_sc$OD_2_I),
                                  "yl" = c(merged_bulk_sc$ID_0.1_counts,
                                           merged_bulk_sc$ID_0.3_counts,
                                           merged_bulk_sc$ID_1.0_counts,
                                           merged_bulk_sc$ID_2.0_counts),
                                  "OD" = rep(c("0.1","0.3","1","2"), each=length(merged_bulk_sc$ID))
)
melted_contour_plots$Xl = log2(melted_contour_plots$Xl+1)
melted_contour_plots$yl = log2(melted_contour_plots$yl+1)

cairo_pdf(file="correlations_sc_bulk_contour.pdf", width=12, height=12, bg = "white")
ggplot(melted_contour_plots, aes(x=Xl, y=yl)) + 
  geom_bin2d(bins = 70) + scale_fill_continuous(limits=c(0, 23), breaks=seq(0,23,by=5), type = "viridis")+ facet_wrap(~OD) +
  theme_ipsum() + ylim(c(0,20))+ xlim(c(0,20)) +
  xlab("log2(merged SC TPM)") + ylab("log2(bulk TPM)")+ ggtitle("dd")+ geom_abline(intercept =0 , slope = 1) 
dev.off()







#---
#
# Heatmap of SPI genes
#
#---

sp1_genes = c("sitA","sitB","sitC","sitD","avrA","sprB","hilC","STM2868","orgC","orgB","orgA","prgK","prgJ","prgI","prgH","hilD","hilA","iagB","sptP","sicP","iacP","sipA","sipD","sipC","sipB","sicA","spaS","spaR","spaQ","spaP","spaO","invJ","invI","invC","invB","invA","invE","invG","invF","invH","invR")
sp2_genes = c("ssrB","ssrA","ssaB","ssaC","ssaD","ssaE","sseA","sseB","sscA","sseC","sseD","sseE","sscB","sseF","sseG","ssaG","ssaH","ssaI","ssaJ","STM1410","ssaK","ssaL","ssaM","ssaV","ssaN","ssaO","ssaP","ssaQ","ssaR","ssaS","ssaT","ssaU")
sp4_genes = c("siiA","siiB","siiC","siiD","siiE","siiF")

flag_genes = c("flhD","flhC","fliE","fliB","fliF","fliG","fliH","fliI","fliJ","fliK",
               "fliL","fliM","fliN","fliO","fliP","fliQ","fliR","fliA","fliZ","fliY",
               "fliD","fliS","fliT","fliC","flhB","flhA","flhE","flgB","flgC","flgD",
               "flgE","flgF","flgG","flgH","flgU","flgJ","flgK","flgL","flgA","flgM",
               "flgN","flgI","fljB","fljA","cheZ","cheY","cheB","cheR","cheM","cheW",
               "cheA","motA","motB","hin","aer","tsr","flaG","tar")

sp_genes_df = counts_all_df_final_tpm2
rownames(sp_genes_df) = sp_genes_df$ID

#put in order of timepoints (for the heatmap)
sp_genes_df_v2 = sp_genes_df[,c("ID",ID_0.1_v2,ID_0.3_v2,ID_1.0_v2,ID_2.0_v2)]

sp_genes_df_v2[,!colnames(sp_genes_df_v2) %in%  cell_pheno_data_df$`Cell ID`]
#remove III_D1_S37_R1_001 to match pheno data
sp_genes_df_v3 = sp_genes_df_v2[,!colnames(sp_genes_df_v2) %in% c("III_D1_S37_R1_001")]

hvg_phenotypes[hvg_phenotypes$`Cell ID` %in% colnames(sp_genes_df_v3[c(2:4)]),]$Condition

#HM annotation
hm_anno_sp_genes = HeatmapAnnotation(Condition = hvg_phenotypes[hvg_phenotypes$`Cell ID` %in% colnames(sp_genes_df_v3[c(2:370)]),]$Condition,
                                     col = list(Condition = c("OD 0.1" = "#C2C2C2", "OD 0.3" = "#FAE10F","OD 1.0" = "#469BCD", "OD 2.0" = "#F56437"))
)

col_fun_spi_genes = circlize::colorRamp2(c(0, 1.5,3,4.5,6,7.5,9,10.5,14), RColorBrewer::brewer.pal(9, "PuBu"))
sp_genes_df_v4 = sp_genes_df_v3[rowSums(sp_genes_df_v3[,ID_2.0_v2]) > 1,]


#--
# SPI genes ESP
#--
ComplexHeatmap::Heatmap(log2(as.matrix( sp_genes_df_v4[sp_genes_df_v4$ID %in% c(sp1_genes,sp2_genes,sp4_genes),ID_2.0_v2])+1 ), 
                        cluster_rows = T, #top_annotation = hm_anno_sp_genes,
                        cluster_columns = T, show_column_names = FALSE, col = col_fun_spi_genes, name="log2(TPM)",
                        row_names_gp = gpar(fontsize = 7))



#--
# Flagella genes EEP and MEP
#--

#0.1
ComplexHeatmap::Heatmap(log2(as.matrix( sp_genes_df_v3[sp_genes_df_v3$ID %in% flag_genes,colnames(sp_genes_df_v3) %in% ID_0.1_v2])+1 ), 
                        cluster_rows = T, #top_annotation = hm_anno_sp_genes,
                        cluster_columns = T, show_column_names = FALSE, col = col_fun_spi_genes, name="log2(TPM)",
                        row_names_gp = gpar(fontsize = 7))

#0.3
ComplexHeatmap::Heatmap(log2(as.matrix( sp_genes_df_v3[sp_genes_df_v3$ID %in% flag_genes,colnames(sp_genes_df_v3) %in% ID_0.3_v2])+1 ), 
                        cluster_rows = T, #top_annotation = hm_anno_sp_genes,
                        cluster_columns = T, show_column_names = FALSE, col = col_fun_spi_genes, name="log2(TPM)",
                        row_names_gp = gpar(fontsize = 7))



#---
#
# Downsampling
#
#---



downsampling_df = counts_all_df_final[!counts_all_df_final$type %in% c("5UTR","3UTR"),]

#The numbers for the line plot
reads_to_sample = c(250,500,1000,2000,4000,8000,16000,32000,64000,128000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000)

#the numbers for the PCA plots
reads_to_sample = c(310,500,550, 575,600,1040,1050,2050,5000, 50000) 

#prepare empty df
downsampling_loop_df = data.frame("Reads" = reads_to_sample,
                                  "Average" = 0,
                                  "OD0.1" = 0,
                                  "OD0.3" = 0,
                                  "OD1.0" = 0,
                                  "OD2.0" = 0)

#Set seed
set.seed(261)

#loop through the values in reads_to_sample
for (g in 1:length(reads_to_sample)) {
  
  sub_sample_size = reads_to_sample[g]
  print(sub_sample_size)
  

  #--
  #OD 0.1
  #--
  #extract od 0.1 cells
  ds_loop_0.1 = downsampling_df[,colnames(downsampling_df) %in% c("gene_id",ID_0.1_v2)]
   
  all_sc_genes_df = data.frame("Gene" = downsampling_df$gene_id)
  for (i in 1:(dim(ds_loop_0.1)[2]-1)) {
    table_result = table(sample(ds_loop_0.1$gene_id, size=sub_sample_size, replace=T, prob=ds_loop_0.1[,i+1]))
    loop_results = data.frame("Genes" = names(table_result), "cell"=as.numeric(unname(table_result)))
    colnames(loop_results)[2] = paste0("Cell_0.1_",i)
    #merge df
    all_sc_genes_df = merge(x=all_sc_genes_df, by.x="Gene", y=loop_results, by.y="Genes", all.x=T, sort = F)
    
  }
  #replace all NA with zero
  all_sc_genes_0.1_df = all_sc_genes_df
  all_sc_genes_0.1_df[is.na(all_sc_genes_0.1_df)] <- 0
  
  #--
  #OD 0.3
  #--
  #extract od 0.1 cells
  ds_loop_0.3 = downsampling_df[,colnames(downsampling_df) %in% c("gene_id",ID_0.3_v2)]
  
  all_sc_genes_df = data.frame("Gene" = downsampling_df$gene_id)
  for (i in 1:(dim(ds_loop_0.3)[2]-1)) {
    table_result = table(sample(ds_loop_0.3$gene_id, size=sub_sample_size, replace=T, prob=ds_loop_0.3[,i+1]))
    loop_results = data.frame("Genes" = names(table_result), "cell"=as.numeric(unname(table_result)))
    colnames(loop_results)[2] = paste0("Cell_0.3_",i)
    #merge df
    all_sc_genes_df = merge(x=all_sc_genes_df, by.x="Gene", y=loop_results, by.y="Genes", all.x=T, sort = F)
    
  }
  #replace all NA with zero
  all_sc_genes_0.3_df = all_sc_genes_df
  all_sc_genes_0.3_df[is.na(all_sc_genes_0.3_df)] <- 0
  
  #--
  #OD 1.0
  #--
  #extract od 1.0 cells
  ds_loop_1.0 = downsampling_df[,colnames(downsampling_df) %in% c("gene_id",ID_1.0_v2)]
  
  all_sc_genes_df = data.frame("Gene" = downsampling_df$gene_id)
  for (i in 1:(dim(ds_loop_1.0)[2]-1)) {
    table_result = table(sample(ds_loop_1.0$gene_id, size=sub_sample_size, replace=T, prob=ds_loop_1.0[,i+1]))
    loop_results = data.frame("Genes" = names(table_result), "cell"=as.numeric(unname(table_result)))
    colnames(loop_results)[2] = paste0("Cell_1.0_",i)
    #merge df
    all_sc_genes_df = merge(x=all_sc_genes_df, by.x="Gene", y=loop_results, by.y="Genes", all.x=T, sort = F)
    
  }
  #replace all NA with zero
  all_sc_genes_1.0_df = all_sc_genes_df
  all_sc_genes_1.0_df[is.na(all_sc_genes_1.0_df)] <- 0
  
  
  #--
  #OD 2.0
  #--
  #extract od 1.0 cells
  ds_loop_2.0 = downsampling_df[,colnames(downsampling_df) %in% c("gene_id",ID_2.0_v2)]
  
  all_sc_genes_df = data.frame("Gene" = downsampling_df$gene_id)
  for (i in 1:(dim(ds_loop_2.0)[2]-1)) {
    table_result = table(sample(ds_loop_2.0$gene_id, size=sub_sample_size, replace=T, prob=ds_loop_2.0[,i+1]))
    loop_results = data.frame("Genes" = names(table_result), "cell"=as.numeric(unname(table_result)))
    colnames(loop_results)[2] = paste0("Cell_2.0_",i)
    #merge df
    all_sc_genes_df = merge(x=all_sc_genes_df, by.x="Gene", y=loop_results, by.y="Genes", all.x=T, sort = F)
    
  }
  #replace all NA with zero
  all_sc_genes_2.0_df = all_sc_genes_df
  all_sc_genes_2.0_df[is.na(all_sc_genes_2.0_df)] <- 0
  
  
  
  #Merge all conditions from above
  ds_merged_v1 = merge(x=all_sc_genes_0.1_df, by.x="Gene", y=all_sc_genes_0.3_df, by.y="Gene", all.x = T)
  ds_merged_v1 = merge(x=ds_merged_v1, by.x="Gene", y=all_sc_genes_1.0_df, by.y="Gene", all.x = T)
  ds_merged_v1 = merge(x=ds_merged_v1, by.x="Gene", y=all_sc_genes_2.0_df, by.y="Gene", all.x = T)
  
  #Remove rRNA
  rRNA_genes = downsampling_df[downsampling_df$type %in% c("rRNA"),"gene_id"]
  ds_merged_v2 = ds_merged_v1[!ds_merged_v1$Gene %in% rRNA_genes,]
  

  #look at summary
  dim(ds_merged_v2)
  mean(colSums(ds_merged_v2[,2:370] > 1))
  
  
  #create metadata
  pheno_ds_meta = data.frame("ID" = colnames(ds_merged_v2)[2:371])
  pheno_ds_meta$Condition = "OD 0.1"
  pheno_ds_meta[grep("Cell_0.3_",pheno_ds_meta$ID),]$Condition <- "OD 0.3"
  pheno_ds_meta[grep("Cell_1.0_",pheno_ds_meta$ID),]$Condition <- "OD 1.0"
  pheno_ds_meta[grep("Cell_2.0_",pheno_ds_meta$ID),]$Condition <- "OD 2.0"
  rownames(pheno_ds_meta) = pheno_ds_meta$ID
  
  
  #Perform PCA
  pca_ds = pca(log2(ds_merged_v2[,2:371]+1), metadata = pheno_ds_meta ,removeVar = 0.1)
  pca_cols2 =  c('OD0.1' = '#BCBCBB', 'OD0.3' = '#F3E226', 'OD1.0' = '#1298CE', 'OD2.0' = '#F26944')
  
  plot_save_name = paste0("Bottleneck_",sub_sample_size,".pdf")
  cairo_pdf(file=plot_save_name, width=6, height=5, bg = "white")
  print(biplot(pca_ds, colby = "Condition", showLoadings = F, lab = NULL, legendPosition = 'right', colkey = pca_cols))
  dev.off()
  
  mean(colSums(ds_merged_v2[,2:371] > 1))
  median(colSums(ds_merged_v2[,2:371] > 1))
  
  #for each OD
  downsampling_loop_df[g,2] <- mean(colSums(ds_merged_v2[,2:371] > 1))
  downsampling_loop_df[g,3] <- mean(colSums(ds_merged_v2[,pheno_ds_meta[pheno_ds_meta$Condition %in% "OD 0.1","ID"]] > 1))
  downsampling_loop_df[g,4] <- mean(colSums(ds_merged_v2[,pheno_ds_meta[pheno_ds_meta$Condition %in% "OD 0.3","ID"]] > 1))
  downsampling_loop_df[g,5] <- mean(colSums(ds_merged_v2[,pheno_ds_meta[pheno_ds_meta$Condition %in% "OD 1.0","ID"]] > 1))
  downsampling_loop_df[g,6] <- mean(colSums(ds_merged_v2[,pheno_ds_meta[pheno_ds_meta$Condition %in% "OD 2.0","ID"]] > 1))
  print(g)
}

downsampling_loop_df

#Should run the above loop twie, saving results into different variables 
#based upon the two variables reads_to_sample

#uncomment line relevant to input values above
#downsampling_loop_df__figure_PCA = downsampling_loop_df
#downsampling_loop_df__figure_line = downsampling_loop_df

downsampling_loop_df__figure_PCA
downsampling_loop_df__figure_line

cairo_pdf(file="Bottleneck_all.pdf", width=8, height=5, bg = "white")
ggplot(reshape2::melt(downsampling_loop_df__figure_line, id.vars="Reads"), aes(Reads,value, colour=variable))+ 
  geom_point() + xlab("Reads per cell") + ylab("Detected Genes") + geom_line()+labs(color='') + hrbrthemes::theme_ipsum()
dev.off()

cairo_pdf(file="Bottleneck_sub.pdf", width=8, height=5, bg = "white")
ggplot(reshape2::melt(downsampling_loop_df__figure_line[1:9,], id.vars="Reads"), aes(Reads,value, colour=variable))+ 
  geom_point() + xlab("Reads per cell") + ylab("Detected Genes") + geom_line() +labs(color='')+ hrbrthemes::theme_ipsum() 
dev.off()

