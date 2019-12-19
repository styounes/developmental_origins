library("Biobase")
library("GEOquery")
library("limma")
library("ggplot2")
library("hgu133plus2.db")
library("genefilter")
library("tidyverse")
library("topGO")
library("ReactomePA")
library("Rtsne")
library("amap")
library("pheatmap")
library("ggrepel")

# Working directory for MacOS workstation:
#setwd("/Users/styounes/Dropbox/Bioinformatics_Projects/Developmental_origins/")

# Working directory for server:
setwd("/Volumes/Data/Dropbox/Bioinformatics_Projects/Developmental_origins/")

# Working directory for Ubuntu workstation:
# setwd("~/Dropbox/Bioinformatics_Projects/Developmental_origins/")

# Function for linearizing the 6-depth list which contains the patient characteristics.
# It moves those characteristics to the end of the list each as a new element and removes the original element.
# It also removes the "relation" element from the control datasets.  This helps ensure the names are aligned in subsequent analysis.
remove_long_elements <- function(mylist) {
  index <- 1
  for(sample in mylist){
    for(element in sample){
      if(length(element) == 1){
        next
      } else {
        for(char in element){
          sample <- append(sample, char)
        }
      }
    }
    sample <- sample[names(sample) != "characteristics_ch1"]
    sample <- sample[names(sample) != "relation"]
    mylist[[index]] <- sample
    index <- index + 1
  }
  return(mylist)
}

# Generates the GSE dataset, complete with annotation data
GSE <- getGEO(filename = "GSE44971_family.soft")
gsmplatforms <- lapply(GSMList(GSE), function(x) {Meta(x)$platform_id})
gsmlist <- Filter(function(GSE) {Meta(GSE)$platform_id == "GPL570"}, GSMList(GSE))
probesets <- Table(GPLList(GSE)[[1]])$ID
data.matrix <- data.matrix <- do.call('cbind',lapply(gsmlist, function(x) {tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix, 2, function(x) {as.numeric(as.character(x))})
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)

# This correctly formats the phenotype data for subsequent insertion into the ExpressionSet object
pheno <- list()
index <- 1
for(i in GSE@gsms) {
  pheno[[index]] <- i@header
  index <- index + 1
}
names(pheno) <- names(gsmlist)
pheno <- remove_long_elements(pheno)
new_names <- as.list(names(pheno$GSM1094814))
new_names[32:37] <- c("tissue", "brain_location", "age_at_diagnosis", "gender", "dna_methylation_subgroup", "mapk_hit")
index <- 1
for(thing in pheno){
  names(pheno[[index]]) <- new_names
  index <- index + 1
}
pheno_matrix <- matrix(unlist(pheno), nrow = length(unlist(pheno[1])))
pheno_matrix <- t(pheno_matrix)
colnames(pheno_matrix) <- new_names
rownames(pheno_matrix) <- pheno_matrix[, "geo_accession"]
pheno_matrix[, "tissue"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "tissue"]))
pheno_matrix[, "brain_location"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "brain_location"]))
pheno_matrix[, "age_at_diagnosis"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "age_at_diagnosis"]))
pheno_matrix[, "gender"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "gender"]))
pheno_matrix[, "dna_methylation_subgroup"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "dna_methylation_subgroup"]))
pheno_matrix[, "mapk_hit"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "mapk_hit"]))
pheno_matrix[, "mapk_hit"] <- trimws(gsub("^.*?:", "", pheno_matrix[, "mapk_hit"]))
pheno_df <- as.data.frame(pheno_matrix)
AgeBracket <- ifelse(as.integer(pheno_df$age_at_diagnosis) >= 13, "Less than 13", "Greater than 13")
index <- 1
Location <- list()
for (location in pheno_df$brain_location) {
  if(location == "Cerebellum") {
    Location[[index]] <- "Infratentorial"
  } else if (location == "Fourth Ventricle"){
    Location[[index]] <- "Infratentorial"
  } else if (location == "Brainstem") {
    Location[[index]] <- "Infratentorial"
  } else {
    Location[[index]] <- "Supratentorial"
  }
  index <- index + 1
}
Mutation <- list()
index <- 1
for(mutation in pheno_df$mapk_hit) {
  if(mutation == "KIAA1549:BRAF fusion") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "KIAA1549:BRAF fusion_15:9") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "KIAA1549:BRAF fusion_16:9") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "KIAA1549:BRAF fusion_16:11") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "KIAA1549:BRAF fusion_breakpoint not determined") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "FAM131B:BRAF fusion") {
    Mutation[[index]] <- "BRAF Fusion"
  } else if (mutation == "NF1") {
    Mutation[[index]] <- "NF1"
  } else if (mutation == "NF1; KIAA1549:BRAF fusion_15:9; BRAF V600E") {
    Mutation[[index]] <- "NF1"
  } else if (mutation == "NF1; KIAA1549:BRAF fusion_breakpoint not determined") {
    Mutation[[index]] <- "NF1"
  } else if (mutation == "BRAF V600E") {
    Mutation [[index]] <- "BRAF Point Mutation"
  } else if (mutation == "BRAF ins598T") {
    Mutation [[index]] <- "BRAF Point Mutation"
  } else {
    Mutation[[index]] <- "Not identified"
  }
  index <- index + 1
}
pheno_df <- cbind(pheno_df, Location = as.character(Location))
pheno_df <- cbind(pheno_df, Mutation = as.character(Mutation))
pheno_df <- cbind(pheno_df, Age_Bracket = AgeBracket)
pheno_df <- as(pheno_df,"AnnotatedDataFrame")

# Finally, builds the expression set and subsets for easy access
eset <- new('ExpressionSet', exprs = data.matrix, phenoData = pheno_df)

# Builds a clustering heatmap of the data
distance_matrix <- as.matrix(Dist(t(exprs(eset)), method = "pearson"))
rownames(distance_matrix) <- row.names(pData(eset))
colnames(distance_matrix) <- NULL
diag(distance_matrix) <- NA
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
annotation_for_heatmap <- data.frame(Mutation = phenoData(eset)$Mutation, Location = phenoData(eset)$Location)
row.names(annotation_for_heatmap) <- row.names(phenoData(eset))
pheatmap::pheatmap(distance_matrix, col = (hmcol), annotation_row = annotation_for_heatmap, 
                   legend = TRUE, treeheight_row = 0, 
                   legend_breaks = c(min(distance_matrix, na.rm = TRUE), max(distance_matrix, na.rm = TRUE)), 
                   legend_labels = (c("Small Distance", "Large Distance")), main = "Clustering Heatmap of All Samples")

tsne_dist_matrix <- as.dist(Dist(t(exprs(eset)), method = "pearson"))
set.seed(24)
tsne_out <- Rtsne(tsne_dist_matrix, perplexity = 19, is_distance = TRUE, theta = 0, max_iter = 2500, pca = FALSE, num_threads = 2)
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = phenoData(eset)$Location)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

# Builds a PCA plot
exp_raw <- exprs(eset)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Location = pData(eset)$Location,
                     Methylation = pData(eset)$dna_methylation_subgroup)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Methylation, colour = Location)) +
  ggtitle("PCA plot of all samples (including normal cerebellum)") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio)

#Plot just the tumors
tumors <- eset[,phenoData(eset)$source_name_ch1 != "Normal cerebellum RNA (commercial source)"]
distance_matrix <- as.matrix(Dist(t(exprs(tumors)), method = "pearson"))
rownames(distance_matrix) <- row.names(pData(tumors))
colnames(distance_matrix) <- NULL
diag(distance_matrix) <- NA
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
annotation_for_heatmap <- data.frame( Location = phenoData(tumors)$Location, Mutation = phenoData(tumors)$Mutation)
row.names(annotation_for_heatmap) <- row.names(phenoData(tumors))
pheatmap::pheatmap(distance_matrix, col = (hmcol), annotation_row = annotation_for_heatmap, 
                   legend = TRUE, treeheight_row = 0, 
                   legend_breaks = c(min(distance_matrix, na.rm = TRUE), max(distance_matrix, na.rm = TRUE)), 
                   legend_labels = (c("Small Distance", "Large Distance")), main = "Clustering Heatmap of All Tumors")

tsne_dist_matrix <- as.dist(Dist(t(exprs(tumors)), method = "pearson"))
set.seed(24)
tsne_out <- Rtsne(tsne_dist_matrix, perplexity = 16, is_distance = TRUE, theta = 0, max_iter = 2500, pca = FALSE, num_threads = 2)
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = phenoData(tumors)$Location)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

exp_raw <- exprs(tumors)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Location = pData(tumors)$Location,
                     Methylation = pData(tumors)$dna_methylation_subgroup)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Methylation, colour = Location)) +
  ggtitle("PCA Plot of All Tumors") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio)

# Creates annotations and filters multiple probes
anno_eset <- AnnotationDbi::select(hgu133plus2.db, keys = (featureNames(eset)), columns = c("SYMBOL", "GENENAME"), 
                                   keytype = "PROBEID")
anno_eset <- subset(anno_eset, !is.na(SYMBOL))
anno_group <- group_by(anno_eset, PROBEID)
anno_summarized <- dplyr::summarize(anno_group, no_of_matches = n_distinct(SYMBOL))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(eset) %in% anno_filtered$PROBEID)
eset <- subset(eset, !ids_to_exclude)
fData(eset)$PROBEID <- rownames(fData(eset))
fData(eset) <- left_join(fData(eset), anno_eset)
rownames(fData(eset)) <- fData(eset)$PROBEID

# Extracts solely BRAF fusion positive tumors and graphs the result
Fusion_Positive <- (phenoData(eset)$Mutation %in% "BRAF Fusion")
BRAF_Fusion_eset <- eset[,Fusion_Positive]
TrueCerebellum <- BRAF_Fusion_eset$brain_location != "Brainstem"
BRAF_Fusion_eset <- BRAF_Fusion_eset[, TrueCerebellum]
TrueCerebellum <- BRAF_Fusion_eset$Location != "Fourth Ventricle"
BRAF_Fusion_eset <- BRAF_Fusion_eset[, TrueCerebellum]

distance_matrix <- as.matrix(amap::Dist(t(exprs(BRAF_Fusion_eset)), method = "pearson"))
rownames(distance_matrix) <- row.names(pData(BRAF_Fusion_eset))
colnames(distance_matrix) <- row.names(pData(BRAF_Fusion_eset))
diag(distance_matrix) <- NA
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
#christmas <- c("#BB2528", "#146B3A")
#hmcol <- rev(colorRampPalette(christmas)(255))
annotation_for_heatmap <- data.frame(Location = phenoData(BRAF_Fusion_eset)$Location)
row.names(annotation_for_heatmap) <- row.names(phenoData(BRAF_Fusion_eset))

my_colors <- list(Location = c(Infratentorial = "#00BFC4", Supratentorial = "#F8766D"))
#my_colors <- list(Location = c(Infratentorial = "#146B3A", Supratentorial = "#BB2528"))
pheatmap::pheatmap(distance_matrix, col = (hmcol), annotation_col = annotation_for_heatmap, 
                   annotation_colors = my_colors, annotation_legend = TRUE,
                   legend = TRUE, treeheight_row = 0, show_colnames = FALSE, show_rownames = FALSE,
                   legend_breaks = c(min(distance_matrix, na.rm = TRUE), max(distance_matrix, na.rm = TRUE)), 
                   legend_labels = (c("Small Distance", "Large Distance")),
                   main = "Clustering Heatmap of BRAF Fusion-Positive Tumors", fontsize = 18, height = 50, width = 50)

#annotation_for_heatmap <- HeatmapAnnotation(Location = phenoData(BRAF_Fusion_eset)$Location, col = my_colors)
#Heatmap(t(exprs(BRAF_Fusion_eset)), name = "Clustering Heatmap of BRAF Fusion-Positive Tumors",
#        clustering_distance_columns = "pearson", col = hmcol, top_annotation = annotation_for_heatmap)

set.seed(42)
tsne_out <- Rtsne(t(exprs(BRAF_Fusion_eset)), perplexity = 12, theta = 0, max_iter = 2500, num_threads = 2)
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], Location = phenoData(BRAF_Fusion_eset)$Location)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=Location), size = 3.5) + 
  scale_color_manual(values = c(Infratentorial = "#00BFC4", Supratentorial = "#F8766D")) + 
  ggtitle("t-SNE Dimensionality Reduction of BRAF Fusion-Positive Tumors") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
  labs(color = "Location")

# PCA Plot
exp_raw <- exprs(BRAF_Fusion_eset)
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Location = pData(BRAF_Fusion_eset)$Location)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Location)) +
  ggtitle("PCA plot of BRAF Fusion Positive Tumors") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio)

# Make the design matrix
Location <- pData(BRAF_Fusion_eset)$Location
design_matrix <- model.matrix(~ 0 + Location)
colnames(design_matrix) <- c("Infratentorial", "Supratentorial")
rownames(design_matrix) <- pData(BRAF_Fusion_eset)$geo_accession

# Run the contrasts
contrast_matrix <- makeContrasts(Infratentorial-Supratentorial, levels = design_matrix)
InfraVsSupra <- eBayes(contrasts.fit(lmFit(BRAF_Fusion_eset, design = design_matrix), contrast_matrix))
InfraVsSupra_Table <- topTable(InfraVsSupra, number = Inf)

# Make a volcano plot of the differentially expressed genes, naming the top 100
threshold <- (InfraVsSupra_Table$P.Value < 0.0001 & abs(InfraVsSupra_Table$logFC) > 0.5)
InfraVsSupra_Table$threshold <- threshold

ggplot(InfraVsSupra_Table) + geom_point(aes(x=logFC, y=-log10(P.Value), color = threshold), size = 3) + 
  ggtitle("Infratentorial vs. Supratentorial Tumors") + xlab("Log2 Fold Change") + ylab("-log10 P-Value") +
  scale_color_manual(values = c("#A9A9A9", "#009051")) +
  geom_vline(xintercept = 0.5, col = "black", linetype = "dotted", size = 1) +
  geom_vline(xintercept = -0.5, col = "black", linetype = "dotted", size = 1) +
  geom_hline(yintercept = 4, col = "black", linetype = "dotted", size = 1) +
  geom_label_repel(aes(x=logFC, y=-log10(P.Value), 
                       label = ifelse(threshold == TRUE & SYMBOL %in% 
                                        c("PAX3", "IRX5", "ASCL1", "IRX2", "MEIS1", "KLF15", "PBX3", "MSX2"), SYMBOL, "")),
                       point.padding = 0.2, size = 5) +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5), axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15), legend.position = "none")

# Gene Ontology Analysis
DE_IvsS <- subset(InfraVsSupra_Table, adj.P.Val < 0.05)$PROBEID
back_genes <- genefinder(BRAF_Fusion_eset, as.character(DE_IvsS), method = "euclidean", scale = "none")
back_genes <- sapply(back_genes, function(x)x$indices)
back_genes <- featureNames(BRAF_Fusion_eset)[back_genes]
back_genes <- setdiff(back_genes, DE_IvsS)
gene_IDs <- rownames(InfraVsSupra_Table)
in_universe <- gene_IDs %in% c(DE_IvsS, back_genes)
in_selection <- gene_IDs %in% DE_IvsS
all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 
top_GO_data <- new("topGOdata", ontology = "BP", allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db")
result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")
res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)
genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hgu133plus2.db", geneCutOff = 1000)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), collapse = "")})
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 1, useInfo = 'def')

# Reactome Pathway Analysis
#entrez_ids <- mapIds(hgu133plus2.db, keys = rownames(InfraVsSupra_Table), keytype = "PROBEID", column = "ENTREZID")
#reactome_enrich <- enrichPathway(gene = entrez_ids[DE_IvsS], universe = entrez_ids[c(DE_IvsS, back_genes)], organism = "human",
#                                 pvalueCutoff = 0.05, qvalueCutoff = 0.9, readable = TRUE)
#reactome_enrich@result$Description <- paste0(str_sub(reactome_enrich@result$Description, 1, 20), "...")
#barplot(reactome_enrich)

# Further curation of gene expression data (i.e. only those which are significant and upregulated)
InfraVsSupra_Table_SigUp <- subset (InfraVsSupra_Table, InfraVsSupra_Table$adj.P.Val < 0.05 & InfraVsSupra_Table$logFC > 0)
InfraVsSupra_Table_SigDown <- subset (InfraVsSupra_Table, InfraVsSupra_Table$adj.P.Val < 0.05 & InfraVsSupra_Table$logFC < 0)
