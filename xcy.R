getwd()

#load packages ####
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)
library(SingleR)
library(celldex)
library(dittoSeq)
library(clusterProfiler)
library(org.Mm.eg.db)   
library(enrichplot)
library(msigdbr)
library(ggupset)
library(EnhancedVolcano)
library(scDblFinder)
library(GOSemSim)
library(tidyr)
library(DESeq2)
library(CellChat)
library(ggalluvial)
library(ComplexHeatmap)
library(patchwork)
library(presto)
library(NMF)
library(SeuratWrappers)
library(monocle3)
library(BiocManager)

options(future.globals.maxSize = 8000 * 1024^2) 

#load files ####
sc.data = Read10X(data.dir = "Mouse_Chondrocytes_5K/filtered_feature_bc_matrix")
sc.data = CreateSeuratObject(sc.data)

anno = read.table("Mouse_Chondrocytes_5K/annotations.csv", header = TRUE, row.names = 1, sep = "\t")
head(anno)
colnames(anno)

#add columns ####
meta.use = anno[, c("Sample", "Sample_type", "rownames.idh")]
colnames(meta.use) = c("Sample", "Sample_type", "idh")

sc.data = AddMetaData(sc.data, metadata = meta.use)

head(sc.data@meta.data[, c("Sample", "Sample_type", "idh")])
table(sc.data$Sample_type)
table(sc.data$idh)

#make violin plot
VlnPlot(sc.data, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0) + NoLegend()

#make feature vs feature scatterplot 
png("countvsfeature.png", height = 1500, width = 1500, res = 300)
FeatureScatter(sc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

#Add the percent_mito to the final plot ####
sc.data = PercentageFeatureSet(sc.data, pattern = "^mt-", col.name = "percent_mito")
VlnPlot(sc.data, features = c("nCount_RNA", "nFeature_RNA", "percent_mito"), pt.size = 0) + NoLegend()
FeatureScatter(sc.data, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(sc.data, feature1 = "nFeature_RNA", feature2 = "percent_mito")

#delete dead cells and subset data ####
sc.data = subset(sc.data,subset = nCount_RNA > 1000 & nFeature_RNA > 1000 & percent_mito<10 )
png("violin_count_feature_mito.png", height = 800, width = 2000, res = 300)
VlnPlot(sc.data, features = c("nCount_RNA","nFeature_RNA","percent_mito"), pt.size = 0) + NoLegend()
dev.off()

FeatureScatter(sc.data, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(sc.data, feature1 = "nFeature_RNA", feature2 = "percent_mito")

#SC Transform ####
sc.data = SCTransform(sc.data, vars.to.regress = "percent_mito", verbose = TRUE)

#RunPCA before integrated
sc.data = RunPCA(object = sc.data, assay = "SCT", reduction.name = "pca_unintegrated")

#VizDimLoadings
png("PCA_loadings_unintegrated.png", height = 1000, width = 1000)
VizDimLoadings(sc.data, dims = 1:2, reduction = "pca_unintegrated")
dev.off()

#PCA plot
png("PCA_plot_unintegrated.png", height = 800, width = 800)
DimPlot(sc.data, reduction = "pca_unintegrated", dims = c(1,2))
dev.off()

#heatmap
png("PCA_heatmap_unintegrated.png", height = 1000, width = 1000)
DimHeatmap(sc.data, dims = 1:5, cells = 5000, balanced = TRUE, reduction = "pca_unintegrated")
dev.off()

#Elbow ####
png("Elbow_plot_unintegrated.png", height = 800, width = 800)
ElbowPlot(sc.data, ndims = 50, reduction = "pca_unintegrated")
dev.off()

dims_to_use = 30

#UMAP before integrated
sc.data = FindNeighbors(sc.data, reduction = "pca_unintegrated", dims = 1:dims_to_use)
sc.data = FindClusters(sc.data, resolution = 0.5)
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, 
                   reduction = "pca_unintegrated", 
                   reduction.name = "umap_unintegrated")

png("UMAP_unintegrated.png", height = 1000, width = 2000)
p1 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "orig.ident")
p2 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "seurat_clusters")
p1 + p2
dev.off()

#Harmony ####
library(harmony)
sc.data = RunHarmony(sc.data, 
                      assay.use = "SCT", 
                      group.by.vars = c("orig.ident"), 
                      theta = c(3), 
                      reduction = "pca_unintegrated", 
                      reduction.save = "pca_integrated")

#usec PCA to draw UMAP
sc.data = FindNeighbors(sc.data, reduction = "pca_integrated", dims = 1:dims_to_use)
sc.data = FindClusters(sc.data, resolution = 0.5)
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, 
                   reduction = "pca_integrated", 
                   reduction.name = "umap_integrated")

png("UMAP_integrated.png", height = 1000, width = 2000)
p1 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "orig.ident")
p2 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "seurat_clusters")
p1 + p2
dev.off()

png("umap.png", height = 800, width = 800)
DimPlot(sc.data, reduction = "umap_integrated")
dev.off()

#### find doublets cells ####
dbl.dens = computeDoubletDensity(sc.data@assays$SCT@data)
sc.data$DoubletScore = dbl.dens

png("DoubletScore.png", height = 700, width = 700)
FeaturePlot(sc.data, "DoubletScore")
dev.off()

dbl.calls = doubletThresholding(data.frame(score = dbl.dens), 
                                 method = "griffiths", 
                                 returnType = "call")
summary(dbl.calls)
sc.data$DoubletCalls = dbl.calls

# UMAP of Doublets
png("Doublet_calls_UMAP.png", height = 800, width = 800)
DimPlot(sc.data, group.by = "DoubletCalls", reduction = "umap_integrated")
dev.off()

# remove doublets
sc.data = subset(sc.data, subset = DoubletCalls == "singlet")

# remove most non-Chondrocyte cells
sc.data = subset(sc.data, subset = Col2a1 > 2)

# sct
sc.data = SCTransform(sc.data, vars.to.regress = "percent_mito", verbose = TRUE)

# PCA
sc.data = RunPCA(object = sc.data, assay = "SCT", reduction.name = "pca_unintegrated")

# rethink the dim to use
# Elbow
png("Elbow_plot_unintegrated2.png", height = 800, width = 800)
ElbowPlot(sc.data, ndims = 50, reduction = "pca_unintegrated")
dev.off()
dims_to_use = 25

#FindCluster and runumap after remove dloublets ####
sc.data = FindNeighbors(sc.data, reduction = "pca_unintegrated", dims = 1:dims_to_use)
sc.data = FindClusters(sc.data, resolution = 0.4)
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, 
                   reduction = "pca_unintegrated", 
                   reduction.name = "umap_unintegrated")

png("UMAP_unintegrated.png", height = 1000, width = 2000)
p1 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "orig.ident")
p2 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "seurat_clusters")
p1 + p2
dev.off()

#redo Harmony
sc.data = RunHarmony(sc.data, 
                      assay.use = "SCT", 
                      group.by.vars = c("orig.ident"), 
                      theta = c(3), 
                      reduction = "pca_unintegrated", 
                      reduction.save = "pca_integrated")

#usec PCA to draw UMAP
sc.data = FindNeighbors(sc.data, reduction = "pca_integrated", dims = 1:dims_to_use)
sc.data = FindClusters(sc.data, resolution = 0.4)
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, 
                   min.dist = 0.2, 
                   spread = 0.5,
                   reduction = "pca_integrated", 
                   reduction.name = "umap_integrated")

png("UMAP_integrated_afterremove.png", height = 700, width = 1400)
p1 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "orig.ident")
p2 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "seurat_clusters")
p1 + p2
dev.off()

saveRDS(sc.data, file = "sc.data.rds")

#find markers gene ####
sc.data = PrepSCTFindMarkers(sc.data)
sc.data.markers = FindAllMarkers(sc.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top15_markers = sc.data.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup()

#save heatmap to a PNG
png("Top15Markers_Heatmap.png", height = 3000, width = 1500, res = 300)
DoHeatmap(sc.data, features = top15_markers$gene) + NoLegend()
dev.off()

#### redo cell type identification ####

# Make sure the active meta.data slot is the seurat_clusters
sc.data = SetIdent(sc.data, value = sc.data$seurat_clusters)

# get the markers
sc.data.markers = FindAllMarkers(sc.data, only.pos = F)

top10= sc.data.markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1) %>% 
  slice_head(n = 10) %>% 
  ungroup()
print(top10)

# manual cell typing
new.cluster.ids = c("resting_notsure","resting","Neuron",
                    "articular","prehypertrophic",
                    "proliferating",
                    "enchondromas/mutant-specific","Stress-related cells")

# creat cell.types
names(new.cluster.ids) = levels(sc.data)
sc.data = RenameIdents(sc.data, new.cluster.ids)
sc.data = AddMetaData(sc.data, sc.data@active.ident, col.name = "cell.types")

# UMAP by celltypes
png("UMAP_cell_types.png", height = 1000, width = 1500)
p1 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "seurat_clusters", label = TRUE)
p2 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "cell.types", label = TRUE)
p1 + p2
dev.off()

# Cell frequencies
png("Cell_frequencies.png", height = 800, width = 1200)
dittoBarPlot(object = sc.data, var = "cell.types", group.by = "orig.ident")
dev.off()


#### pseduotime ####

colours_clusters = c("#FA0202", "#FAD902", "#70CF5B", "#080707","#0FD4C0",
                      "#D40FC4", "#0A3AAB", "#A8A9B3","#A67A0F", "#E6A1B5")
#subset
selected_idents = c("resting", "resting_notsure","articular","prehypertrophic", 
                     "proliferating", "enchondromas/mutant-specific")
table(Idents(sc.data))
cells.for.traj = subset(sc.data, idents = selected_idents)
cds = as.cell_data_set(cells.for.traj)

# Seurat UMAP
p1 = DimPlot(cells.for.traj, reduction = "umap_integrated", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE, cols = colours_clusters) +
  ggtitle("Seurat clusters") + theme(plot.title = element_text(hjust = 0.5))

# UMAP
p2 = DimPlot(cells.for.traj, reduction = "umap_integrated", group.by = "cell.types",
              label = TRUE, repel = TRUE, cols = colours_clusters) +
  ggtitle("Cell types") + theme(plot.title = element_text(hjust = 0.5))

gene.meta.table = colData(cds)
head(gene.meta.table)

# gene
fData(cds)$gene_short_name = rownames(fData(cds))
head(fData(cds))

# assign partitions (“create seurat_clusters” column in Seurat)
recreate.partition = c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) = cds@colData@rownames
recreate.partition = as.factor(recreate.partition)

# save
cds@clusters$UMAP$partitions = recreate.partition

# Assign the cluster info
list_cluster = cells.for.traj@active.ident
cds@clusters$UMAP$clusters = list_cluster

#  Assign UMAP coordinate - cell embeddings 
cds@int_colData@listData$reducedDims$UMAP = cells.for.traj@reductions$umap_integrated@cell.embeddings

#  minimum spanning tree
cds = learn_graph(cds, use_partition = FALSE)

# draw UMAP tree
p3 = plot_cells(cds,
                 color_cells_by = 'cell.types',
                 label_cell_groups = F,
                 label_groups_by_cluster = F,
                 label_branch_points = F,
                 label_roots = F,
                 label_leaves = F,
                 cell_size = 0.5)
p3
ggsave("UMAP_with_trajectory.png", width = 8, height = 6)


#choose root
cds = order_cells(cds, reduction_method = "UMAP")

#add root into umap
p1 = plot_cells(cds,
                 color_cells_by = 'pseudotime',
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = TRUE,
                 label_leaves = FALSE)
p1
ggsave("UMAP_with_root_trajectory.png", width = 8, height = 6)

#check Pseudotime informations
Pseudotime = pseudotime(cds)
cds$monocle3_pseudotime = pseudotime(cds)
data.pseudo = as.data.frame(colData(cds))

# boxplot：pseudotime
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(cell.types, monocle3_pseudotime, median), fill = cell.types)) +
  geom_boxplot()+
  theme_classic()+
  xlab("Pseudotime")
ggsave("Pseudotime_boxplot.png", width = 8, height = 6)

# test pseudotime related gene
deg_cells = graph_test(cds,
                         neighbor_graph = 'principal_graph',
                         cores = 4)

# add pseudotime back to Seurat
cells.for.traj$pseudotime = pseudotime(cds)
Idents(cells.for.traj) = cells.for.traj$cell.types
FeaturePlot(cells.for.traj, features = "pseudotime", label = T) 

# significant gene
deg_cells = deg_cells[order(deg_cells$p_value), ]
head(deg_cells)

# FeaturePlot
FeaturePlot(cells.for.traj, features = "pseudotime", label = TRUE, repel = TRUE) +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Pseudotime on Seurat UMAP")
ggsave("Seurat_FeaturePlot_pseudotime.png", width = 8, height = 6)


####expression analysis in Seurat ####
genes_to_test = rownames(GetAssayData(cells.for.traj, layer = "data"))[1:100]

# Spearman for every gene
pseudotime_cor = sapply(genes_to_test, function(gene) {
  expr_data = GetAssayData(cells.for.traj, layer = "data")[gene, ]
  cor(cells.for.traj$pseudotime, as.numeric(expr_data), method = "spearman")
})

# 10 genes most relate to time
pseudotime_cor = sort(pseudotime_cor, decreasing = TRUE)
head(pseudotime_cor, 10)
#outcome：
# Mybl1  Eya1  Cpa6  Prim2   Mcm3   Bend6  Terf1  Dst   Prex2   Fam135a 


#use FindMarkers compare earlyer and later cells
cells.for.traj$pseudotime_group = ifelse(cells.for.traj$pseudotime < median(cells.for.traj$pseudotime, na.rm = TRUE),
                                   "Early", "Late")
Idents(cells.for.traj) = cells.for.traj$pseudotime_group
early_vs_late_markers = FindMarkers(cells.for.traj, ident.1 = "Late", ident.2 = "Early")
head(early_vs_late_markers)
#outcome:
#           p_val    avg_log2FC  pct.1  pct.2     p_val_adj
#Col3a1 9.376750e-241  -2.145966 0.509 0.906 1.840750e-236
#Ebf1   9.441882e-239   2.639399 0.851 0.423 1.853536e-234
#Runx2  7.167512e-212   2.574415 0.845 0.453 1.407054e-207
#Meg3   1.767231e-211  -1.014811 0.987 0.993 3.469251e-207
#Dlk1   1.188568e-186  -1.431994 0.984 0.985 2.333278e-182
#Nfib   2.256310e-162  -1.408235 0.571 0.905 4.429362e-158


#### Cell to Cell Interactions ####
DefaultAssay(sc.data) = "RNA"
sc.data = NormalizeData(sc.data)

sc.data$samples = sc.data$orig.ident 
cellchat_data = createCellChat(object = sc.data,
                               group.by = "cell.types",
                               assay = "RNA")
CellChatDB = CellChatDB.mouse 
CellChatDB = subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat_data@DB = CellChatDB

cellchat_data = subsetData(cellchat_data)

cellchat_data = identifyOverExpressedGenes(cellchat_data)
cellchat_data = identifyOverExpressedInteractions(cellchat_data)
cellchat_data = computeCommunProb(cellchat_data, type = "triMean")

cellchat_data = filterCommunication(cellchat_data, min.cells = 10)
cellchat_data = computeCommunProbPathway(cellchat_data)
interactions = subsetCommunication(cellchat_data)

head(interactions)

# network
cellchat_data = aggregateNet(cellchat_data)
cellchat_data = netAnalysis_computeCentrality(cellchat_data, slot.name = "netP")

#### circos of all cell types ####
groupSize = as.numeric(table(cellchat_data@idents))

png("netVisual_circle.png", height = 1000, width = 1000)
netVisual_circle(cellchat_data@net$count,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of interactions")
dev.off()

netVisual_circle(cellchat_data@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction strength")


####  each cells interactions individually ####
mat = cellchat_data@net$weight
pdf("cellchat_circle_plots.pdf", width = 12, height = 9)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) 
{
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
} 
dev.off()
par(mfrow=c(1,1))

#hearmap
png("netVisual_heatmap_count.png", height = 1200, width = 1200)
netVisual_heatmap(cellchat_data, measure = "count", 
                  font.size = 20)
dev.off()
png("netVisual_heatmap_weight.png", height = 1000, width = 1000)
netVisual_heatmap(cellchat_data, measure = "weight", 
                  font.size = 20)
dev.off()

# split by sample group
p6 = DimPlot(sc.data, 
              reduction = "umap_integrated",
              group.by = "cell.types", 
              split.by = "Sample_type", 
              ncol = 2)
print(p6)
ggsave("ControlvsMutant.png", width = 12, height = 6, dpi = 300)

# subset by sample group
data.con = subset(sc.data, subset = Sample_type == "Control")
data.mut = subset(sc.data, subset = Sample_type == "Mutant")

#### passway analysis ####
# Marker gene of enchondromas/mutant-specific cells
mut_specific_markers = FindMarkers(sc.data, ident.1 = "enchondromas/mutant-specific")
sig_genes = mut_specific_markers %>% 
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  rownames()
gene_ids = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)

# GO
go_result = enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")

# dotplot for enchondromas cells
p7 = dotplot(go_result, showCategory = 15)
print(p7)
ggsave("GO_Enrichment_Mutant_Specific.png", plot = p7, width = 8, height = 6, dpi = 300)
