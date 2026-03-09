

# Define input directory
datadir <- "0data/private/Xenium_B6CAST/"

# Command args
args <- commandArgs(trailingOnly = T)
sample <- args[1] ## Sample to analyse
# sample <- "M12"

print(paste("Sample:", sample))

# Load packages
suppressPackageStartupMessages(library(dplyr, quietly = T))
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(Seurat, quietly = T))
suppressPackageStartupMessages(library(future, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))
suppressPackageStartupMessages(library(patchwork, quietly = T))

print("Seurat version:")
print(packageVersion(pkg="Seurat"))

# Do this
options(future.globals.maxSize = 8000 * 1024^2)

# Other options
remove_cells_low_expr <- T
low_expr_threshold <- 0

### ======  IMPORT DATA ====== ###

# Import sample_info.csv
sample_info <- data.table::fread("0data/private/Xenium_B6CAST/sample_info.csv")

# Import codewords.csv
codewords <- data.table::fread("0data/private/Xenium_B6CAST/codewords.csv")

# Define data dir for each sample
indirs <- list.dirs(datadir,full.names = T, recursive = F)
names(indirs) <- sub("__.*","",sub(".*output-XETG00289__0052...__","",indirs))

# Check samples and dir names
if(!all(names(indirs) %in% sample_info$Name)) { stop("Names of input directories are not equal to sample names") }

# Define available samples
samples <- sample_info$Name

# Check that specified sample is in the available samples
if(!sample %in% samples) { stop("Specified sample is not one of the available samples") }

# Define new sample name
sample_name <- sample_info %>% dplyr::filter(Name==sample) %>% dplyr::pull(Name2)
sample_name_full <- paste0(sample_name," - ",sample) 

print(paste("Full sample name:", sample_name_full))

# Define input dir to inpurt
indir <- indirs[sample]

# Load the Xenium data
xenium.obj <- Seurat::LoadXenium(indir, fov = "fov", cell.centroids = T, segmentations = "cell")
# xenium.obj <- Seurat::LoadXenium(indir, fov = "fov")

# Get total number of cells
total_num_cells <- ncol(xenium.obj)

# Remove cells with 0 counts and if applicable, low expressed, and get the numbers
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
num_cells_no_zero_counts <- ncol(xenium.obj)
if (remove_cells_low_expr){
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > low_expr_threshold)
  num_cells_low_expr <- ncol(xenium.obj)
  num_cells_no_zero_counts.df <- data.frame(N_CellsNoZeroCounts=num_cells_no_zero_counts, N_CellsLowExpr=num_cells_low_expr, LowExprThres=low_expr_threshold)
} else {
  num_cells_no_zero_counts.df <- data.frame(N_CellsNoZeroCounts=num_cells_no_zero_counts)
}



# Add Sample ID to cells before merging seurat objects
print("Adding sample name to cell identifiers and metadata...")

# Add info to metadata
xenium.obj$sample <- sample # sample id
xenium.obj$sample_name <- sample_name # sample name
xenium.obj$strain <- sample_info$Strain[which(sample_info$Name==sample)] # strain
xenium.obj$cond <- sample_info$Cond[which(sample_info$Name==sample)] # condition (WT)
xenium.obj$barcodes <- rownames(xenium.obj@meta.data) # cell barcodes

# Add sample name to cell identifiers (rownames of metadata)
colnames(xenium.obj) <- paste0(xenium.obj$barcodes,"__",sample,"__",sample_name)

# Add new cell identifiers with sample names to metadata
xenium.obj$cell.ids <- rownames(xenium.obj@meta.data)

# Add an X before gene names, since many gene names start with a number and give error in some plots
rownames(xenium.obj@assays$Xenium) <- paste0("X",rownames(xenium.obj@assays$Xenium)) %>% sub("X*","X",.)

# Do violin plot of genes per cell and transcript counts per cell
vioplot <- VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 1, pt.size = 0)+
  patchwork::plot_annotation(title = sample_name_full)

### ===== DO IMAGE PLOTS WITH MARKER GENES ===== ###

# Define markers for each cell type... add an X
genes_mt <- c("mt-Co2","mt-Atp6")
genes_piwi <- c("Piwil1","Piwil2")
genes_general <- c("Ddx4","Gata4")
genes_sertoli <- c("Clu","Ctsl","Rhox8","Amhr2")
genes_spgonia <- c("Sycp1","Crabp1","Stra8","Uchl1")
genes_sptcyte <- c("Piwil1","Insl6","Tbpl1")
genes_sptcyte2 <- c("Insl6","Tbpl1")
genes_rsptid  <- c("Tsga8","Acrv1","Spaca1","Tssk1")
genes_esptid  <- c("1700015G11Rik","Spata3")

# Do dimensional plot with the axis only
img_dim_plot <- Seurat::ImageDimPlot(xenium.obj, axes = T, molecules = "", cols = "white", nmols = 5000, group.by = "orig.ident") + ggmitji::add_grid(.5,.5,.5,.5,"red") + labs(title = sample_name_full)

# Stablish default boundary as segmentation
SeuratObject::DefaultBoundary(xenium.obj[["fov"]]) <- "segmentation"

### ===== NORMALIZATION AND CLUSTERING ====== ###

# Normalizaton, UMAP and clustering
xenium.obj <- Seurat::NormalizeData(xenium.obj, assay = "Xenium", normalization.method = "LogNormalize")
xenium.obj <- Seurat::FindVariableFeatures(xenium.obj)
xenium.obj <- Seurat::ScaleData(xenium.obj)
xenium.obj <- Seurat::RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- Seurat::RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- Seurat::FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
for (res in c(.1,.3,.5,1)){ xenium.obj <- Seurat::FindClusters(xenium.obj, resolution = res) }

# Do dim plot with umap
dim_plot_umap <- Seurat::DimPlot(xenium.obj, group.by = "Xenium_snn_res.0.5", reduction = "umap", label = T) + labs(title=sample_name_full, caption="Resolution = 0.5")

# Define neigboring resolution in object's identity
Seurat::Idents(xenium.obj) <- "Xenium_snn_res.0.5"

# Find markers for all clusters
all_markers <- Seurat::FindAllMarkers(object = xenium.obj, only.pos = T)
all_markers <- all_markers %>% as.data.frame() %>% dplyr::mutate(gene=sub("^X","",gene), sample=sample, sample_name=sample_name)

### ===== DO FEATURE PLOTS on UMAP ===== ###

# Do feature plot for each marker set
# ftr_plot_mito      <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_mt), max.cutoff = "q99", min.cutoff = "q10", combine = T, ncol = 1) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Mitochondrial marker genes", caption = "max.cutoff = q99, min.cutoff = q10")
ftr_plot_piwi      <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_piwi), max.cutoff = "q90", min.cutoff = "q10", combine = T, ncol = 1) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Piwi genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_soma_germ <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_general), max.cutoff = "q90", min.cutoff = "q10", combine = T, ncol = 1) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Somatic/Germ cell genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_sertoli   <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_sertoli), max.cutoff = "q90", min.cutoff = "q10", combine = T) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Sertoli cell marker genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_spgonia   <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_spgonia), max.cutoff = "q90", min.cutoff = "q10", combine = T) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Spermatogonia marker genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_sptcyte   <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_sptcyte), max.cutoff = "q90", min.cutoff = "q10", combine = T) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Spermatocyte marker genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_rsptid    <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_rsptid), max.cutoff = "q90", min.cutoff = "q10", combine = T) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Round spermatid marker genes", caption = "max.cutoff = q90, min.cutoff = q10")
ftr_plot_esptid    <- Seurat::FeaturePlot(xenium.obj, features = paste0("X",genes_esptid), max.cutoff = "q90", min.cutoff = "q10", combine = T, ncol = 1) + patchwork::plot_annotation(title = sample_name_full, subtitle = "Elognating spermatid marker genes", caption = "max.cutoff = q90, min.cutoff = q10")


### ===== SAVE TO FILE ===== ###


# Create directory for seurat objects that will be updated on each analysis analysis
outdir_objects <- paste0("5nxf1_iap_B6CAST/3xenium_B6CAST/output/samples/seurat_objects/")
dir.create(paste0(outdir_objects),F,T)

# Create directories for initial pre-integration analysis 
outdir <- paste0("5nxf1_iap_B6CAST/3xenium_B6CAST/output/samples/pre_integration/",sample)
dir.create(paste0(outdir,"/plots"),F,T)
# dir.create(paste0(outdir,"/plots/imageDimPlots/"),F,T)
# dir.create(paste0(outdir,"/plots/imageFeaturePlots/"),F,T)

# Save Xenium object
saveRDS(xenium.obj, paste0(outdir_objects,"/",sample,".seuratObject.RDS"))

# Save all markers
data.table::fwrite(all_markers,paste0(outdir,"/",sample,".all_markers.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

# Save number of cells
num_cells <- data.frame(Sample=sample, SampleName=sample, N_TotalCells=total_num_cells) %>% dplyr::bind_cols(num_cells_no_zero_counts.df)
data.table::fwrite(num_cells,paste0(outdir,"/",sample,".nCells.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

# Save violin plot
ggsave(filename = paste0(outdir,"/plots/",sample,".violinPlot.png"),  plot = vioplot, width = 7, height = 4, units = "in")

# Save UMAP 
ggsave(filename = paste0(outdir,"/plots/",sample,".UMAP.png"),  plot = dim_plot_umap, width = 7, height = 7, units = "in")

# Save feature plots with markers on UMAP
pdf(paste0(outdir,"/plots/",sample,".FeaturePlots.pdf"))
dim_plot_umap
ftr_plot_soma_germ
ftr_plot_piwi
ftr_plot_sertoli
ftr_plot_spgonia
ftr_plot_sptcyte
ftr_plot_rsptid
ftr_plot_esptid
dev.off()