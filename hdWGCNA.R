library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 9)
seurat_obj <- readRDS('./Pec_final.rds')

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
   umap_theme() + ggtitle('Zhou et al Control Cortex') + NoLegend()

p

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "tutorial" 
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"), 
  reduction = 'harmony', 
  k = 25, 
  max_shared = 10, 
  ident.group = 'cell_type' 
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "INH", 
  group.by='cell_type', 
  assay = 'RNA', 
  slot = 'data' 
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'INH' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='INH hdWGCNA Dendrogram')

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="Sample"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'Plume'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "INH-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA_object.rds')
