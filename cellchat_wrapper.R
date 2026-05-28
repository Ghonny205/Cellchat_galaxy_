suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(future)
})

# Silenciamos todos los avisos para que Galaxy no bloquee las cajas
options(warn = -1) 
options(future.globals.maxSize = 8000 * 1024^2)

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input Seurat RDS"),
  make_option(c("-s", "--species"), type="character", default="human", help="human/mouse"),
  make_option(c("-g", "--group_by"), type="character", default="ident", help="Columna de tipos celulares"),
  make_option(c("-o", "--output_table"), type="character", default="interactions.tsv"),
  make_option(c("-p", "--output_plot"), type="character", default="network.pdf"),
  make_option(c("-r", "--output_rds"), type="character", default="cellchat_result.rds")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!file.exists(opt$input)) stop("El archivo de entrada no existe.")
seurat_obj <- readRDS(opt$input)

# 1. Aseguramos que los datos estén normalizados (si no lo están, CellChat explota)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
ensayo <- DefaultAssay(seurat_obj)

# 2. EXTRAEMOS LA MATRIZ (El "Hack" definitivo)
# Así evitamos que CellChat intente leer el objeto Seurat y falle por versiones
data.input <- tryCatch({
  GetAssayData(seurat_obj, assay = ensayo, layer = "data")
}, error = function(e) {
  GetAssayData(seurat_obj, assay = ensayo, slot = "data")
})

# 3. EXTRAEMOS LOS METADATOS
meta <- seurat_obj@meta.data
if (opt$group_by != "ident") {
  if (!opt$group_by %in% colnames(meta)) {
    stop(paste("La columna", opt$group_by, "no existe en los metadatos."))
  }
  meta$cellchat_labels <- meta[[opt$group_by]]
} else {
  meta$cellchat_labels <- Idents(seurat_obj)
}
meta$cellchat_labels <- as.factor(meta$cellchat_labels)

if (length(levels(meta$cellchat_labels)) < 2) {
  stop("Se necesitan al menos 2 tipos celulares para analizar comunicación.")
}

# 4. CREAMOS CELLCHAT DESDE LA MATRIZ EN CRUDO (¡Infalible!)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_labels")

# 5. Asignamos la base de datos
if (tolower(opt$species) == "mouse") {
  cellchat@DB <- CellChatDB.mouse
} else {
  cellchat@DB <- CellChatDB.human
}

# 6. Ejecutamos todo el pipeline matemático
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 7. Guardamos los tres resultados para Galaxy
df.net <- subsetCommunication(cellchat)
write.table(df.net, file = opt$output_table, sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(cellchat, file = opt$output_rds)

pdf(opt$output_plot, width = 12, height = 12)
par(mfrow = c(1,1))
netVisual_circle(cellchat@net$count, weight.scale = TRUE, label.edge = FALSE, 
                 title.name = paste("Red de comunicación -", opt$species))
dev.off()

message("Proceso finalizado con éxito.")