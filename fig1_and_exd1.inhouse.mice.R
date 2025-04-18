### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(stringr)
library(scater)
library(DropletUtils)
library(harmony)
library(pheatmap)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/2.AT2_NFkB/outs/2.1.AT2_PND1vsPND7/'
so <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/1.normalization_and_annotation/outs/normalization_with_detected_filter/so.f.annotated.rds')
so.AT2 <- subset(so, cells = Cells(so)[so$annotation == 'AT2'])
so.AT2 <- log_normalize(so.AT2, save_path, 'AT2', 5, 1000)
# GO
# progeny
# scenic

### 0. DEG
Idents(so.AT2) <- so.AT2$PND
deg <- FindMarkers(so.AT2, ident.1 = 'PND7', ident.2 = 'PND1', logfc.threshold = 0, min.pct = 0)

### 1. GO
library(GO.db)
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")

# load hallmark geneset
geneset_name = 'hallmark'
deg.f <- subset(deg, p_val_adj < 0.05)

hallmark_file = '/home/sjcho/datas/public_data/geneset/mh.all.v2023.2.Mm.symbols_hamllmark_genesets.gmt'
token = upload_GMT_file(gmtfile = hallmark_file)
# Your custom annotations ID is gp__52Qh_oRyW_jfk.
# You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
# Just use: gost(my_genes, organism = 'gp__52Qh_oRyW_jfk')
gp = gost(list("up-regulated" = subset(deg.f, avg_log2FC > 0) %>% rownames, "down-regulated" = subset(deg.f, avg_log2FC < 0) %>% rownames), organism = "mmusculus")
gp_bp = subset(gp$result, query == 'up-regulated')
gp[gp$source == 'GO:BP', ]

custom_gp = gost(list("up-regulated" = subset(deg.f, avg_log2FC > 0) %>% rownames, "down-regulated" = subset(deg.f, avg_log2FC < 0) %>% rownames), 
            organism = "gp__52Qh_oRyW_jfk",
            significant = FALSE)
gp_result <- custom_gp$result
gp_up <- gp_result[gp_result$query == 'up-regulated', ]
# gp_result[['log_p']] <- ifelse(gp_result$query == 'up-regulated', -log10(gp_result$p_value), log10(gp_result$p_value))
gp_up[['log_p']] <- -log10(gp_up$p_value)
gp_up$term_id <- factor(gp_up$term_id, levels = unique(gp_result$term_id))

# query is categorical value
# fill the bars by query
p <- ggplot(gp_up[1:4, ], aes(x = reorder(term_id, log_p), y = log_p, fill = query)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c('up-regulated' = "#d72b2b", 'down-regulated' = "darkblue")) +
    labs(title = "", x = "pathway", y = "-log10(p value)") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    ) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
p <- p + theme_classic() + coord_flip()
p <- p + theme(
    axis.title = element_text(color = "black"),
    axis.text  = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )
ggsave(p, file= paste0(save_path, paste0('go_HM_PND7_enriched.png')), width = 8, height = 6)
ggsave(p, file= paste0(save_path, paste0('go_HM_PND7_enriched.pdf')), width = 8, height = 6)
ggplot2pptx(p, file= paste0(save_path, paste0('go_HM_PND7_enriched.pptx')), width = 8, height = 6)

### 2. progeny
# progeny
save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/2.AT2_NFkB/outs/2.1.AT2_PND1vsPND7/progeny/'

library(decoupleR)
library(tidyr)
library(tibble)
library(ggsignif)
net <- get_progeny(organism = 'mouse') # top 500 genes -> default
# use the top 500 responsive genes ranked by p-value

# Extract the normalized log-transformed counts
mat <- as.matrix(so.AT2@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
so.AT2[['pathways_mlm']] <- acts %>%
pivot_wider(id_cols = 'source', names_from = 'condition',
            values_from = 'score') %>%
            column_to_rownames('source') %>%
            Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = so.AT2) <- "pathways_mlm"

# Scale the data
so.AT2 <- ScaleData(so.AT2)
so.AT2@assays$pathways_mlm@data <- so.AT2@assays$pathways_mlm@scale.data

# p1 <- DimPlot(so.AT2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = row) + 
# NoLegend() + ggtitle(row)
# p2 <- (FeaturePlot(so.AT2, features = c("Trail")) & 
# scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
# ggtitle('Trail activity')
# p1 | p2

# Extract activities from object as a long dataframe
df <- t(as.matrix(so.AT2@assays$pathways_mlm@data)) %>%
    as.data.frame() %>%
    mutate(Condition = so.AT2[['PND']][[1]]) %>%
    pivot_longer(cols = -Condition, names_to = "source", values_to = "score") %>%
    group_by(Condition, source) %>%
    summarise(mean = mean(score))

p1 <- VlnPlot(so.AT2, features = c("NFkB"), pt.size = 0, group.by = 'PND', cols = c('gray60', '#d72b2b')) + NoLegend() + geom_boxplot(fill = 'white', width = 0.2)
p1 <- p1 + ylim(so.AT2@assays$pathways_mlm@data['NFkB', ] %>% min, so.AT2@assays$pathways_mlm@data['NFkB', ] %>% max * 1.15)
p1 <- p1 + geom_signif(comparisons = list(c('PND1', 'PND7')), map_signif_level = TRUE, textsize = 3, vjust = -0.5)
p1 <- p1 + ylab('Pathway activity')
# remove x title and text
p1 <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
ggsave(p1, file = paste0(save_path, 'progeny_NFkB_vlnplot.png'), width = 4, height = 4)
ggsave(p1, file = paste0(save_path, 'progeny_NFkB_vlnplot.pdf'), width = 4, height = 4)
ggplot2pptx(p1, file = paste0(save_path, 'progeny_NFkB_vlnplot.pptx'), width = 4, height = 4)

# signature score

DefaultAssay(so.AT2) <- 'RNA'
so.AT2 <- AddModuleScore(so.AT2, features = list(subset(h_gene_sets, gs_name == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')$gene_symbol), name = 'HM_NFkB')
p2 <- VlnPlot(so.AT2, features = c('HM_NFkB1'), pt.size = 0, group.by = 'PND', cols = c('gray60', '#d72b2b')) + NoLegend() + geom_boxplot(fill = 'white', width = 0.2)
p2 <- p2 + ylim(so.AT2$HM_NFkB1 %>% min, so.AT2$HM_NFkB1 %>% max * 1.15)
p2 <- p2 + geom_signif(comparisons = list(c('PND1', 'PND7')), map_signif_level = TRUE, textsize = 3, vjust = -0.5)
p2 <- p2 + ylab('Signature score')
ggsave(p2, file = paste0(save_path, 'HM_NFkB_vlnplot.png'), width = 4, height = 4)
ggsave(p2, file = paste0(save_path, 'HM_NFkB_vlnplot.pdf'), width = 4, height = 4)
ggplot2pptx(p2, file = paste0(save_path, 'HM_NFkB_vlnplot.pptx'), width = 4, height = 4)

p <- p1/p2
ggsave(p, file = paste0(save_path, 'progeny_HM_NFkB.png'), width = 4, height = 6)
ggsave(p, file = paste0(save_path, 'progeny_HM_NFkB.pdf'), width = 4, height = 6)
ggplot2pptx(p, file = paste0(save_path, 'progeny_HM_NFkB.pptx'), width = 4, height = 6)

# nichenet

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/2.AT2_NFkB/outs/2.2.nichenet/'

so <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/1.normalization_and_annotation/outs/normalization_with_detected_filter/so.f.annotated.rds')
so.7 <- subset(so, PND == 'PND7')

### 3. nichenet

# based on https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
# files from https://zenodo.org/records/7074291
library(nichenetr)

organism <- "mouse"

lr_network <- readRDS('/home/sjcho/datas/CCI_db/nichenet_db/lr_network_mouse_21122021.rds')
ligand_target_matrix <- readRDS('/home/sjcho/datas/CCI_db/nichenet_db/ligand_target_matrix_nsga2r_final_mouse.rds')
weighted_networks <- readRDS('/home/sjcho/datas/CCI_db/nichenet_db/weighted_networks_nsga2r_final_mouse.rds')

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

### 1. Define a set of potential ligands for both the sender-agnostic and sender-focused approach
Idents(so.7) = so.7$annotation

receiver = "AT2"
expressed_genes_receiver <- get_expressed_genes(receiver, so.7, pct = 0.1) # 0.1 is default value

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

## sender-foucsed approach
sender_celltypes <- Idents(so.7) %>% unique

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, so.7, 0.15)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
length(expressed_genes_sender)
## [1] 14139
length(potential_ligands)
## [1] 778
length(potential_ligands_focused) # only the ligands that are expressed in the sender cell types
## [1] 449

### 2. Define the gene set of interest
so.AT2 <- subset(so, cells = Cells(so)[so$annotation == 'AT2'])
so.AT2 <- log_normalize(so.AT2, save_path, 'AT2', 5, 1000)
Idents(so.AT2) <- so.AT2$PND
deg <- FindMarkers(so.AT2, ident.1 = 'PND7', ident.2 = 'PND1', logfc.threshold = 0, min.pct = 0)
deg.7.enriched <- deg %>% filter(p_val_adj <= 0.05 & avg_log2FC > log2(1)) # same deg with GO analysis

geneset_oi <- rownames(deg.7.enriched)
### 3. Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
# All expressed genes in the receiver cell population (that are also in the ligand-target matrix) is defined as the ‘background set’
# generally 5000 ~ 10000 genes
length(background_expressed_genes)
## [1] 7983
length(geneset_oi)
## [1] 92

### 4. Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands_focused)
save(ligand_activities, file = paste0(save_path, 'ligand_activities.Rdata'))
p <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange", bins = 100, linewidth = 0.1)  + 
  geom_vline(xintercept = 0.025, color="black", linetype="dashed", size=1.5) + 
  labs(x="ligand activity (corrected AUPR)", y = "# of ligands") +
  theme_classic()
ggsave(p, file = paste0(save_path, 'exd2.X.ligand_activities_hist.png'), width = 6, height = 6)
ggsave(p, file = paste0(save_path, 'exd2.X.ligand_activities_hist.pdf'), width = 6, height = 6)
ggplot2pptx(p, file = paste0(save_path, 'exd2.X.ligand_activities_hist.pptx'), width = 6, height = 6)

ligands_top <- subset(ligand_activities, aupr_corrected > 0.025)$test_ligand
save(ligands_top, file = paste0(save_path, 'ligands_top.Rdata'))
so.7 <- subset(so, PND == 'PND7')
so.7 <- log_normalize(so.7, save_path, 'PND7', 10, 2000)
so.7 <- AddModuleScore(so.7, features = list(ligands_top), name = 'ligand_activity')
# , 
p <- VlnPlot(so.7, features = c('ligand_activity1'), group.by = 'annotation2', cols = PND.cols, pt.size = 0) + geom_boxplot(fill = 'white', width = 0.2) + NoLegend()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + ylab('Ligand activity signature score') + xlab('') + ggtitle('')
ggsave(p, file = paste0(save_path, 'test.png'), width = 6, height = 6)
ggsave(p, file = paste0(save_path, 'fig2.C.ligand_activity_vln.png'), width = 6, height = 6)
ggsave(p, file = paste0(save_path, 'fig2.C.ligand_activity_vln.pdf'), width = 6, height = 6)
ggplot2pptx(p, file = paste0(save_path, 'fig2.C.ligand_activity_vln.pptx'), width = 6, height = 6)

### plot ligand activity
library(RColorBrewer)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% ligands_top) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
vis_ligand_aupr_rev <- as.data.frame(rev(vis_ligand_aupr)); rownames(vis_ligand_aupr_rev) <- rev(rownames(vis_ligand_aupr)); colnames(vis_ligand_aupr_rev) <- colnames(vis_ligand_aupr)
p1 <- make_heatmap_ggplot(vis_ligand_aupr_rev %>% as.matrix(ncol = 1),
                     "Prioritized ligands", "", 
                     legend_title = "ligand activity\n(corrected AUPR)", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"))
ggsave(p1, file = paste0(save_path, 'exd2.X.ligand_activities_heatmap.png'), width = 2, height = 6)

p2 <- DotPlot(so.7, features = rownames(vis_ligand_aupr) %>% rev(), group.by = 'annotation2', assay = 'RNA', dot.scale = 5) + theme_classic()
p2 <- p2 + theme(axis.text.x = element_text(size = 0),
                axis.title = element_text(size = 0))
p2 <- p2 + scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")))
p2 <- p2 + theme(axis.text = element_text(color = 'black'),
                axis.title = element_text(color = 'black'))

p <- egg::ggarrange(p2, p1 + coord_flip(), ncol = 1, heights = c(10, 1))
ggsave(p, file = paste0(save_path, 'exd2.X.top.priority.ligands.expression.png'), width = 9, height = 8)
ggsave(p, file = paste0(save_path, 'exd2.X.top.priority.ligands.expression.pdf'), width = 9, height = 8)
ggplot2pptx(as.ggplot(p), file = paste0(save_path, 'exd2.X.top.priority.ligands.expression.pptx'), width = 9, height = 8)
