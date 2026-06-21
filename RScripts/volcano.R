# 火山图
vol_plot = ggVolcano::ggvolcano(
  data = DEG,
  x = "log2FoldChange",
  y = "padj",
  legend_position = "UR",
  label = "gene_symbol",
  label_number = 10,
  log2FC_cut = logFC_cutoff,
  fills = c("#e94234", "#b4b4d8", "#269846"),
  colors = c("#e94234", "#b4b4d8", "#269846"),
  output = TRUE,
  filename = glue("{plot_dir}/{des}_vol")
)
# vol_plot