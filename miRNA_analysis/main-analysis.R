#setwd("//")

suppressMessages({
  library("tidyverse")
  library("ggpubr")
  library("pheatmap")
  library("matrixStats")
  library("readxl")
  library("data.table")
  library("svglite")
  library("factoextra")
  library("ggrepel")
  library("enrichplot")
  library("VennDiagram")
  library("ggvenn")
  library("igraph")
  library("ggraph")
  library("ComplexHeatmap")
  library("RColorBrewer")
  library("circlize")
  library("patchwork")
  library("scales")
  library("grid")
  library("openxlsx")
  library("purrr")
  library("variancePartition")
  library("limma")
  library("eulerr")
})

set.seed(42)


# ---------------------- Define Color Schemes (Global) ----------------------
point_colors <- c(Control_M = "#d4751a", Control_F = "#d4751a",
                  PLX_M = "#8b3a7a", PLX_F = "#8b3a7a")

treatment_colors <- c("Control" = "#d4751a", "PLX" = "#8b3a7a")

sex_colors <- c("Male" = "#254e7a", "Female" = "#a42e30")

region_colors <- c("Bregma -2.06" = "#b07d2a",
                   "Olfactory Bulb" = "#9e4d8c",
                   "Caudate Putamen" = "#2e7b82",
                   "Cerebellum" = "#2e8c61")

point_labels <- c("Control_M" = "Control M", "Control_F" = "Control F",
                  "PLX_M" = "PLX M", "PLX_F" = "PLX F")

ncRNA_colors <- c(
  "#a3834c",  # tRNA
  "#a0cbe8",  # miRNA
  "#1D6B5A",  # snRNA
  "#bab0ac",  # snoRNA
  "#DDCC77",  # scaRNA
  "#4A4A4A",  # rRNA
  "#ff9999",  # miscRNA
  "#589f70"   # lncRNA
)

regions <- list(S = "Bregma -2.06", OB = "Olfactory Bulb",
                CP = "Caudate Putamen", CE = "Cerebellum")

# ---------------------- Load data miRNA ----------------------
### Depletion study data load
data <- read.table("Input files/sncRNA_complete_quantification_rpmm_norm.csv",
                   header = TRUE, sep = "\t")

sample_names <- c(
  "RNA",
  "S15", "S4", "S5", "S6", "S7", "S8", "S19", "S20", "S16",
  "S9", "S10", "S11", "S17", "S12", "S13", "S14", "S18", "S1", "S2", "S3",
  "OB17", "OB4", "OB5", "OB8", "OB9", "OB10", "OB18", "OB11", "OB12",
  "OB13", "OB19", "OB14", "OB15", "OB16", "OB20", "OB1", "OB2", "OB3",
  "CP16", "CP4", "CP5", "CP6", "CP7", "CP8", "CP20", "CP17", "CP10", "CP11",
  "CP12", "CP13", "CP18", "CP14", "CP15", "CP19", "CP1", "CP2", "CP3",
  "CE15", "CE4", "CE5", "CE6", "CE7", "CE8", "CE19", "CE20", "CE16",
  "CE9", "CE10", "CE11", "CE12", "CE13", "CE17", "CE14", "CE18",
  "CE1", "CE2", "CE3")

colnames(data) <- sample_names
data_mmu <- data[grep("mmu-miR-", data$RNA), ]

miRNAs_of_interest <- c("mmu-miR-146a-5p", "mmu-miR-223-3p")

diff_express_gender_specific_75per <- read.table(
  "snc-pipeline/plx_microglia_gender_specific_75%/diff_exp/miRNA/rpmm_norm_log2.detection_rate_75p_per_group/diff_exp_gender_specific_75per.csv",
  header = TRUE, sep = "\t"
)

# Pomper database
xlsx_pomper_database <- "Input files/LoessM_norm_miRNA_array_WEBSITE_final.xlsx"
out_dir   <- "Figures and plots/miRNA_brainstem_plots"
bar_color <- "#2E5560"

# External validation datasets
diff_exp_external <- read.table(
  "Input files/results_plx_walsh_microglia/diff_exp/miRNA/rpmm_norm_log2.detection_rate_75p_per_group/diff_exp.csv",
  header = TRUE, sep = "\t", check.names = FALSE
)

strip_mmu <- function(x) sub("^mmu-", "", x)

# ---------------------- Sample Groups Definition ----------------------
### Sample-to-group (treatment x sex) assignment within each region
all_samples_columns_small_rna <- list(
  Control_M = c("S17", "S18", "S19", "S20"),
  PLX_M = c("S7", "S8", "S9", "S10", "S11", "S12"),
  Control_F = c("S13", "S14", "S15", "S16"),
  PLX_F = c("S1", "S2", "S3", "S4", "S5", "S6")
)

all_samples_columns_olfactory_bulb <- list(
  Control_M = c("OB17", "OB18", "OB19", "OB20"),
  PLX_M = c("OB8", "OB9", "OB10", "OB11", "OB12"),
  Control_F = c("OB13", "OB14", "OB15", "OB16"),
  PLX_F = c("OB1", "OB2", "OB3", "OB4", "OB5")
)

all_samples_columns_caudate_putamen <- list(
  Control_M = c("CP17", "CP18", "CP19", "CP20"),
  PLX_M = c("CP7", "CP8", "CP10", "CP11", "CP12"),
  Control_F = c("CP13", "CP14", "CP15", "CP16"),
  PLX_F = c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6")
)

all_samples_columns_cerebellum <- list(
  Control_M = c("CE17", "CE18", "CE19", "CE20"),
  PLX_M = c("CE7", "CE8", "CE9", "CE10", "CE11", "CE12"),
  Control_F = c("CE13", "CE14", "CE15", "CE16"),
  PLX_F = c("CE1", "CE2", "CE3", "CE4", "CE5", "CE6")
)

treatment_levels <- c("Control_M", "PLX_M", "Control_F", "PLX_F")
comparison_list_mirna_evs <- list(c("Control_M", "PLX_M"), c("Control_F", "PLX_F"))

### Helper function to process the brain regions
process_brain_region <- function(data, sample_groups) {
  for (name in names(sample_groups)) {
    cols <- sample_groups[[name]]
    data[[paste0(name, "_Mean")]] <- rowMeans(data[, cols, drop = FALSE])
    data[[paste0(name, "_SD")]] <- rowSds(as.matrix(data[, cols, drop = FALSE]))
  }
  return(data)
}

# ---------------------- QC Brain region correlation analysis ----------------------
corr_data <- read.csv("Input files/sample_correlation.spearman.csv", sep = "\t")

get_region <- function(sample_id) {
  if (grepl("^S[0-9]", sample_id))  return("Bregma -2.06")
  if (grepl("^OB", sample_id))      return("Olfactory Bulb")
  if (grepl("^CP", sample_id))      return("Caudate Putamen")
  if (grepl("^CE", sample_id))      return("Cerebellum")
  "Unknown"
}

assign_treatment_sex <- function(region, sample_num) {
  treatment <- if (sample_num <= 12) "PLX" else "Control"
  sex <- if (sample_num <= 6 || (sample_num >= 13 && sample_num <= 16)) "Female" else "Male"
  c(treatment, sex)
}

get_sample_num <- function(sample_id) {
  as.numeric(gsub("^[A-Z]+", "", sample_id))
}

# Sanity check on group assignment
for (range_label in list(c(1, 6, "PLX Female"), c(7, 12, "PLX Male"),
                         c(13, 16, "Control Female"), c(17, 20, "Control Male"))) {
  cat(sprintf("Sample %s-%s (expect %s):\n", range_label[1], range_label[2], range_label[3]))
  for (i in seq(as.integer(range_label[1]), as.integer(range_label[2]))) {
    r <- assign_treatment_sex("Any", i)
    cat(sprintf("  %d: %s %s\n", i, r[1], r[2]))
  }
}

all_samples <- colnames(corr_data)
sample_regions <- sapply(all_samples, get_region)

corr_matrix_full <- as.matrix(corr_data)

# Color scheme
col_fun <- colorRamp2(c(0.85, 0.90, 0.95, 1.0),
                      c("#9ECAE1", "#4292C6", "#08519C", "#08306B"))

all_metadata <- data.frame(
  Sample = all_samples,
  Region = sample_regions,
  Sample_Num = sapply(all_samples, get_sample_num),
  stringsAsFactors = FALSE
)

all_metadata$Treatment <- NA
all_metadata$Sex <- NA

for (i in 1:nrow(all_metadata)) {
  assignment <- assign_treatment_sex(all_metadata$Region[i], all_metadata$Sample_Num[i])
  all_metadata$Treatment[i] <- assignment[1]
  all_metadata$Sex[i] <- assignment[2]
}

# Top annotation
ha_top_all <- HeatmapAnnotation(
  Region = all_metadata$Region,
  Treatment = all_metadata$Treatment,
  Sex = all_metadata$Sex,
  col = list(
    Region = region_colors,
    Treatment = treatment_colors,
    Sex = sex_colors
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Region = list(
      title = "Region",
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 13)
    ),
    Treatment = list(
      title = "Treatment",
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 13)
    ),
    Sex = list(
      title = "Sex",
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 13)
    )
  )
)

# Left annotation
ha_left_all <- rowAnnotation(
  Sex = all_metadata$Sex,
  Treatment = all_metadata$Treatment,
  Region = all_metadata$Region,
  col = list(
    Sex = sex_colors,
    Treatment = treatment_colors,
    Region = region_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Create heatmap with all regions
ht_all <- Heatmap(
  corr_matrix_full,
  name = "Spearman\nCorrelation",
  col = col_fun,
  column_title = "Brain Regions - Sample Correlation",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_top_all,
  left_annotation = ha_left_all,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Spearman\nCorrelation",
    title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica"),
    direction = "vertical",
    legend_height = unit(6, "cm"),
    at = c(0.85, 0.90, 0.95, 1.0)
  )
)

#pdf("Figures and plots/PDFs/All_Regions_Single_correlation_heatmap.pdf", width = 11, height = 13)
#draw(ht_all, heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()

# ---------------------- QC Number of mapped lncRNAs ----------------------
qc_data <- read.csv("Input files/qc_table.biol.csv", sep = "\t")

ncRNA_names <- c("tRNA", "miRNA", "snRNA", "snoRNA", "scaRNA", "rRNA", "miscRNA", "lncRNA")
qc_data$Region <- sapply(qc_data$Sample_ID, get_region)

# Extract numeric part from Sample_ID
get_sample_number <- function(sample_id) {
  num <- gsub("^[A-Z]+", "", sample_id)
  as.numeric(num)
}

qc_data$Sample_Number <- sapply(qc_data$Sample_ID, get_sample_num)

### Assign treatment and sex
assignments <- t(mapply(assign_treatment_sex, qc_data$Region, qc_data$Sample_Number))
qc_data$Treatment <- assignments[, 1]
qc_data$Sex <- assignments[, 2]
qc_data$Group <- paste(qc_data$Treatment, qc_data$Sex, sep = "_")

ncRNA_vars <- c(
  "total_mapped_tRNA",
  "total_mapped_to_precursor_or_miRNA",
  "total_mapped_snRNA",
  "total_mapped_snoRNA",
  "total_mapped_scaRNA",
  "total_mapped_rRNA",
  "total_mapped_miscRNA",
  "total_mapped_lncRNA"
)

ncRNA_matrix <- as.matrix(qc_data[, ncRNA_vars])
rownames(ncRNA_matrix) <- qc_data$Sample_ID

# Build a 20-position (20 samples) layout for a region
create_complete_region_data <- function(region_name, qc_data, ncRNA_matrix) {
  region_idx <- which(qc_data$Region == region_name)
  region_qc <- qc_data[region_idx, ]
  
  complete_data <- data.frame(
    Sample_Number = 1:20,
    Treatment = NA,
    Sex = NA,
    stringsAsFactors = FALSE
  )
  
  complete_matrix <- matrix(0, nrow = ncol(ncRNA_matrix), ncol = 20)
  rownames(complete_matrix) <- colnames(ncRNA_matrix)
  colnames(complete_matrix) <- paste0("Sample", 1:20)
  
  for (i in 1:nrow(region_qc)) {
    sample_num <- region_qc$Sample_Number[i]
    complete_data$Treatment[sample_num] <- region_qc$Treatment[i]
    complete_data$Sex[sample_num] <- region_qc$Sex[i]
    complete_matrix[, sample_num] <- ncRNA_matrix[region_idx[i], ]
  }
  
  return(list(metadata = complete_data, matrix = complete_matrix))
}

# Shared y-axis across all regions
global_ymax <- 0
for (region in regions) {
  idx <- which(qc_data$Region == region)
  if (length(idx) > 0) {
    region_data <- t(ncRNA_matrix[idx, , drop = FALSE])
    region_max <- max(colSums(region_data))
    global_ymax <- max(global_ymax, region_max)
  }
}

global_ymax <- global_ymax * 1.15
cat("Global y-axis maximum:", global_ymax, "\n")

#pdf("Figures and plots/PDFs/ncRNA_mapping_by_region_improved.pdf", width = 8, height = 11)

layout(matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1),
       heights = c(1, 1, 1, 1, 0.25))

par(mar = c(4, 6, 3, 2),
    oma = c(0, 0, 0, 0),
    xpd = FALSE,
    family = "Helvetica")

for (i in seq_along(regions)) {
  
  region <- regions[[i]]
  
  complete_region <- create_complete_region_data(region, qc_data, ncRNA_matrix)
  region_data <- complete_region$matrix
  region_meta <- complete_region$metadata
  
  # Missing samples get a gray placeholder
  empty_bars <- which(is.na(region_meta$Treatment))
  filled_bars <- which(!is.na(region_meta$Treatment))
  
  test_bp <- barplot(region_data[1,], plot = FALSE, space = 0.35)
  bar_width <- diff(test_bp)[1] / 1.35
  
  plot.new()
  plot.window(xlim = range(test_bp) + c(-0.5, 0.5), ylim = c(0, global_ymax))
  
  # Horizontal gridlines every 10 million reads
  y_interval <- 1e+07
  y_gridlines <- seq(0, ceiling(global_ymax / y_interval) * y_interval, by = y_interval)
  abline(h = y_gridlines, col = "gray90", lty = "solid", lwd = 0.5)
  
  title(main = region, cex.main = 1.5, font.main = 2, line = 1.5)
  axis(2, las = 1, cex.axis = 1)
  mtext("Total Mapped Reads", side = 2, line = 4.5, cex = 1, font = 2)
  
  # Overlay the stacked barplot on top of the gridlines
  par(new = TRUE)
  bp <- barplot(
    region_data,
    col = ncRNA_colors,
    border = NA,
    las = 1,
    space = 0.35,
    ylab = "",
    main = "",
    ylim = c(0, global_ymax),
    axes = FALSE,
    names.arg = rep("", ncol(region_data)),
    add = TRUE
  )
  
  axis(1, at = bp, labels = 1:20, las = 1, cex.axis = 0.8, tick = FALSE, line = -0.5)
  
  # Annotation strip positions, just below the x-axis (treatment row above sex row)
  annotation_bottom <- -0.05 * global_ymax
  treatment_bottom <- annotation_bottom
  treatment_top <- treatment_bottom + (0.03 * global_ymax)
  sex_bottom <- treatment_bottom - (0.05 * global_ymax)
  sex_top <- sex_bottom + (0.03 * global_ymax)
  
  par(xpd = TRUE)
  
  # Treatment color strip (filled bars only)
  for (bar_idx in filled_bars) {
    rect(
      xleft  = bp[bar_idx] - bar_width/2,
      xright = bp[bar_idx] + bar_width/2,
      ybottom = treatment_bottom,
      ytop    = treatment_top,
      col = treatment_colors[region_meta$Treatment[bar_idx]],
      border = NA
    )
  }
  
  # Sex color strip (filled bars only)
  for (bar_idx in filled_bars) {
    rect(
      xleft  = bp[bar_idx] - bar_width/2,
      xright = bp[bar_idx] + bar_width/2,
      ybottom = sex_bottom,
      ytop    = sex_top,
      col = sex_colors[region_meta$Sex[bar_idx]],
      border = NA
    )
  }
  
  # Gray out positions with no sample
  if (length(empty_bars) > 0) {
    for (bar_idx in empty_bars) {
      rect(
        xleft  = bp[bar_idx] - bar_width/2,
        xright = bp[bar_idx] + bar_width/2,
        ybottom = 0,
        ytop    = global_ymax,
        col = rgb(0.95, 0.95, 0.95, 0.5),
        border = "gray80",
        lwd = 0.5
      )
    }
  }
  
  par(xpd = FALSE)
}

### Legend panel at the bottom
par(mar = c(0.5, 1, 0.5, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))

par(xpd = TRUE)

legend(
  x = "right",
  inset = c(0, 0),
  legend = c("lncRNA", "miscRNA", "rRNA", "scaRNA",
             "snoRNA", "snRNA", "miRNA", "tRNA"),
  fill = rev(ncRNA_colors),
  title = "ncRNA Type",
  cex = 1.1,
  border = NA,
  bty = "n",
  horiz = FALSE,
  title.adj = 0,
  ncol = 4,
  title.cex = 1.4,
  title.font = 2
)

legend(
  x = "left",
  inset = c(0.35, 0),
  legend = names(sex_colors),
  fill = sex_colors,
  title = "Sex",
  cex = 1.1,
  border = NA,
  title.adj = 0,
  bty = "n",
  horiz = FALSE,
  title.cex = 1.4,
  title.font = 2
)

legend(
  x = "center",
  inset = c(-0.35, 0),
  legend = names(treatment_colors),
  fill = treatment_colors,
  title = "Treatment",
  cex = 1.1,
  border = NA,
  bty = "n",
  horiz = FALSE,
  title.adj = 0,
  title.cex = 1.4,
  title.font = 2
)

#dev.off()

# ---------------------- Process Brain Regions ----------------------
small_rna_main_results <- process_brain_region(data_mmu[, c(1, 2:21)],
                                               all_samples_columns_small_rna)
olfactory_bulb_results <- process_brain_region(data_mmu[, c(1, 22:39)],
                                               all_samples_columns_olfactory_bulb)
caudate_putamen_results <- process_brain_region(data_mmu[, c(1, 40:58)],
                                                all_samples_columns_caudate_putamen)
cerebellum_results <- process_brain_region(data_mmu[, c(1, 59:78)],
                                           all_samples_columns_cerebellum)

# ---------------------- Bar Plots - Small RNA (Bregma -2.06) ----------------------
# Define sample columns to keep for data
columns_of_interest_with_data_points_small_rna <- c(
  "RNA", "S15", "S4", "S5", "S6", "S7", "S8", "S19", "S20", "S16",
  "S9", "S10", "S11", "S17", "S12", "S13", "S14", "S18", "S1", "S2", "S3"
)

sample_to_treatment_small_rna <- do.call(rbind, lapply(names(all_samples_columns_small_rna), function(grp) {
  data.frame(Sample = all_samples_columns_small_rna[[grp]],
             Treatment = grp,
             stringsAsFactors = FALSE)
}))

# Re-shape data to long format for plotting
subset_df_for_box_plots_with_data_points_small_rna <- small_rna_main_results %>%
  filter(RNA %in% miRNAs_of_interest) %>%
  dplyr::select(all_of(columns_of_interest_with_data_points_small_rna)) %>%
  pivot_longer(cols = -RNA, names_to = "Sample", values_to = "Value") %>%
  left_join(sample_to_treatment_small_rna, by = "Sample") %>%
  mutate(Treatment = factor(Treatment, levels = treatment_levels))

df_with_outliers_small_rna <- subset_df_for_box_plots_with_data_points_small_rna %>%
  group_by(RNA, Treatment) %>%
  mutate(
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    is_outlier = Value < (Q1 - 1.5 * IQR) | Value > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup() %>%
  mutate(Treatment = factor(Treatment, levels = treatment_levels))

### Build one bar plot per miRNA
plots_with_stats_small_rna <- df_with_outliers_small_rna %>%
  group_split(RNA) %>%
  lapply(function(df) {
    gene <- unique(df$RNA)
    
    df$Treatment <- factor(df$Treatment, levels = treatment_levels)
    df$Treatment_Simple <- ifelse(grepl("Control", df$Treatment), "Control", "PLX")
    df$Treatment_Simple <- factor(df$Treatment_Simple, levels = c("Control", "PLX"))
    
    treatments_in_df <- unique(df$Treatment)
    comparisons_filtered <- comparison_list_mirna_evs[
      sapply(comparison_list_mirna_evs, function(x) all(x %in% treatments_in_df))
    ]
    
    ggplot(df, aes(x = Treatment, y = Value, fill = Treatment_Simple)) +
      geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), width = 0.6) +
      geom_errorbar(stat = "summary", fun.data = mean_se,
                    position = position_dodge(width = 0.9), width = 0.2) +
      geom_vline(xintercept = 2.5, linetype = "dashed", color = "#dddddd", linewidth = 0.7) +
      geom_jitter(color = "black",
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85),
                  alpha = 0.6, size = 2, show.legend = FALSE) +
      scale_fill_manual(values = treatment_colors, name = "Treatment") +
      ggsignif::geom_signif(
        comparisons = comparisons_filtered,
        test = "t.test",
        map_signif_level = function(p) sprintf("p = %.2g", p),
        tip_length = 0.01,
        textsize = 4
      ) +
      theme_minimal() +
      ylim(0, 22) +
      labs(title = gene, x = "", y = "rpmm normalized counts") +
      theme_minimal() +
      theme(text        = element_text(family = "Helvetica", size = 14, color = "black"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title    = element_text(size = 16, face = "bold"),
            legend.text     = element_text(size = 14),
            legend.box.spacing = unit(0.1, "cm"),
            plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title.y    = element_text(size = 14, face = "bold"),
            axis.text.y     = element_text(size = 12, color = "black"),
            axis.text.x     = element_text(size = 1, color = "white"),
            axis.ticks.x    = element_blank()
      )
  })

plots_with_stats_small_rna

#pdf("Figures and plots/miRNA_barplots_small_rna.pdf", width = 5, height = 7)
#for(i in seq_along(plots_with_stats_small_rna)) {
#   print(plots_with_stats_small_rna[[i]])
# }
#dev.off()

# ---------------------- Volcano plot and paired plot for raw and adjusted p-values ----------------------
### Function that creates the volcano plot
make_mirna_regional_plots <- function(diff_express_df,
                                      regions,
                                      pval_type = c("raw", "adj"),
                                      fc_threshold = 1,
                                      pval_threshold = 0.05,
                                      sex_colors,
                                      save_plots = FALSE) {
  
  pval_type <- match.arg(pval_type)
  pval_prefix <- if (pval_type == "raw") "ttest_rawp" else "ttest_adjp"
  pval_label  <- if (pval_type == "raw") "raw p-value" else "adj p-value"
  file_tag    <- if (pval_type == "raw") "rawp" else "adjp"
  
  ### Prepare region/sex data
  prepare_region_sex_data <- function(region_code, region_name) {
    build_df <- function(sex_code, sex_label) {
      fc_col   <- paste0("fc_Comparison__Treatment_", region_code, "_", sex_code,
                         "_vs_Control_", region_code, "_", sex_code)
      pval_col <- paste0(pval_prefix, "_Comparison__Treatment_", region_code, "_", sex_code,
                         "_vs_Control_", region_code, "_", sex_code)
      
      diff_express_df %>%
        dplyr::select(RNA, fc = all_of(fc_col), pval = all_of(pval_col)) %>%
        mutate(
          log2FC      = log2(fc),
          log10P      = -log10(pval),
          Sex         = sex_label,
          Significant = (abs(log2FC) > fc_threshold & pval < pval_threshold),
          Region      = region_name
        )
    }
    bind_rows(build_df("M", "Male"), build_df("F", "Female"))
  }
  
  all_region_data <- map_dfr(names(regions), function(rc) {
    prepare_region_sex_data(rc, regions[[rc]])
  })
  
  ### Volcano plot helper
  create_regional_volcano <- function(data, region_name, show_y_label = TRUE, show_y_text = TRUE) {
    region_data <- data %>%
      filter(Region == region_name) %>%
      mutate(
        plot_category = case_when(
          log10P < -log10(pval_threshold)        ~ "Not significant p-value",
          Significant                            ~ "Significant (both FC & p-value)",
          TRUE                                   ~ "Significant p-value only"
        )
      )
    
    p <- ggplot(region_data, aes(x = log2FC, y = log10P)) +
      geom_point(data = region_data %>% filter(plot_category == "Not significant p-value"),
                 color = "gray", size = 1.5, alpha = 0.4) +
      geom_point(data = region_data %>% filter(plot_category == "Significant p-value only"),
                 aes(color = Sex), size = 2, alpha = 0.3) +
      geom_point(data = region_data %>% filter(plot_category == "Significant (both FC & p-value)"),
                 aes(color = Sex), size = 3.5, alpha = 0.8) +
      scale_color_manual(values = sex_colors) +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray60") +
      theme_minimal() +
      ylim(0, 5) +
      xlim(-4, 4) +
      labs(
        x     = "log2(Fold Change)",
        y     = if (show_y_label) paste0("-log10(", pval_label, ")") else "",
        title = region_name,
        color = "Sex"
      ) +
      theme(
        text             = element_text(family = "Helvetica", size = 12, color = "black"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x      = element_text(size = 14, color = "black"),
        axis.text.y      = if (show_y_text) element_text(size = 14, color = "black") else element_blank(),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "none",
        panel.grid.minor = element_blank()
      )
    
    label_data <- region_data %>%
      filter(Significant) %>%
      mutate(RNA_clean = gsub("^mmu-", "", RNA))
    
    if (nrow(label_data) > 0) {
      p <- p + geom_text_repel(
        data = label_data, aes(label = RNA_clean),
        size = 4, fontface = "italic", color = "gray30",
        max.overlaps = Inf, min.segment.length = 0,
        box.padding = 1.0, point.padding = 0.8,
        force = 5, force_pull = 0.1,
        max.iter = 10000, max.time = 2,
        segment.color = "gray30", segment.size = 0.25, segment.alpha = 0.5,
        direction = "both", xlim = c(-5, 5), ylim = c(-0.5, 6),
        seed = 42, show.legend = FALSE
      )
    }
    p
  }
  
  ### Combined volcano
  plot_S  <- create_regional_volcano(all_region_data, "Bregma -2.06",    TRUE,  TRUE)
  plot_OB <- create_regional_volcano(all_region_data, "Olfactory Bulb",  FALSE, FALSE)
  plot_CP <- create_regional_volcano(all_region_data, "Caudate Putamen", FALSE, FALSE)
  plot_CE <- create_regional_volcano(all_region_data, "Cerebellum",      FALSE, FALSE)
  
  combined_regional_plot <- (plot_S | plot_CP | plot_CE | plot_OB) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Regional miRNA Deregulation: Sex-Specific Patterns",
      theme = theme(
        plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "bottom"
      )
    ) &
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 18, face = "bold"),
      legend.text     = element_text(size = 16)
    )
  
  mirna_counts <- all_region_data %>%
    filter(Significant) %>%
    group_by(Region, Sex) %>%
    summarise(
      Upregulated   = sum(log2FC >  fc_threshold),
      Downregulated = sum(log2FC < -fc_threshold),
      Total         = n(),
      .groups = "drop"
    )
  
  complete_grid <- expand.grid(
    Region    = unique(all_region_data$Region),
    Sex       = c("Male", "Female"),
    Direction = c("Downregulated", "Upregulated"),
    stringsAsFactors = FALSE
  )
  
  mirna_counts_long <- mirna_counts %>%
    pivot_longer(cols = c(Upregulated, Downregulated),
                 names_to = "Direction", values_to = "Count") %>%
    right_join(complete_grid, by = c("Region", "Sex", "Direction")) %>%
    mutate(
      Count     = replace_na(Count, 0),
      Direction = factor(Direction, levels = c("Downregulated", "Upregulated")),
      Region    = factor(Region, levels = sort(unique(Region)))
    )
  
  paired_plot <- ggplot(mirna_counts_long, aes(x = Region, y = Count, fill = Sex)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = ifelse(Count > 0, Count, "")),
              position = position_dodge(width = 0.8),
              vjust = -0.4, size = 5, fontface = "bold", color = "black") +
    scale_fill_manual(values = sex_colors) +
    facet_grid(. ~ Direction) +
    labs(
      subtitle = paste0("Thresholds: |log2(FC)| > ", fc_threshold,
                        ", ", pval_label, " < ", pval_threshold),
      x = "Brain Region",
      y = "Number of Deregulated miRNAs (PLX vs Control)",
      fill = "Sex"
    ) +
    ylim(0, 10) +
    theme_minimal() +
    theme(
      text             = element_text(family = "Helvetica", size = 12, color = "black"),
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 18, hjust = 0.5, margin = margin(b = 60)),
      axis.title       = element_text(size = 16, face = "bold"),
      axis.text.x      = element_text(hjust = 0.5, size = 14),
      axis.text.y      = element_text(size = 14),
      strip.text       = element_blank(),
      panel.spacing    = unit(1, "lines"),
      legend.position  = "bottom",
      legend.title     = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 16)
    ) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1), expand = expansion(mult = c(0, 0.05)))
  
  list(
    data          = all_region_data,
    counts        = mirna_counts_long,
    volcano_plot  = combined_regional_plot,
    paired_plot   = paired_plot
  )
}

# Run with raw p-values
res_raw <- make_mirna_regional_plots(
  diff_express_df = diff_express_gender_specific_75per,
  regions    = regions,
  pval_type       = "raw",
  sex_colors      = sex_colors
)

print(res_raw$volcano_plot)
print(res_raw$paired_plot)

# Run with adjusted p-values
res_adj <- make_mirna_regional_plots(
  diff_express_df = diff_express_gender_specific_75per,
  regions    = regions,
  pval_type       = "adj",
  sex_colors      = sex_colors
)

print(res_adj$volcano_plot)
print(res_adj$paired_plot)

#ggsave("Figures and plots/PDFs/combined_regional_volcano_rawp.pdf", res_raw$volcano_plot, width = 14, height = 6)
#ggsave("Figures and plots/PDFs/combined_regional_volcano_adjp.pdf", res_adj$volcano_plot, width = 14, height = 6)

#ggsave("Figures and plots/PNGs/combined_regional_volcano_rawp.png", res_raw$volcano_plot, width = 14, height = 6)
#ggsave("Figures and plots/PNGs/combined_regional_volcano_adjp.png", res_adj$volcano_plot, width = 14, height = 6)

#ggsave("Figures and plots/PDFs/mirna_count_paired_rawp.pdf",  res_raw$paired_plot,  width = 13, height = 8)
#ggsave("Figures and plots/PDFs/mirna_count_paired_adjp.pdf",  res_adj$paired_plot,  width = 13, height = 8)

#ggsave("Figures and plots/PNGs/mirna_count_paired_rawp.png",  res_raw$paired_plot,  width = 13, height = 8)
#ggsave("Figures and plots/PNGs/mirna_count_paired_adjp.png",  res_adj$paired_plot,  width = 13, height = 8)

# ---------------------- PCA Analysis - All Brain Regions ----------------------
### Prepare expression matrix for PCA
pca_data <- data_mmu %>% as.data.frame()
rownames(pca_data) <- pca_data$RNA
pca_data <- pca_data[, -1] %>% t() %>% as.data.frame()

col_vars  <- apply(pca_data, 2, var, na.rm = TRUE)
col_means <- colMeans(pca_data, na.rm = TRUE)
pca_data_filtered <- pca_data[, (col_vars > 0) & (col_means > 0)]

sample_ids <- rownames(pca_data_filtered)
sample_metadata <- data.frame(
  Sample     = sample_ids,
  Region     = sapply(sample_ids, get_region),
  Sample_Num = sapply(sample_ids, get_sample_num)
)

# Assign treatment and sex from region and sample number
ts <- t(mapply(assign_treatment_sex, sample_metadata$Region, sample_metadata$Sample_Num))
sample_metadata$Treatment <- ts[, 1]
sample_metadata$Sex       <- ts[, 2]
sample_metadata_pca <- sample_metadata %>%
  mutate(Group = paste(Treatment, Sex, sep = "_"),
         Region_Treatment = paste(Region, Treatment, sep = "_"))

### Run the PCA and extract variance explained by each PC
pca_result <- prcomp(pca_data_filtered, scale. = TRUE, center = TRUE)
variance_explained <- summary(pca_result)$importance[2, ] * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)

plot_data <- cbind(sample_metadata_pca,
                   data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2]))

pca_theme <- theme_minimal() +
  theme(text         = element_text(family = "Helvetica", size = 12, color = "black"),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 16, face = "bold"),
        plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 16),
        legend.position = "right")

### PCA 1: colored by brain region, shaped by treatment
pca_region <- ggplot(plot_data, aes(PC1, PC2, color = Region, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = c(Control = 16, PLX = 17)) +
  labs(title = "PCA: miRNA Expression Across Brain Regions",
       x = paste0("PC1 (", pc1_var, "% variance)"),
       y = paste0("PC2 (", pc2_var, "% variance)"),
       color = "Brain Region", shape = "Treatment") +
  pca_theme

### PCA 2: colored by treatment, shaped by brain region
pca_treatment <- ggplot(plot_data, aes(PC1, PC2, color = Treatment, shape = Region)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = c("Bregma -2.06" = 16, "Olfactory Bulb" = 17,
                                "Caudate Putamen" = 15, "Cerebellum" = 18)) +
  labs(title = "PCA: Treatment Effect Across Regions",
       x = paste0("PC1 (", pc1_var, "% variance)"),
       y = paste0("PC2 (", pc2_var, "% variance)"),
       color = "Treatment", shape = "Region") +
  pca_theme

### PCA 3: colored by sex, shaped by treatment
pca_sex <- ggplot(plot_data, aes(PC1, PC2, color = Sex, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = sex_colors) +
  scale_shape_manual(values = c(Control = 16, PLX = 17)) +
  labs(title = "PCA: Sex Differences in miRNA Expression",
       x = paste0("PC1 (", pc1_var, "% variance)"),
       y = paste0("PC2 (", pc2_var, "% variance)"),
       color = "Sex", shape = "Treatment") +
  pca_theme

#ggsave("Figures and plots/PCAs/PCA_by_region.pdf", pca_region, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_treatment.pdf", pca_treatment, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_sex.pdf", pca_sex, width = 10, height = 8, dpi = 300)

# ---------------------- Exporting all the tables for supplement ----------------------
OUT_DIR <- "Figures and plots/Generated matrices/For supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXCEL_FONT        <- "Helvetica"
EXCEL_HEADER_FILL <- "#4F81BD"
EXCEL_HEADER_SIZE <- 13
EXCEL_BODY_SIZE   <- 12

make_header_style <- function() {
  createStyle(fontName = EXCEL_FONT, fontSize = EXCEL_HEADER_SIZE,
              fontColour = "white", halign = "center", valign = "center",
              textDecoration = "bold", fgFill = EXCEL_HEADER_FILL,
              border = "Bottom", borderColour = "white")
}

make_body_style <- function() {
  createStyle(fontName = EXCEL_FONT, fontSize = EXCEL_BODY_SIZE,
              halign = "left", valign = "center")
}

apply_sheet_styles <- function(wb, sheet_name, n_rows, n_cols,
                               header_style = make_header_style(),
                               body_style   = make_body_style()) {
  addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0)
    addStyle(wb, sheet_name, body_style,
             rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1:n_cols, widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}

add_index_sheet <- function(wb, sheet_names, descriptions, counts = NULL) {
  summary_df <- data.frame(Sheet_Number = seq_along(sheet_names),
                           Sheet_Name = sheet_names, Description = descriptions)
  if (!is.null(counts)) summary_df$N_entries <- counts
  addWorksheet(wb, "Index", gridLines = TRUE)
  writeData(wb, "Index", summary_df)
  apply_sheet_styles(wb, "Index", nrow(summary_df), ncol(summary_df))
  invisible(summary_df)
}

# Supp Table 2: all miRNAs with raw + adj p-values across regions and sexes
prepare_region_sex_data_full <- function(diff_express_df, region_code, region_name) {
  make_block <- function(sex_code, sex_label) {
    fc_col   <- paste0("fc_Comparison__Treatment_",         region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    rawp_col <- paste0("ttest_rawp_Comparison__Treatment_", region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    adjp_col <- paste0("ttest_adjp_Comparison__Treatment_", region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    
    diff_express_df %>%
      dplyr::select(miRNA = RNA, FoldChange = all_of(fc_col),
                    raw_pvalue = all_of(rawp_col), adj_pvalue = all_of(adjp_col)) %>%
      mutate(
        Region          = region_name,
        Sex             = sex_label,
        log2FC          = log2(FoldChange),
        Significant_raw = (abs(log2FC) > 1 & raw_pvalue < 0.05),
        Significant_adj = (abs(log2FC) > 1 & adj_pvalue < 0.05),
        Direction_raw   = case_when(Significant_raw & log2FC > 0 ~ "Upregulated",
                                    Significant_raw & log2FC < 0 ~ "Downregulated",
                                    TRUE ~ "Not significant"),
        Direction_adj   = case_when(Significant_adj & log2FC > 0 ~ "Upregulated",
                                    Significant_adj & log2FC < 0 ~ "Downregulated",
                                    TRUE ~ "Not significant")
      ) %>%
      dplyr::select(miRNA, Region, Sex, log2FC, FoldChange,
                    raw_pvalue, Direction_raw, Significant_raw,
                    adj_pvalue, Direction_adj, Significant_adj)
  }
  bind_rows(make_block("M", "M"), make_block("F", "F"))
}

all_mirna_full <- map_dfr(names(regions), function(reg_code) {
  prepare_region_sex_data_full(diff_express_gender_specific_75per, reg_code, regions[[reg_code]])
})

### Export
export_table2_all_mirnas <- function() {
  fname <- file.path(OUT_DIR, "Supplementary_Table_2_Detected_miRNAs_DEA.xlsx")
  wb <- createWorkbook()
  hs <- make_header_style()
  bs <- make_body_style()
  
  region_order <- c("Bregma -2.06", "Olfactory Bulb", "Caudate Putamen", "Cerebellum")
  sex_order    <- c("M", "F")
  
  for (reg in region_order) {
    for (sx in sex_order) {
      sheet_name <- paste0(substr(gsub(" ", "_", gsub("-", "", reg)), 1, 18), "_", sx)
      
      df <- all_mirna_full %>% filter(Region == reg, Sex == sx) %>% arrange(raw_pvalue)
      if (nrow(df) == 0) next
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df)
      apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df), hs, bs)
      
      rows_up   <- which(df$Direction_raw == "Upregulated")   + 1
      rows_down <- which(df$Direction_raw == "Downregulated") + 1
      
      if (length(rows_up) > 0)
        addStyle(wb, sheet_name, createStyle(fgFill = "#eebfc3"),
                 rows = rows_up, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
      if (length(rows_down) > 0)
        addStyle(wb, sheet_name, createStyle(fgFill = "#bad2d2"),
                 rows = rows_down, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
    }
  }
  
  summary_raw <- all_mirna_full %>%
    filter(Significant_raw) %>%
    group_by(Region, Sex, Direction_raw) %>%
    summarise(N_miRNAs = n(), .groups = "drop") %>%
    dplyr::rename(Direction = Direction_raw) %>%
    mutate(Threshold = "raw p < 0.05")
  
  summary_adj <- all_mirna_full %>%
    filter(Significant_adj) %>%
    group_by(Region, Sex, Direction_adj) %>%
    summarise(N_miRNAs = n(), .groups = "drop") %>%
    dplyr::rename(Direction = Direction_adj) %>%
    mutate(Threshold = "adj p < 0.05")
  
  summary_df <- bind_rows(summary_raw, summary_adj) %>%
    dplyr::select(Threshold, Region, Sex, Direction, N_miRNAs) %>%
    arrange(Threshold, match(Region, region_order), Sex, Direction)
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_df)
  apply_sheet_styles(wb, "Summary", nrow(summary_df), ncol(summary_df), hs, bs)
  
  summary_idx <- which(names(wb) == "Summary")
  worksheetOrder(wb) <- c(summary_idx, setdiff(seq_along(names(wb)), summary_idx))
  
  saveWorkbook(wb, fname, overwrite = TRUE)
}

#export_table2_all_mirnas()

# ---------------------- miRNA expression across brainstem (PMID 32188880) ----------------------
mirna_set <- c(
  "miR-142a-3p"  = "mmu_miR_142_3p_000464",
  "miR-142a-5p"  = "mmu_miR_142_5p_002248",
  "miR-146a-5p"  = "mmu_miR_146a_000468",
  "miR-200a-3p"  = "mmu_miR_200a_000502",
  "miR-205-5p"   = "mmu_miR_205_000509",
  "miR-223-3p"   = "mmu_miR_223_002295",
  "miR-409-3p"   = "mmu_miR_409_3p_002332",
  "miR-433-5p"   = "mmu_miR_433_5p_001078"
)

cell_type_map <- tibble(
  cell_type = c("Brainstem", "All Neurons", "Motor Neurons", "Astrocytes", "Microglial/Immune"),
  ct_col    = c("mean_bs.ago", "mean_bs.Syn", "mean_bs.ChAT", "mean_bs.GFAP", "mean_bs.Lyz2"),
  det_col   = c("detectable_BS.ago", "detectable_BS.Syn", "detectable_BS.ChAT",
                "detectable_BS.GFAP", "detectable_BS.Lyz2")
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

array_pomper <- read_excel(xlsx_pomper_database)

### Helper function: per-cell-type fold change vs Brainstem
extract_mirna <- function(primer_id) {
  row <- array_pomper %>% filter(`Taqman miRNA primer name` == primer_id)
  if (nrow(row) == 0) stop("Primer not found: ", primer_id)
  
  cell_type_map %>%
    mutate(
      mean_ct  = as.numeric(row[1, ct_col]),
      detected = as.logical(row[1, det_col])
    ) %>%
    mutate(
      ref_ct    = mean_ct[cell_type == "Brainstem"],
      fold      = 2^(ref_ct - mean_ct),
      fold_plot = ifelse(detected, fold, 0),
      nd_label  = ifelse(detected, "", "ND")
    )
}

### Horizontal bar plot
plot_mirna <- function(df, mirna_name) {
  df$cell_type <- factor(df$cell_type, levels = rev(cell_type_map$cell_type))
  
  ggplot(df, aes(x = fold_plot, y = cell_type)) +
    geom_col(fill = bar_color, width = 0.7) +
    geom_text(aes(label = nd_label),
              hjust = -0.3, size = 2.8, color = "grey45", fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = mirna_name,
      subtitle = "Brainstem versus:",
      x        = "Fold change",
      y        = NULL
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title         = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle      = element_text(size = 11, hjust = 0, margin = margin(b = 6)),
      panel.grid         = element_blank(),
      axis.line.x        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x       = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y       = element_blank(),
      axis.text.y        = element_text(size = 10, color = "black"),
      axis.text.x        = element_text(size = 10, color = "black"),
      axis.title.x       = element_text(size = 11, margin = margin(t = 6)),
      plot.margin        = margin(10, 14, 8, 8)
    )
}

### Loop (one plot per miRNA)
walk2(names(mirna_set), unname(mirna_set), function(mirna_name, primer_id) {
  df <- extract_mirna(primer_id)
  cat("\n", mirna_name, "\n", sep = "")
  print(df %>% dplyr::select(cell_type, mean_ct, detected, fold))
  
  p <- plot_mirna(df, mirna_name)
  print(p)
  
  fname <- gsub("[^A-Za-z0-9]+", "_", mirna_name)
  
  #ggsave(file.path(out_dir, paste0(fname, "_brainstem.pdf")), p, width = 3.3, height = 3.0, device = cairo_pdf)
  #ggsave(file.path(out_dir, paste0(fname, "_brainstem.png")), p, width = 3.3, height = 3.0, dpi = 300, bg = "white")
})

# ---------------------- miRNA expression (overlap 3 datasets - PMID: 37217597; 40161726) ----------------------
### Detection threshold for external sorted-cell datasets (rpmm)
RPMM_DETECTION_THRESHOLD <- 5

rpmm_external <- read.delim(
  "Input files/miRNA_complete_quantification_rpmm_norm.csv",
  header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)

annot_external <- read.delim(
  "Input files/annotation.csv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

### Define sample groups from annotation
engel_young_ids <- annot_external %>%
  dplyr::filter(Comparison_Group == "12W-microglia") %>% dplyr::pull(ID)
engel_aged_ids  <- annot_external %>%
  dplyr::filter(Comparison_Group == "84W-microglia") %>% dplyr::pull(ID)

walsh_cd45_F_ids <- annot_external %>%
  dplyr::filter(Comparison_Group == "8W-F-CD45") %>% dplyr::pull(ID)
walsh_cd45_M_ids <- annot_external %>%
  dplyr::filter(Comparison_Group == "8W-M-CD45") %>% dplyr::pull(ID)

rpmm_mir <- rpmm_external %>% dplyr::filter(grepl("^mmu-miR-", miRNA))

rpmm_mir$engel_young_mean <- rowMeans(rpmm_mir[, engel_young_ids, drop = FALSE], na.rm = TRUE)
rpmm_mir$engel_aged_mean  <- rowMeans(rpmm_mir[, engel_aged_ids,  drop = FALSE], na.rm = TRUE)

engel_detected <- rpmm_mir %>%
  dplyr::filter(engel_young_mean >= RPMM_DETECTION_THRESHOLD |
                  engel_aged_mean  >= RPMM_DETECTION_THRESHOLD) %>%
  dplyr::pull(miRNA) %>%
  strip_mmu() %>%
  sort() %>%
  unique()

### Walsh: mean rpmm per sex in CD45+ fraction, detected if >= threshold in either sex
rpmm_mir$walsh_cd45_F_mean <- rowMeans(rpmm_mir[, walsh_cd45_F_ids, drop = FALSE], na.rm = TRUE)
rpmm_mir$walsh_cd45_M_mean <- rowMeans(rpmm_mir[, walsh_cd45_M_ids, drop = FALSE], na.rm = TRUE)

walsh_detected <- rpmm_mir %>%
  dplyr::filter(walsh_cd45_F_mean >= RPMM_DETECTION_THRESHOLD |
                  walsh_cd45_M_mean >= RPMM_DETECTION_THRESHOLD) %>%
  dplyr::pull(miRNA) %>%
  strip_mmu() %>%
  sort() %>%
  unique()

region_codes_depleted <- c(S  = "Bregma -2.06",
                           OB = "Olfactory Bulb",
                           CP = "Caudate Putamen",
                           CE = "Cerebellum")

get_depleted_set <- function(df, region_codes,
                             fc_threshold = 1, p_threshold = 0.05) {
  hits <- character(0)
  for (rc in names(region_codes)) for (sx in c("M", "F")) {
    fc_col <- paste0("fc_Comparison__Treatment_", rc, "_", sx, "_vs_Control_", rc, "_", sx)
    p_col  <- paste0("ttest_rawp_Comparison__Treatment_", rc, "_", sx, "_vs_Control_", rc, "_", sx)
    if (!all(c(fc_col, p_col) %in% names(df))) next
    fc <- df[[fc_col]]; p <- df[[p_col]]
    l2 <- ifelse(fc > 0, log2(fc), NA_real_)
    down <- !is.na(l2) & !is.na(p) & p < p_threshold & l2 < -fc_threshold
    hits <- c(hits, strip_mmu(df$RNA[down]))
  }
  sort(unique(hits))
}

depleted_this_study <- get_depleted_set(diff_express_gender_specific_75per, region_codes_depleted)

cat("Set sizes —",
    "Depleted:", length(depleted_this_study),
    "| Engel detected:", length(engel_detected),
    "| Walsh detected:", length(walsh_detected), "\n")

### Compute set intersections
A <- depleted_this_study
B <- engel_detected
C <- walsh_detected

triple  <- Reduce(intersect, list(A, B, C))
AB_only <- setdiff(intersect(A, B), C)
AC_only <- setdiff(intersect(A, C), B)
BC_only <- setdiff(intersect(B, C), A)

cat("Triple overlap:", paste(triple, collapse = ", "), "\n")
if (length(AB_only) > 0) cat("Depleted ∩ Engel only:", paste(AB_only, collapse = ", "), "\n")
if (length(AC_only) > 0) cat("Depleted ∩ Walsh only:", paste(AC_only, collapse = ", "), "\n")
cat("Engel ∩ Walsh only:", length(BC_only), "miRNAs\n")

### Euler plot
euler_input <- c(
  "Depleted"       = length(setdiff(setdiff(A, B), C)),
  "Engel"          = length(setdiff(setdiff(B, A), C)),
  "Walsh"          = length(setdiff(setdiff(C, A), B)),
  "Depleted&Engel" = length(AB_only),
  "Depleted&Walsh" = length(AC_only),
  "Engel&Walsh"    = length(BC_only),
  "Depleted&Engel&Walsh" = length(triple)
)

fit <- euler(euler_input)

label_fn <- function(v, max_names = 4) {
  if (length(v) == 0) return("")
  if (length(v) > max_names) return(as.character(length(v)))
  paste(v, collapse = "\n")
}

region_counts <- c(
  label_fn(setdiff(setdiff(A, B), C)),
  label_fn(setdiff(setdiff(B, A), C)),
  label_fn(setdiff(setdiff(C, A), B)),
  label_fn(AB_only),
  label_fn(AC_only),
  label_fn(BC_only),
  label_fn(triple)
)

p_euler <- plot(fit,
                fills = list(fill = c("#ad75a1", "#c7a469", "#6ca2a7"), alpha = 0.45),
                edges = list(col = "grey40", lwd = 0.5),
                labels = list(
                  labels = c(
                    "Downregulated\nafter PLX\n(this study)",
                    "Detected in\nFACS-sorted microglia\n(Engel et al.)",
                    "Detected in\nCD45+ sorted cells\n(Walsh et al.)"
                  ),
                  fontfamily = "Helvetica",
                  fontsize = 10
                ),
                quantities = list(
                  labels = region_counts,
                  fontfamily = "Helvetica",
                  fontsize = 8,
                  fontface = "italic"
                ))

print(p_euler)

# ---------------------- miRNA biological variability ----------------------
# Expression matrix
mirna_mat <- data_mmu %>%
  as.data.frame() %>%
  { rownames(.) <- .$RNA; . } %>%
  dplyr::select(-RNA) %>%
  as.matrix()

mirna_mat <- log2(mirna_mat + 1)

expressed_mirnas <- rowMeans(mirna_mat) > log2(1 + 1) &
  matrixStats::rowVars(mirna_mat) > 0
mirna_mat <- mirna_mat[expressed_mirnas, , drop = FALSE]
cat(sprintf("Variance partition: %d miRNAs x %d samples retained\n",
            nrow(mirna_mat), ncol(mirna_mat)))

sample_info <- data.frame(
  Sample     = colnames(mirna_mat),
  Region     = sapply(colnames(mirna_mat), get_region),
  Sample_Num = sapply(colnames(mirna_mat), get_sample_num),
  stringsAsFactors = FALSE
)

group_assignments      <- t(mapply(assign_treatment_sex,
                                   sample_info$Region,
                                   sample_info$Sample_Num))
sample_info$Treatment  <- group_assignments[, 1]
sample_info$Sex        <- group_assignments[, 2]
sample_info$Animal     <- factor(paste0("M", sample_info$Sample_Num))
sample_info$Region     <- factor(sample_info$Region)
sample_info$Treatment  <- factor(sample_info$Treatment, levels = c("Control", "PLX"))
sample_info$Sex        <- factor(sample_info$Sex,       levels = c("Female", "Male"))

rownames(sample_info) <- sample_info$Sample
sample_info <- sample_info[colnames(mirna_mat), ]
stopifnot(identical(rownames(sample_info), colnames(mirna_mat)))

cat("\nAnimal x Region table (should be ~1 per cell):\n")
print(table(sample_info$Animal, sample_info$Region))
cat("\nAnimal -> Sex/Treatment (must be 1:1):\n")
print(unique(sample_info[, c("Animal", "Sex", "Treatment")]) %>% arrange(Animal))

VPART_COLORS <- c(
  Region    = "#8b3a7a",
  Animal    = "#538f8f",
  Treatment = "#b07d2a",
  Sex       = "#a42e30",
  Residuals = "#d3d3d3"
)

VPART_THEME <- theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title        = element_text(size = 14, face = "bold",   hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.text         = element_text(size = 12, color = "black"),
    axis.title.y      = element_text(size = 13, face = "bold"),
    panel.border      = element_blank(),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    legend.position   = "none"
  )

### Global partition: all regions
formula_all  <- ~ (1|Region) + (1|Sex) + (1|Treatment) + (1|Animal)
vpart_all    <- fitExtractVarPartModel(mirna_mat, formula_all, sample_info)
vpart_all    <- sortCols(vpart_all)

cat("\n=== Global variance partition (median % across miRNAs) ===\n")
print(round(sort(colMeans(vpart_all) * 100, decreasing = TRUE), 2))

p_vpart_all <- as.data.frame(unclass(vpart_all), check.names = FALSE) %>%
  pivot_longer(cols = everything(),
               names_to  = "Component",
               values_to = "Variance") %>%
  mutate(
    Component = factor(Component, levels = colnames(vpart_all)),
    Variance  = Variance * 100
  ) %>%
  ggplot(aes(x = Component, y = Variance, fill = Component)) +
  geom_violin(scale = "width", trim = FALSE, colour = "white") +
  geom_boxplot(width = 0.07, outlier.shape = NA,
               fill = "white", colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = VPART_COLORS) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "miRNA variance partition - all brain regions",
       y     = "Variance explained (%)") +
  VPART_THEME +
  theme(axis.title.x = element_blank())

print(p_vpart_all)

#ggsave("Figures and plots/PNGs/VarPart_miRNA_global.png", p_vpart_all, width = 8, height = 6, dpi = 300)
#ggsave("Figures and plots/PDFs/VarPart_miRNA_global.pdf", p_vpart_all, width = 8, height = 6)

# ---------------------- Combined FC plot: study depletion effect x Pomper cell-type enrichment ----------------------
mirna_order <- names(mirna_set)

### Build Pomper matrix
pomper_matrix_df <- map_dfr(names(mirna_set), function(mirna_name) {
  primer_id <- mirna_set[[mirna_name]]
  df <- tryCatch(extract_mirna(primer_id), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df %>% mutate(miRNA = mirna_name)
})

pomper_matrix <- pomper_matrix_df %>%
  mutate(log2FC = ifelse(detected, log2(fold), NA_real_)) %>%
  dplyr::select(miRNA, cell_type, log2FC) %>%
  pivot_wider(names_from = cell_type, values_from = log2FC) %>%
  column_to_rownames("miRNA")

pomper_matrix <- pomper_matrix[mirna_order, cell_type_map$cell_type]
pomper_matrix <- as.matrix(pomper_matrix)

region_sex_cols <- expand.grid(
  Region = c("Bregma -2.06", "Olfactory Bulb", "Caudate Putamen", "Cerebellum"),
  Sex = c("M", "F"),
  stringsAsFactors = FALSE
) %>% mutate(col_label = paste0(Region, "_", Sex))

study_matrix_df <- all_mirna_full %>%
  mutate(miRNA_clean = strip_mmu(miRNA)) %>%
  filter(miRNA_clean %in% mirna_order) %>%
  mutate(col_label = paste0(Region, "_", Sex))

study_matrix <- study_matrix_df %>%
  dplyr::select(miRNA_clean, col_label, log2FC) %>%
  pivot_wider(names_from = col_label, values_from = log2FC) %>%
  column_to_rownames("miRNA_clean")

study_matrix <- study_matrix[mirna_order, region_sex_cols$col_label]
study_matrix <- as.matrix(study_matrix)

sig_matrix_df <- study_matrix_df %>%
  dplyr::select(miRNA_clean, col_label, Significant_raw) %>%
  pivot_wider(names_from = col_label, values_from = Significant_raw)

sig_matrix <- sig_matrix_df %>%
  column_to_rownames("miRNA_clean")
sig_matrix <- sig_matrix[mirna_order, region_sex_cols$col_label]
sig_matrix <- as.matrix(sig_matrix)
sig_matrix[is.na(sig_matrix)] <- FALSE

study_col_ha <- HeatmapAnnotation(
  Region = region_sex_cols$Region,
  Sex    = ifelse(region_sex_cols$Sex == "M", "Male", "Female"),
  col = list(
    Region = region_colors,
    Sex    = sex_colors
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Region = list(title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 10)),
    Sex    = list(title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 10))
  )
)

pomper_col_fun <- colorRamp2(
  c(0, max(pomper_matrix, na.rm = TRUE) / 2, max(pomper_matrix, na.rm = TRUE)),
  c("white", "#6ca2a7", "#2E5560")
)

study_max <- max(abs(study_matrix), na.rm = TRUE)
study_col_fun <- colorRamp2(
  c(-study_max, 0, study_max),
  c("#254e7a", "white", "#8b3a7a")
)

### Build the two heatmaps
ht_pomper <- Heatmap(
  pomper_matrix,
  name = "log2FC\n(vs Brainstem)",
  col = pomper_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontface = "italic", fontsize = 10),
  column_names_rot = 40,
  column_names_gp = gpar(fontsize = 10),
  column_title = "Cell-type enrichment\n(Pomper et al.)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  na_col = "grey90",
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
)

ht_study <- Heatmap(
  study_matrix,
  name = "log2FC\n(PLX vs Control)",
  col = study_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = study_col_ha,
  column_title = "Depletion effect (this study)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  na_col = "grey90",
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_matrix[i, j]) && sig_matrix[i, j]) {
      grid.text("*", x, y, gp = gpar(fontsize = 16, fontface = "bold", col = "black"))
    }
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                              labels_gp = gpar(fontsize = 8))
)

combined_ht <- ht_pomper + ht_study

#pdf("Figures and plots/PDFs/miRNA_celltype_depletion_heatmap.pdf", width = 6, height = 8)
draw(combined_ht,
     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
     merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
#dev.off()

# ---------------------- miR-223-3p in Engel FACS-sorted microglia (young vs aged) ----------------------
engel_samples <- annot_external %>%
  dplyr::filter(Dataset == "FACS_sorted_microglia") %>%
  dplyr::select(ID, Comparison_Group, Age)

### Extract miR-223-3p expression for Engel samples
mir223_vals <- rpmm_external %>%
  dplyr::filter(miRNA == "mmu-miR-223-3p") %>%
  dplyr::select(all_of(engel_samples$ID)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "ID", values_to = "rpmm") %>%
  dplyr::left_join(engel_samples, by = "ID") %>%
  dplyr::mutate(
    Age_Group = factor(
      ifelse(Comparison_Group == "12W-microglia", "Young (12W)", "Aged (84W)"),
      levels = c("Young (12W)", "Aged (84W)")
    )
  )

young_vals <- mir223_vals %>% dplyr::filter(Age_Group == "Young (12W)") %>% dplyr::pull(rpmm)
aged_vals  <- mir223_vals %>% dplyr::filter(Age_Group == "Aged (84W)")  %>% dplyr::pull(rpmm)

fc_val    <- mean(aged_vals) / mean(young_vals)
log2fc    <- log2(fc_val)
ttest_res <- t.test(aged_vals, young_vals, var.equal = FALSE)

age_colors <- c("Young (12W)" = "#7ea172", "Aged (84W)" = "#e7a977")

p_mir223_engel <- ggplot(mir223_vals, aes(x = Age_Group, y = rpmm, fill = Age_Group)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, color = NA) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2) +
  geom_jitter(color = "black", width = 0.1, alpha = 0.6, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = age_colors) +
  ggsignif::geom_signif(
    comparisons = list(c("Young (12W)", "Aged (84W)")),
    test = "t.test",
    test.args = list(var.equal = FALSE),
    map_signif_level = function(p) sprintf("p = %.3f", p),
    tip_length = 0.01,
    textsize = 4.5) +
  annotate("text",
           x = 1.5, y = max(aged_vals) * 0.85,
           label = sprintf("FC = %.1f\nlog2FC = %.2f", fc_val, log2fc),
           size = 4, hjust = 0.5, color = "grey20", fontface = "italic",
           family = "Helvetica") +
  labs(
    title = "miR-223-3p in FACS-sorted microglia (Engel et al. 2025)",
    x = "", y = "rpmm normalized counts") +
  ylim(0, 12) +
  theme_minimal() +
  theme(
    text             = element_text(family = "Helvetica"),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y     = element_text(size = 14, face = "bold"),
    axis.text.x      = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    legend.position  = "none",
    panel.grid.major.x = element_blank())

print(p_mir223_engel)

#ggsave("Figures and plots/PNGs/miR223_engel_FACS_barplot.png", p_mir223_engel, width = 7, height = 8)
#ggsave("Figures and plots/PDFs/miR223_engel_FACS_barplot.pdf", p_mir223_engel, width = 7, height = 8)

# ---------------------- Pomper cell-type enrichment for all deregulated miRNAs ----------------------
### Collect all deregulated miRNAs (any region, any sex, any direction)
get_all_dereg_mirnas <- function(df, region_codes,
                                 fc_threshold = 1, p_threshold = 0.05) {
  hits <- list()
  for (rc in names(region_codes)) {
    for (sx in c("M", "F")) {
      fc_col <- paste0("fc_Comparison__Treatment_", rc, "_", sx,
                       "_vs_Control_", rc, "_", sx)
      p_col  <- paste0("ttest_rawp_Comparison__Treatment_", rc, "_", sx,
                       "_vs_Control_", rc, "_", sx)
      if (!all(c(fc_col, p_col) %in% names(df))) next
      
      fc <- df[[fc_col]]; p <- df[[p_col]]
      l2 <- ifelse(fc > 0, log2(fc), NA_real_)
      sig <- !is.na(l2) & !is.na(p) & p < p_threshold & abs(l2) > fc_threshold
      
      if (any(sig)) {
        hits[[length(hits) + 1]] <- data.frame(
          miRNA      = strip_mmu(df$RNA[sig]),
          Region     = region_codes[[rc]],
          Sex        = sx,
          log2FC     = l2[sig],
          raw_pvalue = p[sig],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(hits)
}

all_dereg_detail <- get_all_dereg_mirnas(diff_express_gender_specific_75per, region_codes_depleted)
unique_dereg_mirnas <- sort(unique(all_dereg_detail$miRNA))

### Match to Pomper primer IDs using the existing mirna_set mapping
mirna_to_primer <- setNames(rep(NA_character_, length(unique_dereg_mirnas)), unique_dereg_mirnas)
for (m in unique_dereg_mirnas) {
  if (m %in% names(mirna_set)) mirna_to_primer[m] <- mirna_set[[m]]
}

matched  <- names(mirna_to_primer[!is.na(mirna_to_primer)])
no_match <- names(mirna_to_primer[is.na(mirna_to_primer)])
if (length(no_match) > 0) {
  cat("Not in Pomper array:", paste(no_match, collapse = ", "), "\n")
}

pomper_rows <- lapply(matched, function(m) {
  df <- tryCatch(extract_mirna(mirna_to_primer[m]), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  df %>% dplyr::mutate(miRNA = m)
})
pomper_table <- bind_rows(pomper_rows)

### Wide format: one row per miRNA, one column per cell type# "ND" for not detected
pomper_wide <- pomper_table %>%
  dplyr::mutate(
    fold_display = ifelse(detected, as.character(round(fold, 2)), "ND")
  ) %>%
  dplyr::select(miRNA, cell_type, fold_display) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = fold_display)

### Deregulation summary across regions/sexes
dereg_summary <- all_dereg_detail %>%
  dplyr::group_by(miRNA) %>%
  dplyr::summarise(
    n_region_sex = n(),
    Regions      = paste(unique(Region), collapse = "; "),
    Sexes        = paste(sort(unique(Sex)), collapse = ", "),
    Direction    = ifelse(all(log2FC > 0), "Up",
                          ifelse(all(log2FC < 0), "Down", "Mixed")),
    .groups = "drop"
  )

# Combine
ct_cols <- cell_type_map$cell_type
final_table <- dereg_summary %>%
  dplyr::left_join(pomper_wide, by = "miRNA") %>%
  dplyr::mutate(across(all_of(ct_cols), ~ replace_na(.x, "NA"))) %>%
  dplyr::arrange(Direction, miRNA)

export_table3_pomper_enrichment <- function(df) {
  fname <- file.path(OUT_DIR, "Supplementary_Table_3_Pomper_CellType_Enrichment.xlsx")
  wb <- createWorkbook()
  hs <- make_header_style()
  bs <- make_body_style()
  
  addWorksheet(wb, "Cell-Type Enrichment")
  writeData(wb, "Cell-Type Enrichment", df)
  apply_sheet_styles(wb, "Cell-Type Enrichment", nrow(df), ncol(df), hs, bs)
  
  mic_rows <- which(!df$`Microglial/Immune` %in% c("ND", "NA")) + 1
  if (length(mic_rows) > 0) {
    addStyle(wb, "Cell-Type Enrichment",
             createStyle(fgFill = "#e8d5e0"),
             rows = mic_rows, cols = 1:ncol(df),
             gridExpand = TRUE, stack = TRUE)
  }
  
  nd_style <- createStyle(fontColour = "#999999")
  for (col_idx in which(names(df) %in% ct_cols)) {
    nd_rows <- which(df[[col_idx]] %in% c("ND", "NA")) + 1
    if (length(nd_rows) > 0) {
      addStyle(wb, "Cell-Type Enrichment", nd_style,
               rows = nd_rows, cols = col_idx,
               gridExpand = TRUE, stack = TRUE)
    }
  }
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

#export_table3_pomper_enrichment(final_table)