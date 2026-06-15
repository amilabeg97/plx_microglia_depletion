# Setting the working directory and importing libraries
#setwd("//")

suppressMessages({
  library(tidyverse)
  library(ggpubr)
  library(pheatmap)
  library(matrixStats)
  library(readxl)
  library(data.table)
  library(svglite)
  library(factoextra)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(VennDiagram)
  library(ggvenn)
  library(UpSetR)
  library(ComplexUpset)
  library(igraph)
  library(ggraph)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(patchwork)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(openxlsx)
  library(purrr)
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

regions <- list(
  S = "Bregma -2.06",
  OB = "Olfactory Bulb",
  CP = "Caudate Putamen",
  CE = "Cerebellum"
)

# ---------------------- Load data miRNA ----------------------
### Depletion study data load
data <- read.table("Input files/sncRNA_complete_quantification_rpmm_norm.csv", 
                   header = TRUE, sep = "\t")

colnames(data) <- c(
  "RNA", "S15", "S4", "S5", "S6", "S7", "S8", "S19", "S20", "S16",  
  "S9", "S10", "S11", "S17", "S12", "S13", "S14", "S18", "S1", "S2",    
  "S3", "OB17", "OB4", "OB5", "OB8", "OB9", "OB10", "OB18", "OB11", "OB12",  
  "OB13", "OB19", "OB14", "OB15", "OB16", "OB20", "OB1", "OB2", "OB3", "CP16",
  "CP4", "CP5", "CP6", "CP7", "CP8", "CP20", "CP17", "CP10", "CP11", "CP12",  
  "CP13", "CP18", "CP14", "CP15", "CP19", "CP1", "CP2", "CP3", "CE15", "CE4",
  "CE5", "CE6", "CE7", "CE8", "CE19", "CE20", "CE16", "CE9", "CE10", "CE11", 
  "CE12", "CE13", "CE17", "CE14", "CE18", "CE1", "CE2", "CE3")

data_mmu <- data[grep("mmu-miR-", data$RNA), ]

colnames(data_mmu) <- c(
  "RNA", "S15", "S4", "S5", "S6", "S7", "S8", "S19", "S20", "S16",  
  "S9", "S10", "S11", "S17", "S12", "S13", "S14", "S18", "S1", "S2",    
  "S3", "OB17", "OB4", "OB5", "OB8", "OB9", "OB10", "OB18", "OB11", "OB12",  
  "OB13", "OB19", "OB14", "OB15", "OB16", "OB20", "OB1", "OB2", "OB3", "CP16",
  "CP4", "CP5", "CP6", "CP7", "CP8", "CP20", "CP17", "CP10", "CP11", "CP12",  
  "CP13", "CP18", "CP14", "CP15", "CP19", "CP1", "CP2", "CP3", "CE15", "CE4",
  "CE5", "CE6", "CE7", "CE8", "CE19", "CE20", "CE16", "CE9", "CE10", "CE11", 
  "CE12", "CE13", "CE17", "CE14", "CE18", "CE1", "CE2", "CE3")

miRNAs_of_interest <- c("mmu-miR-146a-5p", "mmu-miR-223-3p")

diff_express_gender_specific_75per <- read.table(
  "snc-pipeline/plx_microglia_gender_specific_75%/diff_exp/miRNA/rpmm_norm_log2.detection_rate_75p_per_group/diff_exp_gender_specific_75per.csv", 
  header = TRUE, sep = "\t"
)

### Pomper database
xlsx_pomper_database <- "Input files/LoessM_norm_miRNA_array_WEBSITE_final.xlsx"
out_dir   <- "Figures and plots/miRNA_brainstem_plots"
bar_color <- "#2E5560"

### Load the diff_exp table from two external datasets
diff_exp_external <- read.table(
  "Input files/results_plx_walsh_microglia/diff_exp/miRNA/rpmm_norm_log2.detection_rate_75p_per_group/diff_exp.csv",
  header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE
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

# Helper functions
get_region <- function(sample_id) {
  if (grepl("^S[0-9]", sample_id)) return("Bregma -2.06")
  if (grepl("^OB", sample_id)) return("Olfactory Bulb")
  if (grepl("^CP", sample_id)) return("Caudate Putamen")
  if (grepl("^CE", sample_id)) return("Cerebellum")
  return("Unknown")
}

assign_treatment_sex <- function(region, sample_num) {
  treatment <- if (sample_num <= 12) "PLX" else "Control"
  sex <- if (sample_num <= 6 || (sample_num >= 13 && sample_num <= 16)) "Female" else "Male"
  
  return(c(treatment, sex))
}

### Sanity check that the samples are correctly labeled
cat("Sample 1-6 (should be PLX Female):\n")
for (i in 1:6) {
  result <- assign_treatment_sex("Any Region", i)
  cat(sprintf("  Sample %d: %s %s\n", i, result[1], result[2]))
}

cat("\nSample 7-12 (should be PLX Male):\n")
for (i in 7:12) {
  result <- assign_treatment_sex("Any Region", i)
  cat(sprintf("  Sample %d: %s %s\n", i, result[1], result[2]))
}

cat("\nSample 13-16 (should be Control Female):\n")
for (i in 13:16) {
  result <- assign_treatment_sex("Any Region", i)
  cat(sprintf("  Sample %d: %s %s\n", i, result[1], result[2]))
}

cat("\nSample 17-20 (should be Control Male):\n")
for (i in 17:20) {
  result <- assign_treatment_sex("Any Region", i)
  cat(sprintf("  Sample %d: %s %s\n", i, result[1], result[2]))
}

get_sample_num <- function(sample_id) {
  num <- gsub("^[A-Z]+", "", sample_id)
  return(as.numeric(num))
}

all_samples <- colnames(corr_data)
sample_regions <- sapply(all_samples, get_region)

#cat("Sample distribution across regions:\n")
#print(table(sample_regions))

corr_matrix_full <- as.matrix(corr_data)

# Color scheme
col_fun <- colorRamp2(c(0.85, 0.90, 0.95, 1.0), 
                      c("#9ECAE1", "#4292C6", "#08519C", "#08306B"))


### Get metadata for all samples
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

### Save the heatmap
#pdf("Figures and plots/PDFs/All_Regions_Single_correlation_heatmap.pdf", width = 11, height = 13)
#draw(ht_all, heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()

# ---------------------- QC Number of mapped lncRNAs ----------------------
### Read the QC data
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
qc_data$Treatment <- assignments[,1]
qc_data$Sex <- assignments[,2]
qc_data$Group <- paste(qc_data$Treatment, qc_data$Sex, sep = "_")

# ncRNA mapping variables (one column per ncRNA type)
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
  
  # Place each existing sample at its own position 
  for (i in 1:nrow(region_qc)) {
    sample_num <- region_qc$Sample_Number[i]
    complete_data$Treatment[sample_num] <- region_qc$Treatment[i]
    complete_data$Sex[sample_num] <- region_qc$Sex[i]
    complete_matrix[, sample_num] <- ncRNA_matrix[region_idx[i], ]
  }
  
  return(list(metadata = complete_data, matrix = complete_matrix))
}

# Shared y-axis max across all regions
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

### Create plot by region - vertical layout
#pdf("Figures and plots/PDFs/ncRNA_mapping_by_region_improved.pdf", width = 8, height = 11)

# Four region panels stacked vertically + 1 legend at the bottom
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
  
  # Bar width from a dry-run barplot, needed to size the annotation rectangles below
  test_bp <- barplot(region_data[1,], plot = FALSE, space = 0.35)
  bar_width <- diff(test_bp)[1] / 1.35
  
  # Open an empty plotting region first so gridlines can be drawn behind the bars
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
  
  par(xpd = TRUE) # allow drawing into the margin below the plot
  
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

# Three legends in the bottom panel: ncRNA type, sex, treatment
legend(
  x = "right",
  inset = c(0, 0),
  legend = c("lncRNA","miscRNA","rRNA","scaRNA",
             "snoRNA","snRNA","miRNA","tRNA"),
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
      stat_compare_means(comparisons = comparisons_filtered, method = "t.test",
                         label = "p.signif", tip.length = 0.01) +
      theme_minimal() +
      ylim(0, 25) +
      labs(title = gene, x = "", y = "rpmm normalized counts") +
      theme(
        text = element_text(family = "Helvetica", size = 14, color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  })

plots_with_stats_small_rna

# Save the plots
#pdf("Figures and plots/miRNA_barplots_small_rna.pdf", width = 6, height = 7)
# for(i in seq_along(plots_with_stats_small_rna)) {
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
  
  ### Paired bar plot
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
    ylim(0,10) +
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

### Run the functions
# Raw p-value
res_raw <- make_mirna_regional_plots(
  diff_express_df = diff_express_gender_specific_75per,
  regions    = regions,
  pval_type       = "raw",
  sex_colors      = sex_colors
)

print(res_raw$volcano_plot)
print(res_raw$paired_plot)

# Adjusted p-value 
res_adj <- make_mirna_regional_plots(
  diff_express_df = diff_express_gender_specific_75per,
  regions    = regions,
  pval_type       = "adj",
  sex_colors      = sex_colors
)

print(res_adj$volcano_plot)
print(res_adj$paired_plot)


### Save the volano plots
#ggsave("Figures and plots/PDFs/combined_regional_volcano_rawp.pdf", res_raw$volcano_plot, width = 14, height = 6)
#ggsave("Figures and plots/PDFs/combined_regional_volcano_adjp.pdf", res_adj$volcano_plot, width = 14, height = 6)

#ggsave("Figures and plots/PNGs/combined_regional_volcano_rawp.png", res_raw$volcano_plot, width = 14, height = 6)
#ggsave("Figures and plots/PNGs/combined_regional_volcano_adjp.png", res_adj$volcano_plot, width = 14, height = 6)

### Save the miRNA count summary 
#ggsave("Figures and plots/PDFs/mirna_count_paired_rawp.pdf",  res_raw$paired_plot,  width = 13, height = 8)
#ggsave("Figures and plots/PDFs/mirna_count_paired_adjp.pdf",  res_adj$paired_plot,  width = 13, height = 8)

#ggsave("Figures and plots/PNGs/mirna_count_paired_rawp.png",  res_raw$paired_plot,  width = 13, height = 8)
#ggsave("Figures and plots/PNGs/mirna_count_paired_adjp.png",  res_adj$paired_plot,  width = 13, height = 8)


# ---------------------- PCA Analysis - All Brain Regions ----------------------
### Prepare expression matrix for PCA
pca_data <- data_mmu %>%
  as.data.frame()

rownames(pca_data) <- pca_data$RNA
pca_data <- pca_data[, -1]
pca_data <- t(pca_data) %>% as.data.frame()

# Filter out miRNAs with zero variance or zero mean
col_vars <- apply(pca_data, 2, var, na.rm = TRUE)
col_means <- colMeans(pca_data, na.rm = TRUE)
cols_to_keep <- (col_vars > 0) & (col_means > 0)
pca_data_filtered <- pca_data[, cols_to_keep]

# Build sample metadata
sample_ids <- rownames(pca_data_filtered)
sample_metadata <- data.frame(
  Sample     = sample_ids,
  Region     = sapply(sample_ids, get_region),
  Sample_Num = sapply(sample_ids, get_sample_num),
  stringsAsFactors = FALSE
)

# Assign treatment and sex from region and sample number
ts <- t(mapply(assign_treatment_sex, sample_metadata$Region, sample_metadata$Sample_Num))
sample_metadata$Treatment <- ts[, 1]
sample_metadata$Sex        <- ts[, 2]

sample_metadata_pca <- sample_metadata %>%
  mutate(Group = paste(Treatment, Sex, sep = "_"),
         Region_Treatment = paste(Region, Treatment, sep = "_"))

### Run the PCA and extract variance explained by each PC
pca_result <- prcomp(pca_data_filtered, scale. = TRUE, center = TRUE)
variance_explained <- summary(pca_result)$importance[2,] * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)

pca_coords <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2])
plot_data <- cbind(sample_metadata_pca, pca_coords)

### PCA 1: colored by brain region, shaped by treatment
pca_region <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Region, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = region_colors) +
  scale_shape_manual(values = c("Control" = 16, "PLX" = 17)) +
  theme_minimal() +
  labs(
    title = "PCA: miRNA Expression Across Brain Regions",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Brain Region", shape = "Treatment"
  ) +
  theme(
    text = element_text(family = "Helvetica", size = 12, color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

### PCA 2: colored by treatment, shaped by brain region
pca_treatment <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Treatment, shape = Region)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = c("Bregma -2.06" = 16, "Olfactory Bulb" = 17,
                                "Caudate Putamen" = 15, "Cerebellum" = 18)) +
  theme_minimal() +
  labs(
    title = "PCA: Treatment Effect Across Regions",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Treatment", shape = "Region"
  ) +
  theme(
    text = element_text(family = "Helvetica", size = 12, color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

### PCA 3: colored by sex, shaped by treatment
pca_sex <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Sex, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = sex_colors) +
  scale_shape_manual(values = c("Control" = 16, "PLX" = 17)) +
  theme_minimal() +
  labs(
    title = "PCA: Sex Differences in miRNA Expression",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Sex", shape = "Treatment"
  ) +
  theme(
    text = element_text(family = "Helvetica", size = 12, color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

#ggsave("Figures and plots/PCAs/PCA_by_region.png", pca_region, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_treatment.png", pca_treatment, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_sex.png", pca_sex, width = 10, height = 8, dpi = 300)

#ggsave("Figures and plots/PCAs/PCA_by_region.pdf", pca_region, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_treatment.pdf", pca_treatment, width = 10, height = 8, dpi = 300)
#ggsave("Figures and plots/PCAs/PCA_by_sex.pdf", pca_sex, width = 10, height = 8, dpi = 300)

# ---------------------- Exporting all the tables for supplement ----------------------
### Output directory
OUT_DIR <- "Figures and plots/Generated matrices/For supplement"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared style helpers
EXCEL_FONT        <- "Helvetica"
EXCEL_HEADER_FILL <- "#4F81BD"
EXCEL_HEADER_SIZE <- 13
EXCEL_BODY_SIZE   <- 12

make_header_style <- function() {
  createStyle(
    fontName       = EXCEL_FONT,
    fontSize       = EXCEL_HEADER_SIZE,
    fontColour     = "white",
    halign         = "center",
    valign         = "center",
    textDecoration = "bold",
    fgFill         = EXCEL_HEADER_FILL,
    border         = "Bottom",
    borderColour   = "white"
  )
}

make_body_style <- function() {
  createStyle(
    fontName = EXCEL_FONT,
    fontSize = EXCEL_BODY_SIZE,
    halign   = "left",
    valign   = "center"
  )
}

apply_sheet_styles <- function(wb, sheet_name, n_rows, n_cols,
                               header_style = make_header_style(),
                               body_style   = make_body_style()) {
  addStyle(wb, sheet_name, header_style,
           rows = 1, cols = 1:n_cols, gridExpand = TRUE)
  if (n_rows > 0)
    addStyle(wb, sheet_name, body_style,
             rows = 2:(n_rows + 1), cols = 1:n_cols, gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1:n_cols, widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}

add_index_sheet <- function(wb, sheet_names, descriptions, counts = NULL) {
  summary_df <- data.frame(
    Sheet_Number = seq_along(sheet_names),
    Sheet_Name   = sheet_names,
    Description  = descriptions,
    stringsAsFactors = FALSE
  )
  if (!is.null(counts)) summary_df$N_entries <- counts
  addWorksheet(wb, "Index", gridLines = TRUE)
  writeData(wb, "Index", summary_df)
  apply_sheet_styles(wb, "Index", nrow(summary_df), ncol(summary_df))
  invisible(summary_df)
}


### TABLE 2 — all deregulated miRNAs (raw + adj p-values, sex-specific, all brain regions)
prepare_region_sex_data_full <- function(diff_express_df, region_code, region_name) {
  
  make_block <- function(sex_code, sex_label) {
    fc_col   <- paste0("fc_Comparison__Treatment_",         region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    rawp_col <- paste0("ttest_rawp_Comparison__Treatment_", region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    adjp_col <- paste0("ttest_adjp_Comparison__Treatment_", region_code, "_", sex_code,
                       "_vs_Control_",                      region_code, "_", sex_code)
    
    diff_express_df %>%
      dplyr::select(
        miRNA      = RNA,
        FoldChange = all_of(fc_col),
        raw_pvalue = all_of(rawp_col),
        adj_pvalue = all_of(adjp_col)
      ) %>%
      mutate(
        Region          = region_name,
        Sex             = sex_label,
        log2FC          = log2(FoldChange),
        Significant_raw = (abs(log2FC) > 1 & raw_pvalue < 0.05),
        Significant_adj = (abs(log2FC) > 1 & adj_pvalue < 0.05),
        Direction_raw   = case_when(
          Significant_raw & log2FC > 0 ~ "Upregulated",
          Significant_raw & log2FC < 0 ~ "Downregulated",
          TRUE                         ~ "Not significant"
        ),
        Direction_adj   = case_when(
          Significant_adj & log2FC > 0 ~ "Upregulated",
          Significant_adj & log2FC < 0 ~ "Downregulated",
          TRUE                         ~ "Not significant"
        )
      ) %>%
      dplyr::select(miRNA, Region, Sex,
                    log2FC, FoldChange,
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
      
      df <- all_mirna_full %>%
        filter(Region == reg, Sex == sx) %>%
        arrange(raw_pvalue)
      
      if (nrow(df) == 0) next
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df)
      apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df), hs, bs)
      
      # Colour rows by raw-p significance (matches main figure convention);
      # adj-p significance is available as a column for reviewers.
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
  
  # Summary sheet — counts under both thresholds
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
  
  # Move the Summary sheet to the front of the workbook
  summary_idx <- which(names(wb) == "Summary")
  worksheetOrder(wb) <- c(summary_idx, setdiff(seq_along(names(wb)), summary_idx))
  
  saveWorkbook(wb, fname, overwrite = TRUE)
}

#export_table2_all_mirnas()

# ---------------------- miRNA expression across brainstem (PMID 32188880) ----------------------
mirna_set <- c(
  "miR-223-3p"  = "mmu_miR_223_002295",
  "miR-146a-5p" = "mmu_miR_146a_000468",
  "miR-153-3p"  = "mmu_miR_153_001191",
  "miR-142a-5p"  = "mmu_miR_142_5p_002248",
  "miR-142a-3p"  = "mmu_miR_142_3p_000464",
  "miR-409-3p" = "mmu_miR_409_3p_002332",
  "miR-433-5p" = "mmu_miR_433_5p_001078"
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
  # reverse levels so Brainstem sits at top of horizontal plot
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
  
  ### Save the plots
  #ggsave(file.path(out_dir, paste0(fname, "_brainstem.pdf")), p, width = 3.3, height = 3.0, device = cairo_pdf)
  #ggsave(file.path(out_dir, paste0(fname, "_brainstem.png")), p, width = 3.3, height = 3.0, dpi = 300, bg = "white")
})

# ---------------------- miRNA expression (overlap 3 datasets - PMID: 37217597; 40161726) ----------------------
### Direction-matched significant-set extractor
get_sig_set <- function(df, fc_col, p_col,
                        direction = c("up", "down", "both"),
                        fc_threshold = 1, p_threshold = 0.05) {
  direction <- match.arg(direction)
  fc <- df[[fc_col]]; p <- df[[p_col]]
  l2 <- ifelse(fc > 0, log2(fc), NA_real_)
  keep <- !is.na(l2) & !is.na(p) & p < p_threshold & abs(l2) > fc_threshold
  if (direction == "up")   keep <- keep & l2 >  fc_threshold
  if (direction == "down") keep <- keep & l2 < -fc_threshold
  strip_mmu(df$RNA[keep])
}

### Depleted set (derived from this study)
region_codes_depleted <- c(S = "Bregma -2.06", 
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

###  Build the three sets
depleted_this_study <- get_depleted_set(diff_express_gender_specific_75per, region_codes_depleted)

aging_up <- get_sig_set(diff_exp_external,
                        "fc_Comparison_Group__84W-microglia_vs_12W-microglia",
                        "ttest_rawp_Comparison_Group__84W-microglia_vs_12W-microglia", "up")

enrich_F <- get_sig_set(diff_exp_external,
                        "fc_Comparison_Group__8W-F-CD45_vs_8W-F-Bulk",
                        "ttest_rawp_Comparison_Group__8W-F-CD45_vs_8W-F-Bulk", "up")
enrich_M <- get_sig_set(diff_exp_external,
                        "fc_Comparison_Group__8W-M-CD45_vs_8W-M-Bulk",
                        "ttest_rawp_Comparison_Group__8W-M-CD45_vs_8W-M-Bulk", "up")
microglia_enriched <- union(enrich_F, enrich_M)

cat("Set sizes — Depleted:", length(depleted_this_study),
    "| Aging:", length(aging_up),
    "| Enriched:", length(microglia_enriched), "\n")
cat("All three:",
    paste(Reduce(intersect, list(depleted_this_study, aging_up, microglia_enriched)), collapse = ", "), "\n")


### Labeled Venn — miRNA names written into each region
venn_sets <- list(
  "Depleted microglia\n(PLX, this study)"                  = depleted_this_study,
  "Aged microglia (84W > 12W)\n(Engel et al.)"             = aging_up,
  "Microglia-enriched (8W CD45 > Bulk)\n(Walsh et al.)"    = microglia_enriched
)

# Define sets names and function for Venn diagram
A <- depleted_this_study
B <- aging_up
C <- microglia_enriched

venn_function <- function(v, max_names = 4) {
  if (length(v) == 0) return("")
  if (length(v) > max_names) return(as.character(length(v))) 
  paste(v, collapse = "\n")
}

### Region label positions
region_labels <- data.frame(
  x = c(-0.8,  0.8,  0.0,  0.0, -0.5,  0.5,  0.0),
  y = c( 0.8,  0.8, -0.8,  0.8, -0.1, -0.1,  0.2),
  lab = c(
    venn_function(setdiff(setdiff(A, B), C)),
    venn_function(setdiff(setdiff(B, A), C)),
    venn_function(setdiff(setdiff(C, A), B)),
    venn_function(setdiff(intersect(A, B), C)),
    venn_function(setdiff(intersect(A, C), B)),
    venn_function(setdiff(intersect(B, C), A)),
    venn_function(Reduce(intersect, list(A, B, C)))
  ),
  stringsAsFactors = FALSE)

venn_labeled <- ggvenn(
  venn_sets,
  fill_color      = c("#ad75a1", "#c7a469", "#6ca2a7"),
  stroke_size     = 0.3,
  set_name_size   = 3.5,
  text_size       = 0,
  show_percentage = FALSE) +
  geom_text(data = region_labels, aes(x, y, label = lab),
            inherit.aes = FALSE, size = 2.5, fontface = "italic",
            family = "Helvetica", lineheight = 0.9, color = "black") +
  ggtitle("") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                  family = "Helvetica"))

print(venn_labeled)

#ggsave("Figures and plots/PNGs/mirna_overlap_three_set_labeled.png", venn_labeled, width = 8, height = 6)
#ggsave("Figures and plots/PDFs/mirna_overlap_three_set_labeled.pdf", venn_labeled, width = 8, height = 6)
