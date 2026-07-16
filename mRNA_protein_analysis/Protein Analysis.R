# Setting the working directory and importing libraries
#setwd("//")

suppressMessages({
  library("tidyverse")
  library("ggpubr")
  library("cowplot")
  library("pheatmap")
  library("DESeq2")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("patchwork")
  library("viridis")
  library("scales")
  library("circlize")
  library("ComplexUpset")
  library("rstatix")
  library("effsize")
  library("matrixStats")
  library("igraph")
  library("UpSetR")
  library("readxl")
  library("data.table")
  library("svglite")
  library("factoextra")
  library("ggrepel")
  library("RColorBrewer")
  library("gt")
  library("knitr")
  library("gridExtra")
  library("ggsignif")
  library("openxlsx")
  library("ComplexHeatmap")
  library("grid")
  library("purrr")
  library("DXMarkers")
  library("GSEABase")
  library("msigdbr")
  library("scCustomize")
  library("ggvenn")
  library("tibble")
  library("DirichletReg")
  library("variancePartition")
  library("limma")
})

set.seed(42)

# ---------------------- GLOBAL CONFIGURATION ----------------------
### Directories
DIR_INPUT      <- "Input files"
DIR_PNG        <- "Figures and plots/PNGs"
DIR_PDF        <- "Figures and plots/PDFs"
DIR_SUPP       <- "Figures and plots/Generated matrices/For supplement"
DIR_REGION_F   <- file.path(DIR_INPUT, "brain-region-degs/sex_F")
DIR_REGION_M   <- file.path(DIR_INPUT, "brain-region-degs/sex_M")

### File paths
PATH_DIANN_PG_MATRIX   <- file.path(DIR_INPUT, "report_plx.pg_matrix.tsv")
PATH_DIANN_REPORT      <- file.path(DIR_INPUT, "report_plx.tsv")
PATH_PERSEUS_PROCESSED <- file.path(DIR_INPUT, "t-tests_plx_vs_control_new.xlsx")
PATH_PERSEUS_BEFORE    <- file.path(DIR_INPUT, "perseus_before_processing.xlsx")
PATH_PSEUDOBULK_COUNTS <- file.path(DIR_INPUT, "pseudobulk_counts_celltype.csv")
PATH_PSEUDOBULK_META   <- file.path(DIR_INPUT, "pseudobulk_metadata_celltype.csv")
PATH_CELL_NUMBER       <- file.path(DIR_INPUT, "cell_number_spatial.csv")
PATH_C1QB_CPM          <- file.path(DIR_INPUT, "C1qa_C1qb_C1qc_Csf1r_CPM.csv")
PATH_SPATIAL_META      <- file.path(DIR_INPUT, "combined_data_process_meta_allenmap.csv")
PATH_CELL_COMP         <- file.path(DIR_INPUT, "celltype_composition_after_filtering.csv")
PATH_MIRNA_CONTEXT     <- file.path(DIR_INPUT, "Conserved_Site_Context_Scores.txt")
PATH_MSIGDB_GMT        <- file.path(DIR_INPUT, "m8.all.v2026.1.Mm.symbols.gmt")

MIRTAR_FILES <- c(
  file.path(DIR_INPUT, "miRTarBase/miRTarBase_SE_R.csv"),
  file.path(DIR_INPUT, "miRTarBase/miRTarBase_SE_W.csv"),
  file.path(DIR_INPUT, "miRTarBase/miRTarBase_SE_WR.csv")
)

### Analysis thresholds
FC_THRESHOLD        <- 1.0
PVAL_THRESHOLD      <- 0.05
NEG_LOG10_P_CUTOFF  <- 1.3
COHENS_D_THRESHOLD  <- 1.0
MIRNA_PERCENTILE    <- 75
MIN_CELLS_FILTER    <- 10
MIN_SAMPLE_FRAC     <- 0.20

### Sample column vectors
PLX_F_COLS     <- paste0("PLX_F_",     1:6)
PLX_M_COLS     <- paste0("PLX_M_",     1:6)
CTRL_F_COLS    <- paste0("Control_F_", 1:4)
CTRL_M_COLS    <- paste0("Control_M_", 1:4)
ALL_SAMPLE_COLS <- c(PLX_F_COLS, PLX_M_COLS, CTRL_F_COLS, CTRL_M_COLS)

# Convenience list used by Cohen's d computation
PLX_COLS  <- list(F = PLX_F_COLS,  M = PLX_M_COLS)
CTRL_COLS <- list(F = CTRL_F_COLS, M = CTRL_M_COLS)

# Sample → group/sex mapping (used for pseudobulk metadata)
# PLX1-6  = PLX Female  |  PLX7-12  = PLX Male
# PLX13-16 = Control Female  |  PLX17-20 = Control Male
FEMALE_SAMPLES <- c(paste0("PLX", 1:6),  paste0("PLX", 13:16))
MALE_SAMPLES   <- c(paste0("PLX", 7:12), paste0("PLX", 17:20))
PLX_SAMPLES    <- paste0("PLX", 1:12)
CTRL_SAMPLES   <- paste0("PLX", 13:20)

### Cell types
CELL_TYPES <- c("Astrocyte", "GABAergic", "Glutamatergic", "Microglia", "Oligodendrocyte", "Vascular")
CELL_SHORT_NAMES <- c("astro",     "gaba",      "gluta", "microglia", "oligo",     "vascular")

# Display order for diverging bar plots (males top, females bottom)
CELL_TYPE_ORDER <- c(
  "Microglia_M",    "GABAergic_M",    "Glutamatergic_M",
  "Vascular_M",     "Astrocyte_M",    "Oligodendrocyte_M",
  "Microglia_F",    "GABAergic_F",    "Glutamatergic_F",
  "Vascular_F",     "Astrocyte_F",    "Oligodendrocyte_F")

CELL_TYPES_UPSET_ORDER <- c("Oligodendrocyte", "Vascular", "Astrocyte", "Glutamatergic",   "GABAergic", "Microglia")

CELL_TYPES_EXCLUDE_SPATIAL <- c("Choroid Plexus", "Ependymal", "MB–HB-Sero–Dopa")

### Brain regions
REGION_FILE_NAMES <- c(
  "DESeq2_all_CA1_CA2_Pyramidal_Layer.csv",
  "DESeq2_all_CA3_Pyramidal_Layer.csv",
  "DESeq2_all_Cortex_Deep_Layers.csv",
  "DESeq2_all_Cortex_Upper_Layers.csv",
  "DESeq2_all_Cortical_Subplate.csv",
  "DESeq2_all_Dentate_Gyrus.csv",
  "DESeq2_all_Hippocampus_Cortex_Layer_I.csv",
  "DESeq2_all_Hypothalamus.csv",
  "DESeq2_all_Striatum.csv",
  "DESeq2_all_Thalamus.csv",
  "DESeq2_all_Tissue_Border.csv",
  "DESeq2_all_Ventricle.csv",
  "DESeq2_all_White_Matter.csv"
)

REGION_NAMES <- gsub("DESeq2_all_|.csv", "", REGION_FILE_NAMES) %>% gsub("_", " ", .)

### miRNA identifiers
MIRNA_IDS  <- c("mmu-miR-142a-3p", "mmu-miR-142a-5p", "mmu-miR-146a-5p", "mmu-miR-223-3p", "mmu-miR-409-3p")

MIRNA_KEYS <- c("142a_3p", "142a_5p", "146a_5p", "223_3p", "409_3p")

MIRNA_LABELS <- c("miR-142a-3p", "miR-142a-5p", "miR-146a-5p", "miR-223-3p",  "miR-409-3p")

### Colors
TREATMENT_COLORS <- c(Control = "#d4751a", PLX = "#8b3a7a")
SEX_COLORS       <- c(Male    = "#254e7a", Female = "#a42e30")

CELL_TYPE_COLORS <- c(
  GABAergic        = "#0D5B78",
  Glutamatergic    = "#DA808C",
  "MB-HB-Sero-Dopa"= "#3E9E64",
  Astrocyte        = "#665C47",
  Oligodendrocyte  = "#FFB84F",
  "Choroid Plexus" = "#76464a",
  Ependymal        = "#e15c61",
  Vascular         = "#B0DFE2",
  Microglia        = "#C7AFD4"
)

DEG_SEGMENT_COLORS <- c(
  "Down-regulated"           = "#538f8f",
  "Gm/Rik (down)"            = "#a8a8a8",
  "Microglial marker (down)" = "#D3D3D3",
  "Up-regulated"             = "#D6616B",
  "Gm/Rik (up)"              = "#a8a8a8",
  "Microglial marker (up)"   = "#D3D3D3"
)

MIRNA_TARGET_COLORS <- c(Validated = "#856A67", Predicted = "#CBAFAB")

CORRELATION_COLORS <- colorRamp2(c(0.85, 0.925, 1),
                                 c("#f7fbff", "#6baed6", "#08306b"))

# Proteins highlighted in volcano plots
LABEL_PROTEINS_F <- c("C1qb", "C1qc", "Cx3cl1", "Sipa1")
LABEL_PROTEINS_M <- c("Col12a1", "Kif2b", "Ccn6", "Unc5d", "Cdc25a", "Zyg11a",
                      "Syt1", "Syt2", "Synpr", "Fxr1", "Atp5g", "Atp5l",
                      "Ndufv3", "Cox6a1")

# Genes/proteins to plot in supplement panels
GENES_TO_PLOT    <- c("Ang", "Cd163", "Cd209f", "Lzts3",
                      "Fgfbp1", "Nlrp3", "Shank3", "Vinac1")
PROTEINS_TO_PLOT <- c("C1qc", "Cdc25a", "Col12a1", "Cx3cl1",
                      "Kif2b", "Sipa1", "Unc5d", "Zyg11a")

### Plot axis limits
Y_LIM_CELL_TYPES   <- 250
Y_BREAKS_CELL_TYPES <- seq(-250, 250, 100)
Y_LIM_BRAIN_REGION  <- 400
Y_BREAKS_BRAIN_REGION <- seq(-400, 400, 100)

### Excel export style
EXCEL_FONT        <- "Helvetica"
EXCEL_HEADER_FILL <- "#4F81BD"
EXCEL_HEADER_SIZE <- 13
EXCEL_BODY_SIZE   <- 12

# ---------------------- SHARED ggplot2 THEMES & SCALE OBJECTS ----------------------
diverging_bar_theme <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      text               = element_text(family = "Helvetica", color = "black"),
      plot.title         = element_text(size = 16, face = "bold",   hjust = 0.5),
      plot.subtitle      = element_text(size = 13, face = "italic", hjust = 0.5),
      axis.text          = element_text(size = 14, color = "black"),
      axis.title.x       = element_text(size = 15, face = "bold"),
      axis.title.y       = element_text(size = 15, face = "bold"),
      legend.position    = "bottom",
      legend.title       = element_text(size = 15, face = "bold"),
      legend.text        = element_text(size = 13),
      legend.key.size    = unit(0.5, "cm"),
      legend.spacing.x   = unit(0.3, "cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

deg_fill_scale <- scale_fill_manual(
  values = DEG_SEGMENT_COLORS,
  breaks = c("Down-regulated", "Up-regulated",
             "Gm/Rik (down)", "Microglial marker (down)"),
  labels = c("Down-regulated (bona fide)", "Up-regulated (bona fide)",
             "Gm/Rik (uncharacterized)", "Microglial marker"),
  name = "Legend"
)

deg_fill_scale_simple <- scale_fill_manual(
  values = c(Down_regulated = "#538f8f", Up_regulated = "#D6616B"),
  labels = c("Down-regulated", "Up-regulated"),
  name   = "Legend"
)

# Clean label lookup: internal key → display string
CLEAN_LABELS <- c(
  Astrocyte_F       = "Astrocyte (F)",       Astrocyte_M       = "Astrocyte (M)",
  GABAergic_F       = "GABAergic (F)",       GABAergic_M       = "GABAergic (M)",
  Glutamatergic_F   = "Glutamatergic (F)",   Glutamatergic_M   = "Glutamatergic (M)",
  Microglia_F       = "Microglia (F)",       Microglia_M       = "Microglia (M)",
  Oligodendrocyte_F = "Oligodendrocyte (F)", Oligodendrocyte_M = "Oligodendrocyte (M)",
  Vascular_F        = "Vascular (F)",        Vascular_M        = "Vascular (M)"
)

# ---------------------- HELPER FUNCTIONS ----------------------

### Gene classification helpers
# Test whether gene symbols look like Gm#### or ####Rik uncharacterized genes
is_gm_rik <- function(genes) {
  grepl("^Gm\\d+|\\dRik$|^[0-9].*Rik", genes, ignore.case = FALSE)
}

# Test whether genes are likely microglial markers or uncharacterized genes.
# Used to flag genes that should be excluded from non-microglial DEG sets.
is_artifact <- function(genes, marker_list) {
  toupper(genes) %in% toupper(marker_list) | is_gm_rik(genes)
}

# Partition a gene vector into genuine DEGs, Gm/Rik, and microglial markers
classify_deg_counts <- function(genes, marker_list) {
  is_micro  <- toupper(genes) %in% toupper(marker_list)
  is_gmrik  <- !is_micro & is_gm_rik(genes)
  list(
    Genuine_DEGs    = sum(!is_micro & !is_gmrik),
    GmRik_DEGs      = sum(is_gmrik),
    Microglial_DEGs = sum(is_micro)
  )
}

### Differential expression helpers

# Annotate a DESeq2 results data frame with an Up/Down/Non-Significant column.
# @param suffix Column name suffix, e.g. "astro_f"
add_diff_expression <- function(df, suffix,
                                fc_thresh   = FC_THRESHOLD,
                                pval_thresh = PVAL_THRESHOLD) {
  col_name <- paste0("diff_expressed_", suffix)
  df[[col_name]] <- "Non-Significant"
  df[[col_name]][df$log2FoldChange >  fc_thresh & df$padj < pval_thresh] <- "Up-regulated"
  df[[col_name]][df$log2FoldChange < -fc_thresh & df$padj < pval_thresh] <- "Down-regulated"
  df
}

# Load and annotate DEG CSVs for all cell types and requested sexes.
load_and_process_degs <- function(cell_types, cell_short_names,
                                  sexes = c("M", "F", "Overall")) {
  deg_list <- list()
  for (sex in sexes) {
    for (i in seq_along(cell_types)) {
      cell_type  <- cell_types[i]
      short_name <- cell_short_names[i]
      
      list_name <- if (sex == "Overall") {
        paste0("deg_", cell_type)
      } else {
        paste0("deg_", cell_type, "_", sex)
      }
      
      file_path <- if (sex == "Overall") {
        file.path(DIR_INPUT, "DEGs - Overall",
                  sprintf("deg_results_%s.csv", cell_type))
      } else {
        file.path(DIR_INPUT, "DEGs - M vs F",
                  sprintf("deg_results_%s_%s.csv", cell_type, sex))
      }
      
      if (file.exists(file_path)) {
        deg_data  <- read.table(file_path, header = TRUE, sep = ",")
        col_suffix <- if (sex == "Overall") short_name else paste0(short_name, "_", tolower(sex))
        deg_list[[list_name]] <- add_diff_expression(deg_data, col_suffix)
      }
    }
  }
  deg_list
}

### Validated DEG extraction
# Return a named list of bona-fide gene vectors for a given direction,
# with microglial markers and Gm/Rik entries removed.
get_validated_genes <- function(spatial_deg_list, direction_filter, marker_list) {
  lapply(names(spatial_deg_list), function(ct_sex) {
    entry <- spatial_deg_list[[ct_sex]]
    degs  <- entry$df %>%
      filter(.data[[entry$col]] == direction_filter) %>%
      pull(real_gene_name)
    degs[!is_artifact(degs, marker_list)]
  }) %>%
    setNames(CLEAN_LABELS[names(spatial_deg_list)]) %>%
    .[sapply(., length) > 0]
}

# Print a sex-split shared-gene summary for a validated gene list.
check_shared <- function(gene_lists, direction_label,
                         cell_types = c("Astrocyte", "GABAergic", "Glutamatergic",
                                        "Microglia", "Oligodendrocyte", "Vascular")) {
  cat(sprintf("\n-- %s --\n", direction_label))
  map_dfr(cell_types, function(ct) {
    f      <- gene_lists[[paste0(ct, " (F)")]]; if (is.null(f)) f <- character(0)
    m      <- gene_lists[[paste0(ct, " (M)")]]; if (is.null(m)) m <- character(0)
    shared <- intersect(f, m)
    data.frame(
      Cell_Type    = ct,
      F            = length(f),
      M            = length(m),
      Shared       = length(shared),
      Shared_genes = paste(shared, collapse = ", ")
    )
  }) %>% print(width = 120)
}

### Cohen's d
# Compute Hedges-corrected Cohen's d for female and male comparisons,
# appended as new columns to the input data frame.
compute_cohens_d <- function(df, plx_cols, ctrl_cols) {
  df %>%
    rowwise() %>%
    mutate(
      cohens_d_F = cohen.d(
        c_across(all_of(plx_cols[["F"]])),
        c_across(all_of(ctrl_cols[["F"]])),
        hedges.correction = TRUE
      )$estimate,
      cohens_d_M = cohen.d(
        c_across(all_of(plx_cols[["M"]])),
        c_across(all_of(ctrl_cols[["M"]])),
        hedges.correction = TRUE
      )$estimate
    ) %>%
    ungroup()
}

#' Compute Hedges-corrected Cohen's d for two numeric vectors.
calc_cohens_d <- function(plx_vals, ctrl_vals) {
  cohen.d(plx_vals, ctrl_vals, hedges.correction = TRUE)$estimate
}

### Plotting helpers

# Add count labels to both arms of a diverging bar chart.
add_diverging_labels <- function(data_up, data_down, x_var, label_size = 5) {
  list(
    geom_text(data = data_up,
              aes(x = .data[[x_var]], y = label_y, label = Count),
              inherit.aes = FALSE, size = label_size, fontface = "bold", color = "black"),
    geom_text(data = data_down,
              aes(x = .data[[x_var]], y = label_y, label = Count),
              inherit.aes = FALSE, size = label_size, fontface = "bold", color = "black")
  )
}

### DEG counting for brain region files
# Return simple up/down DEG counts from a single region CSV.
count_degs_simple <- function(file_path,
                              pval_thresh = PVAL_THRESHOLD,
                              fc_thresh   = FC_THRESHOLD) {
  data <- read.csv(file_path)
  list(
    up   = sum(data$log2FoldChange >  fc_thresh & data$padj < pval_thresh, na.rm = TRUE),
    down = sum(data$log2FoldChange < -fc_thresh & data$padj < pval_thresh, na.rm = TRUE)
  )
}

# Return segmented (genuine / Gm/Rik / microglial) up/down counts.
load_brain_deg_segments <- function(file_path, marker_list,
                                    pval_thresh = PVAL_THRESHOLD,
                                    fc_thresh   = FC_THRESHOLD) {
  if (!file.exists(file_path)) return(NULL)
  data     <- read.csv(file_path)
  gene_col <- colnames(data)[1]
  data_sig <- data %>%
    filter(!is.na(log2FoldChange), !is.na(padj),
           abs(log2FoldChange) > fc_thresh, padj < pval_thresh) %>%
    mutate(
      gene       = .data[[gene_col]],
      is_micro   = toupper(gene) %in% toupper(marker_list),
      is_gmrik   = !is_micro & is_gm_rik(gene),
      is_genuine = !is_micro & !is_gmrik
    )
  summarise_segs <- function(df)
    list(genuine = sum(df$is_genuine), gmrik = sum(df$is_gmrik), micro = sum(df$is_micro))
  up   <- summarise_segs(data_sig %>% filter(log2FoldChange >  fc_thresh))
  down <- summarise_segs(data_sig %>% filter(log2FoldChange < -fc_thresh))
  list(
    up_genuine   = up$genuine,   up_gmrik   = up$gmrik,   up_micro   = up$micro,
    down_genuine = down$genuine, down_gmrik = down$gmrik, down_micro = down$micro
  )
}

### UpSet data builder
# Build a binary membership data frame for ComplexUpset from a validated gene list.
build_upset_df <- function(validated_list,
                           all_cts = c("Astrocyte", "GABAergic", "Glutamatergic",
                                       "Microglia", "Oligodendrocyte", "Vascular")) {
  do.call(rbind, lapply(c("F", "M"), function(sex_char) {
    sex_label <- ifelse(sex_char == "M", "Male", "Female")
    sex_paren <- paste0("(", sex_char, ")")
    sex_lists <- validated_list[grepl(sex_paren, names(validated_list), fixed = TRUE)]
    names(sex_lists) <- sub(paste0(" \\(", sex_char, "\\)$"), "", names(sex_lists))
    all_genes <- unique(unlist(sex_lists))
    if (length(all_genes) == 0) return(NULL)
    membership <- as.data.frame(
      sapply(all_cts, function(ct) all_genes %in% (sex_lists[[ct]] %||% character(0)))
    )
    membership$gene <- all_genes
    membership$Sex  <- sex_label
    membership
  }))
}

# Microglial marker formatting helper
fmt_markers <- function(x) {
  unique(paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))
}

### Excel export helpers
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
  addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:n_cols, gridExpand = TRUE)
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

# ---------------------- DATA LOADING & PREPROCESSING ----------------------
### Proteomics
report_diann_before_processing <- read.delim(
  PATH_DIANN_PG_MATRIX, header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE
)
colnames(report_diann_before_processing) <- c(
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
  ALL_SAMPLE_COLS
)

perseus_proc_plx <- read_excel(PATH_PERSEUS_PROCESSED)
colnames(perseus_proc_plx) <- c(
  ALL_SAMPLE_COLS,
  "Student's T-test Significant PLX_F_Control_F",
  "Student's T-test Significant PLX_M_Control_M",
  "Student's T-test significant",
  "Number Of Imputations",
  "-Log Student's T-test p-value PLX_F_Control_F",
  "Student's T-test q-value PLX_F_Control_F",
  "Student's T-test Difference PLX_F_Control_F",
  "Student's T-test Test statistic PLX_F_Control_F",
  "-Log Student's T-test p-value PLX_M_Control_M",
  "Student's T-test q-value PLX_M_Control_M",
  "Student's T-test Difference PLX_M_Control_M",
  "Student's T-test Test statistic PLX_M_Control_M",
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"
)

report_full_diann <- read_tsv(PATH_DIANN_REPORT)
perseus_before    <- read_excel(PATH_PERSEUS_BEFORE)

### Spatial / pseudobulk
counts   <- read.csv(PATH_PSEUDOBULK_COUNTS, row.names = 1) %>%
  rownames_to_column(var = "row_id")

metadata <- read.csv(PATH_PSEUDOBULK_META, row.names = 1) %>%
  rownames_to_column(var = "row_id") %>%
  mutate(Sex = case_when(
    Sample %in% FEMALE_SAMPLES ~ "F",
    Sample %in% MALE_SAMPLES   ~ "M"
  ))

### Spatial cell metadata
cell_number_spatial <- read.csv(PATH_CELL_NUMBER, row.names = 1) %>%
  mutate(anno_0825 = recode(anno_0825,
                            "MB–HB-Sero–Dopa"  = "MB-HB-Sero-Dopa",
                            "Oligodendyocyte"  = "Oligodendrocyte"
  ))

### miRNA target lists
all_mirnas_predicted_tsm <- read.delim(PATH_MIRNA_CONTEXT) %>%
  mutate(miRNA = sub("\\.\\d+$", "", miRNA))

### Build target lists
predicted_targets <- lapply(seq_along(MIRNA_IDS), function(i) {
  all_mirnas_predicted_tsm %>%
    filter(miRNA == MIRNA_IDS[i],
           weighted.context...score.percentile >= MIRNA_PERCENTILE) %>%
    distinct(Gene.Symbol, .keep_all = TRUE)
}) %>% setNames(MIRNA_KEYS)

miRTarBase_data   <- lapply(MIRTAR_FILES, read_csv)

validated_targets <- lapply(seq_along(MIRNA_IDS), function(i) {
  bind_rows(lapply(miRTarBase_data, function(df) filter(df, miRNA == MIRNA_IDS[i]))) %>%
    distinct(`Target Gene`, .keep_all = TRUE)
}) %>% setNames(MIRNA_KEYS)

### Inventory
cat("\n=== miRNA target list inventory ===\n")
inventory <- data.frame(
  miRNA         = MIRNA_LABELS,
  TargetScan_ID = MIRNA_IDS,
  N_predicted   = sapply(predicted_targets, nrow),
  N_validated   = sapply(validated_targets, nrow)
)
print(inventory)

### Microglial marker lists
tmp_panglaodb <- tempfile(fileext = ".tsv.gz")
download.file("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
              tmp_panglaodb, mode = "wb", quiet = TRUE)

markers_panglaodb <- read.delim(tmp_panglaodb, stringsAsFactors = FALSE) %>%
  filter(cell.type == "Microglia", grepl("Mm", species)) %>%
  pull(official.gene.symbol) %>% unique() %>% fmt_markers()

markers_cellmarker <- view_data_source(data_source = "CellMarker2.0") %>%
  mutate(across(c(cell_name, tissue_class, species),
                ~ iconv(.x, to = "UTF-8", sub = ""))) %>%
  filter(
    grepl("icrogli",                          cell_name,    ignore.case = TRUE),
    grepl("Mouse|Mus",                        species,      ignore.case = TRUE),
    grepl("Brain|Central nervous|neural|Cerebr", tissue_class, ignore.case = TRUE) |
      is.na(tissue_class)
  ) %>%
  pull(geneSymbol) %>% strsplit("[,; ]+") %>% unlist() %>%
  trimws() %>% .[nchar(.) > 0] %>% unique() %>% fmt_markers()

tabula_markers <- getGmt(PATH_MSIGDB_GMT) %>%
  .[names(.) %in% c("TABULA_MURIS_SENIS_BRAIN_MYELOID_MICROGLIAL_CELL_AGEING",
                    "TABULA_MURIS_SENIS_BRAIN_MYELOID_MACROPHAGE_AGEING")] %>%
  lapply(geneIds) %>% unlist() %>% unique()

ALL_MICROGLIAL_MARKERS <- unique(c(markers_panglaodb, markers_cellmarker, tabula_markers))

cat("PanglaoDB:", length(markers_panglaodb),
    "| CellMarker2.0:", length(markers_cellmarker),
    "| MSigDB:", length(tabula_markers),
    "| Combined:", length(ALL_MICROGLIAL_MARKERS), "\n")

# ---------------------- SPATIAL DEG ----------------------
# Load all DEG CSVs and annotate with diff_expression columns
deg_data <- load_and_process_degs(CELL_TYPES, CELL_SHORT_NAMES)
list2env(deg_data, envir = .GlobalEnv)

# Named list used throughout: entry = list(df, col)
spatial_deg_list <- list(
  Astrocyte_F       = list(df = deg_Astrocyte_F,       col = "diff_expressed_astro_f"),
  Astrocyte_M       = list(df = deg_Astrocyte_M,       col = "diff_expressed_astro_m"),
  GABAergic_F       = list(df = deg_GABAergic_F,       col = "diff_expressed_gaba_f"),
  GABAergic_M       = list(df = deg_GABAergic_M,       col = "diff_expressed_gaba_m"),
  Glutamatergic_F   = list(df = deg_Glutamatergic_F,   col = "diff_expressed_gluta_f"),
  Glutamatergic_M   = list(df = deg_Glutamatergic_M,   col = "diff_expressed_gluta_m"),
  Microglia_F       = list(df = deg_Microglia_F,       col = "diff_expressed_microglia_f"),
  Microglia_M       = list(df = deg_Microglia_M,       col = "diff_expressed_microglia_m"),
  Oligodendrocyte_F = list(df = deg_Oligodendrocyte_F, col = "diff_expressed_oligo_f"),
  Oligodendrocyte_M = list(df = deg_Oligodendrocyte_M, col = "diff_expressed_oligo_m"),
  Vascular_F        = list(df = deg_Vascular_F,        col = "diff_expressed_vascular_f"),
  Vascular_M        = list(df = deg_Vascular_M,        col = "diff_expressed_vascular_m")
)

validated_up   <- get_validated_genes(spatial_deg_list, "Up-regulated",   ALL_MICROGLIAL_MARKERS)
validated_down <- get_validated_genes(spatial_deg_list, "Down-regulated", ALL_MICROGLIAL_MARKERS)

check_shared(validated_down, "Down-regulated")
check_shared(validated_up,   "Up-regulated")

# ---------------------- PROTEOMICS QC PLOTS ----------------------
### Peptide / protein count per sample (before Perseus)
proteins_per_sample <- report_full_diann %>%
  group_by(File.Name) %>%
  summarise(
    n_proteins   = n_distinct(Protein.Group),
    n_precursors = n_distinct(Precursor.Id),
    n_peptides   = n_distinct(Stripped.Sequence),
    .groups = "drop"
  ) %>%
  mutate(
    Sample    = str_extract(File.Name, "P\\d+"),
    Condition = case_when(
      Sample %in% paste0("P", 1:12)  ~ "PLX",
      Sample %in% paste0("P", 13:20) ~ "Control",
      TRUE                            ~ "Unknown"
    ),
    Sex = case_when(
      Sample %in% paste0("P", c(7:12, 17:20)) ~ "M",
      Sample %in% paste0("P", c(1:6,  13:16)) ~ "F",
      TRUE                                     ~ "Unknown"
    )
  )

p_peptides <- proteins_per_sample %>%
  mutate(Group = factor(paste(Condition, Sex, sep = "_"),
                        levels = c("Control_F", "Control_M", "PLX_F", "PLX_M"))) %>%
  ggplot(aes(x = Group, y = n_peptides, fill = Condition, color = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, color = "black") +
  scale_fill_manual(name = "Treatment",  values = TREATMENT_COLORS) +
  scale_color_manual(name = "Treatment", values = TREATMENT_COLORS) +
  scale_x_discrete(labels = c(Control_F = "Female", Control_M = "Male",
                              PLX_F = "Female",     PLX_M = "Male")) +
  guides(fill  = guide_legend(override.aes = list(color = NA)),
         color = guide_legend(override.aes = list(color = NA))) +
  ylim(15000, 30000) +
  labs(title = "Peptide Count per Sample", x = NULL, y = "Peptide Count",
       fill = "Treatment", color = "Treatment") +
  theme_minimal(base_size = 14) +
  theme(
    text            = element_text(family = "Helvetica", size = 14, color = "black"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    axis.text.y     = element_text(size = 12, color = "black"),
    axis.title.y    = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    legend.margin   = margin(t = 45),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p_peptides)
#ggsave(file.path(DIR_PNG, "Protein QC/Peptide QC.png"),    p_peptides, width = 8, height = 4, dpi = 300)
#ggsave(file.path(DIR_PDF, "Protein QC/Peptide QC.pdf"),    p_peptides, width = 8, height = 4, dpi = 300)

### Spearman correlation heatmap (all samples)
cor_matrix <- report_diann_before_processing %>%
  select(matches("^PLX|^Control")) %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

sample_names      <- colnames(cor_matrix)
sample_annotation <- data.frame(
  Sample    = sample_names,
  Condition = ifelse(grepl("^PLX", sample_names), "PLX", "Control"),
  Sex       = ifelse(grepl("_M_",  sample_names), "Male", "Female"),
  row.names = sample_names
)

col_ha <- HeatmapAnnotation(
  Condition = sample_annotation$Condition,
  Sex       = sample_annotation$Sex,
  col       = list(Condition = TREATMENT_COLORS, Sex = SEX_COLORS),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Condition = list(title = "Condition",
                     title_gp  = gpar(fontsize = 14, fontface = "bold"),
                     labels_gp = gpar(fontsize = 12)),
    Sex       = list(title = "Sex",
                     title_gp  = gpar(fontsize = 14, fontface = "bold"),
                     labels_gp = gpar(fontsize = 12))
  )
)

row_ha <- rowAnnotation(
  Condition = sample_annotation$Condition,
  Sex       = sample_annotation$Sex,
  col       = list(Condition = TREATMENT_COLORS, Sex = SEX_COLORS),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  cor_matrix,
  name                       = "Spearman\nCorrelation",
  col                        = CORRELATION_COLORS,
  cluster_rows               = TRUE,
  cluster_columns            = TRUE,
  clustering_distance_rows   = "euclidean",
  clustering_distance_columns= "euclidean",
  clustering_method_rows     = "complete",
  clustering_method_columns  = "complete",
  top_annotation             = col_ha,
  left_annotation            = row_ha,
  rect_gp                    = gpar(col = "white", lwd = 0.5),
  show_row_names             = FALSE,
  show_column_names          = FALSE,
  column_names_rot           = 45,
  heatmap_legend_param = list(
    title     = "Spearman\nCorrelation",
    title_gp  = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 12, fontfamily = "Helvetica"),
    at        = c(0.85, 0.90, 0.95, 1.00),
    legend_height = unit(4, "cm")
  ),
  column_title    = "Sample Correlation Matrix — Proteomics",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica")
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
#pdf(file.path(DIR_PDF, "Correlations/Correlation Heatmap Proteomics.pdf"), width = 10, height = 9)
#draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()

#png(file.path(DIR_PNG, "Correlations/Correlation Heatmap Proteomics.png"), width  = 10, height = 9, units  = "in", res    = 300)
#draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()

### Intensity histogram
histogram_proteins <- perseus_proc_plx %>%
  select(matches("^PLX|^Control")) %>%
  rowwise() %>%
  mutate(Mean_log2_Intensity = mean(c_across(everything()), na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = Mean_log2_Intensity)) +
  geom_histogram(bins = 50, fill = "#4D5382", alpha = 1.0, color = "white") +
  labs(
    title    = "Distribution of Protein Intensities",
    subtitle = sprintf("n = %d proteins | log2-transformed, filtered and imputed data",
                       nrow(perseus_proc_plx)),
    x = "Mean Log2 Intensity (across all samples)",
    y = "Number of Proteins"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text            = element_text(family = "Helvetica", size = 14, color = "black"),
    plot.title      = element_text(size = 14, face = "bold",   hjust = 0.5),
    plot.subtitle   = element_text(size = 12, face = "italic", hjust = 0.5),
    axis.text       = element_text(size = 12, color = "black"),
    axis.title.x    = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y    = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

#ggsave(file.path(DIR_PNG, "Protein QC/Proteomics QC Histogram.png"), histogram_proteins, width = 8, height = 6)
#ggsave(file.path(DIR_PDF, "Protein QC/Proteomics QC Histogram.pdf"), histogram_proteins, width = 8, height = 6)

### Violin plot (intensity by group/sex)
violin_plot_protein <- perseus_proc_plx %>%
  select(matches("^PLX|^Control")) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "log2_Intensity") %>%
  filter(!is.na(log2_Intensity)) %>%
  mutate(
    Condition = ifelse(grepl("^PLX", Sample), "PLX", "Control"),
    Sex       = ifelse(grepl("_M_",  Sample), "Male", "Female"),
    Group     = paste(Condition, Sex, sep = "_")
  ) %>%
  ggplot(aes(x = Group, y = log2_Intensity, fill = Condition)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white", alpha = 0.8) +
  scale_fill_manual(name = "Treatment", values = TREATMENT_COLORS) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  labs(title = "Intensity Distribution by Condition and Sex",
       x = "Group", y = "Log2 Intensity") +
  theme_minimal(base_size = 12) +
  theme(
    text            = element_text(family = "Helvetica", size = 14, color = "black"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    axis.text.y     = element_text(size = 12, color = "black"),
    axis.title.y    = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    legend.margin   = margin(t = 45),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

#ggsave(file.path(DIR_PNG, "Protein QC/Proteomics QC Violin Plot.png"), violin_plot_protein, width = 8, height = 6)
#ggsave(file.path(DIR_PDF, "Protein QC/Proteomics QC Violin Plot.pdf"), violin_plot_protein, width = 8, height = 6)

### Proteins detected (after Perseus)
proteins_detected <- perseus_proc_plx %>%
  select(matches("^PLX|^Control")) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "N_Proteins") %>%
  mutate(
    Condition = ifelse(grepl("^PLX", Sample), "PLX", "Control"),
    Sex       = ifelse(grepl("_M_",  Sample), "Male", "Female"),
    Group     = paste(Condition, Sex, sep = "_")
  ) %>%
  group_by(Group, Condition) %>%
  summarise(Mean = mean(N_Proteins), .groups = "drop") %>%
  ggplot(aes(x = Group, y = Mean, fill = Condition)) +
  geom_col(alpha = 0.8, color = "white", width = 0.7) +
  scale_fill_manual(name = "Treatment", values = TREATMENT_COLORS) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  labs(title = "Number of Detected Proteins per Group", x = NULL, y = "Number of Proteins") +
  theme_minimal(base_size = 12) +
  theme(
    text            = element_text(family = "Helvetica", size = 14, color = "black"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    axis.text.y     = element_text(size = 12, color = "black"),
    axis.title.y    = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    legend.margin   = margin(t = 45),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

#ggsave(file.path(DIR_PNG, "Protein QC/Proteins detected QC.png"), proteins_detected, width = 8, height = 5)
#ggsave(file.path(DIR_PDF, "Protein QC/Proteins detected QC.pdf"), proteins_detected, width = 8, height = 5)

### Proteins detected (before Perseus)
sample_cols_before <- colnames(perseus_before)[1:20]

detected_per_sample <- perseus_before %>%
  select(all_of(sample_cols_before)) %>%
  summarise(across(everything(), ~ sum(!is.na(.) & . != "NaN"))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Detected") %>%
  mutate(
    Condition = ifelse(grepl("^PLX", Sample), "PLX", "Control"),
    Sex       = ifelse(grepl("_M_",  Sample), "Male", "Female"),
    Group     = paste(Condition, Sex, sep = "_"),
    Missing   = nrow(perseus_before) - Detected
  )

proteins_detected_before_preproc <- detected_per_sample %>%
  ggplot(aes(x = Group, y = Detected, fill = Condition)) +
  geom_boxplot(alpha = 0.8, color = "grey30", width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.6, color = "black", shape = 16,
              show.legend = FALSE) +
  scale_fill_manual(name = "Treatment", values = TREATMENT_COLORS) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  ylim(4500, 5500) +
  labs(title    = "Number of Detected Proteins per Group",
       subtitle = "Before Processing",
       x = NULL, y = "Number of Proteins") +
  theme_minimal(base_size = 12) +
  theme(
    text             = element_text(family = "Helvetica", size = 14, color = "black"),
    plot.title       = element_text(size = 16, face = "bold",   hjust = 0.5),
    plot.subtitle    = element_text(size = 14, face = "italic", hjust = 0.5),
    axis.text.x      = element_text(margin = margin(b = 25), color = "white"),
    axis.title.x     = element_blank(),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title.y     = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 12),
    legend.position  = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(proteins_detected_before_preproc)

#ggsave(file.path(DIR_PNG, "Protein QC/Proteins detected Before Proc.png"), proteins_detected_before_preproc, width = 8, height = 5)
#ggsave(file.path(DIR_PDF, "Protein QC/Proteins detected Before Proc.pdf"), proteins_detected_before_preproc, width = 8, height = 5)

# ---------------------- VOLCANO PLOTS — PROTEOMICS ----------------------
t_test_comparisons_plx <- read_excel(PATH_PERSEUS_PROCESSED)
t_test_comparisons_plx$`N: -Log Student's T-test p-value PLX_F_Control_F` <-
  as.numeric(t_test_comparisons_plx$`N: -Log Student's T-test p-value PLX_F_Control_F`)
t_test_comparisons_plx$`N: -Log Student's T-test p-value PLX_M_Control_M` <-
  as.numeric(t_test_comparisons_plx$`N: -Log Student's T-test p-value PLX_M_Control_M`)

# Annotate DEP direction
t_test_comparisons_plx <- t_test_comparisons_plx %>%
  mutate(
    diff_expressed_F = case_when(
      `N: Student's T-test Difference PLX_F_Control_F` >  FC_THRESHOLD &
        `N: -Log Student's T-test p-value PLX_F_Control_F` > NEG_LOG10_P_CUTOFF ~ "Up-regulated",
      `N: Student's T-test Difference PLX_F_Control_F` < -FC_THRESHOLD &
        `N: -Log Student's T-test p-value PLX_F_Control_F` > NEG_LOG10_P_CUTOFF ~ "Down-regulated",
      TRUE ~ "Non-Significant"
    ),
    diff_expressed_M = case_when(
      `N: Student's T-test Difference PLX_M_Control_M` >  FC_THRESHOLD &
        `N: -Log Student's T-test p-value PLX_M_Control_M` > NEG_LOG10_P_CUTOFF ~ "Up-regulated",
      `N: Student's T-test Difference PLX_M_Control_M` < -FC_THRESHOLD &
        `N: -Log Student's T-test p-value PLX_M_Control_M` > NEG_LOG10_P_CUTOFF ~ "Down-regulated",
      TRUE ~ "Non-Significant"
    )
  )

protein_plot_data <- t_test_comparisons_plx %>%
  mutate(
    log2FC_F      = `N: Student's T-test Difference PLX_F_Control_F`,
    log10P_F      = `N: -Log Student's T-test p-value PLX_F_Control_F`,
    Significant_F = abs(log2FC_F) > FC_THRESHOLD & log10P_F > NEG_LOG10_P_CUTOFF,
    Direction_F   = case_when(
      log2FC_F >  FC_THRESHOLD & log10P_F > NEG_LOG10_P_CUTOFF ~ "Up-regulated",
      log2FC_F < -FC_THRESHOLD & log10P_F > NEG_LOG10_P_CUTOFF ~ "Down-regulated",
      TRUE                                                       ~ "Non-significant"
    ),
    log2FC_M      = `N: Student's T-test Difference PLX_M_Control_M`,
    log10P_M      = `N: -Log Student's T-test p-value PLX_M_Control_M`,
    Significant_M = abs(log2FC_M) > FC_THRESHOLD & log10P_M > NEG_LOG10_P_CUTOFF,
    Direction_M   = case_when(
      log2FC_M >  FC_THRESHOLD & log10P_M > NEG_LOG10_P_CUTOFF ~ "Up-regulated",
      log2FC_M < -FC_THRESHOLD & log10P_M > NEG_LOG10_P_CUTOFF ~ "Down-regulated",
      TRUE                                                       ~ "Non-significant"
    ),
    label_F = ifelse(`T: Genes` %in% LABEL_PROTEINS_F & Direction_F != "Non-significant",
                     `T: Genes`, NA),
    label_M = ifelse(`T: Genes` %in% LABEL_PROTEINS_M & Direction_M != "Non-significant",
                     `T: Genes`, NA)
  )

# Shared volcano aesthetics
VOLCANO_COLORS <- c("Down-regulated" = "#538f8f",
                    "Non-significant" = "gray",
                    "Up-regulated"   = "#D6616B")

VOLCANO_THEME <- theme_minimal() +
  theme(
    text         = element_text(family = "Helvetica", size = 12, color = "black"),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

volcano_plot_f <- ggplot(protein_plot_data, aes(x = log2FC_F, y = log10P_F)) +
  geom_point(aes(color = Direction_F, alpha = Significant_F, size = Significant_F)) +
  scale_color_manual(values = VOLCANO_COLORS, name = "Legend") +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.3)) +
  scale_size_manual( values = c("TRUE" = 3,   "FALSE" = 1.5)) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = NEG_LOG10_P_CUTOFF,
             linetype = "dashed", color = "gray60") +
  geom_text_repel(aes(label = label_F), size = 5, family = "Helvetica",
                  max.overlaps = 50, box.padding = 1, point.padding = 0.3,
                  segment.color = "gray50", min.segment.length = 0, na.rm = TRUE) +
  ylim(0, 10) + xlim(-4, 4) +
  labs(x = "log2(Fold Change)", y = "-log10(raw p-value)", title = NULL) +
  VOLCANO_THEME +
  guides(alpha = "none", size = "none",
         color = guide_legend(override.aes = list(size = 3)))

volcano_plot_m <- ggplot(protein_plot_data, aes(x = log2FC_M, y = log10P_M)) +
  geom_point(aes(color = Direction_M, alpha = Significant_M, size = Significant_M)) +
  scale_color_manual(values = VOLCANO_COLORS, name = "Legend") +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.3)) +
  scale_size_manual( values = c("TRUE" = 3,   "FALSE" = 1.5)) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = NEG_LOG10_P_CUTOFF,
             linetype = "dashed", color = "gray60") +
  geom_text_repel(aes(label = label_M), size = 5, family = "Helvetica",
                  max.overlaps = 50, box.padding = 1, point.padding = 0.3,
                  segment.color = "gray50", min.segment.length = 0, na.rm = TRUE) +
  ylim(0, 10) + xlim(-4, 4) +
  labs(x = "log2(Fold Change)", y = "", title = NULL) +
  VOLCANO_THEME +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(alpha = "none", size = "none",
         color = guide_legend(override.aes = list(size = 3)))

combined_protein_plot <- (volcano_plot_f | volcano_plot_m) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Proteomic Deregulation: Sex-Specific Patterns (PLX vs Control)",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5,
                                            margin = margin(b = 65)))
  ) &
  theme(text         = element_text(family = "Helvetica", size = 12, color = "black"),
        axis.text.x  = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position  = "bottom",
        legend.title     = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 16))

print(combined_protein_plot)

#ggsave(file.path(DIR_PNG, "Volcano plots/Protein/Combined Volcano Plot Protein.png"), combined_protein_plot, width = 13, height = 8)
#ggsave(file.path(DIR_PDF, "Volcano plots/Protein/Combined Volcano Plot Protein.pdf"), combined_protein_plot, width = 13, height = 8)

# ---------------------- DEP SUMMARY BAR PLOT & PROTEIN LISTS ----------------------
# ---------------------- DEP SUMMARY PAIRED BAR PLOT ----------------------
protein_summary_long <- data.frame(
  Sample = c("Female", "Male"),
  Up_regulated = c(
    sum(t_test_comparisons_plx$diff_expressed_F == "Up-regulated",   na.rm = TRUE),
    sum(t_test_comparisons_plx$diff_expressed_M == "Up-regulated",   na.rm = TRUE)
  ),
  Down_regulated = c(
    sum(t_test_comparisons_plx$diff_expressed_F == "Down-regulated", na.rm = TRUE),
    sum(t_test_comparisons_plx$diff_expressed_M == "Down-regulated", na.rm = TRUE)
  )
) %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated),
               names_to = "Regulation", values_to = "Count") %>%
  mutate(
    Regulation = factor(Regulation,
                        levels = c("Down_regulated", "Up_regulated"),
                        labels = c("Down-regulated", "Up-regulated")),
    Sample     = factor(Sample, levels = c("Female", "Male"))
  )

protein_summary_plot <- ggplot(protein_summary_long,
                               aes(x = Sample, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_dodge(width = 0.8),
            vjust = -0.4,
            size = 7, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("Down-regulated" = DEG_SEGMENT_COLORS[["Down-regulated"]],
                               "Up-regulated"   = DEG_SEGMENT_COLORS[["Up-regulated"]]),
                    name = "Legend") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Sex",
       y = "Number of Deregulated Proteins (PLX vs Control)",
       fill = "Legend") +
  theme_minimal() +
  theme(
    text          = element_text(family = "Helvetica", size = 12, color = "black"),
    axis.text.x   = element_blank(),
    axis.text.y   = element_text(size = 16, color = "black"),
    axis.title.x  = element_text(size = 1, margin = margin(b = 5)),
    axis.title.y  = element_text(size = 16, face = "bold"),
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5,
                                 margin = margin(b = 65)),
    legend.position = "bottom",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 16)
  )

print(protein_summary_plot)
#ggsave(file.path(DIR_PNG, "Protein DEG Summary Paired Barplot.png"), protein_summary_plot, width = 8, height = 7)
#ggsave(file.path(DIR_PDF, "Protein DEG Summary Paired Barplot.pdf"), protein_summary_plot, width = 8, height = 7)

# Gene lists for downstream use
protein_up_F   <- t_test_comparisons_plx$`T: Genes`[t_test_comparisons_plx$diff_expressed_F == "Up-regulated"]
protein_down_F <- t_test_comparisons_plx$`T: Genes`[t_test_comparisons_plx$diff_expressed_F == "Down-regulated"]
protein_up_M   <- t_test_comparisons_plx$`T: Genes`[t_test_comparisons_plx$diff_expressed_M == "Up-regulated"]
protein_down_M <- t_test_comparisons_plx$`T: Genes`[t_test_comparisons_plx$diff_expressed_M == "Down-regulated"]

### Proteomics DEP lists (male/female × up/down, sorted by -log10 p-value)
protein_fc_base <- t_test_comparisons_plx %>%
  transmute(
    Gene         = `T: Genes`,
    log2FC_F     = `N: Student's T-test Difference PLX_F_Control_F`,
    neg_log10P_F = `N: -Log Student's T-test p-value PLX_F_Control_F`,
    q_value_F    = `N: Student's T-test q-value PLX_F_Control_F`,
    log2FC_M     = `N: Student's T-test Difference PLX_M_Control_M`,
    neg_log10P_M = `N: -Log Student's T-test p-value PLX_M_Control_M`,
    q_value_M    = `N: Student's T-test q-value PLX_M_Control_M`,
    diff_expressed_F,
    diff_expressed_M
  )

make_dep_table <- function(df, direction, sex) {
  fc_col  <- paste0("log2FC_",     sex)
  p_col   <- paste0("neg_log10P_", sex)
  q_col   <- paste0("q_value_",   sex)
  dir_col <- paste0("diff_expressed_", sex)
  df %>%
    filter(.data[[dir_col]] == direction) %>%
    transmute(
      Gene        = Gene,
      log2FC      = .data[[fc_col]],
      neg_log10P  = .data[[p_col]],
      q_value     = .data[[q_col]]
    ) %>%
    arrange(desc(neg_log10P))
}

male_up_proteins_fc   <- make_dep_table(protein_fc_base, "Up-regulated",   "M")
male_down_proteins_fc <- make_dep_table(protein_fc_base, "Down-regulated", "M")
female_up_proteins_fc <- make_dep_table(protein_fc_base, "Up-regulated",   "F")
female_down_proteins_fc <- make_dep_table(protein_fc_base, "Down-regulated", "F")

# ---------------------- DIVERGENCE PLOT — PROTEOMICS ----------------------
# Overlap categories
female_only_up_protein        <- setdiff(protein_up_F,   protein_up_M)
female_only_down_protein      <- setdiff(protein_down_F, protein_down_M)
male_only_up_protein          <- setdiff(protein_up_M,   protein_up_F)
male_only_down_protein        <- setdiff(protein_down_M, protein_down_F)
divergent_f_up_m_down_protein <- intersect(protein_up_F,   protein_down_M)
divergent_m_up_f_down_protein <- intersect(protein_up_M,   protein_down_F)
overlapping_proteins_up       <- intersect(protein_up_M,   protein_up_F)
overlapping_proteins_down     <- intersect(protein_down_M, protein_down_F)

plot_data_divergence <- perseus_proc_plx %>%
  select(
    Genes      = Genes,
    Female_FC  = `Student's T-test Difference PLX_F_Control_F`,
    Male_FC    = `Student's T-test Difference PLX_M_Control_M`,
    Female_pval= `-Log Student's T-test p-value PLX_F_Control_F`,
    Male_pval  = `-Log Student's T-test p-value PLX_M_Control_M`
  ) %>%
  filter(!is.na(Female_FC) & !is.na(Male_FC)) %>%
  group_by(Genes) %>%
  arrange(desc(abs(Female_FC) + abs(Male_FC))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    Category = case_when(
      Genes %in% overlapping_proteins_up       ~ "Concordant UP",
      Genes %in% overlapping_proteins_down     ~ "Concordant DOWN",
      Genes %in% divergent_f_up_m_down_protein ~ "Divergent: F\u2191 M\u2193",
      Genes %in% divergent_m_up_f_down_protein ~ "Divergent: M\u2191 F\u2193",
      Genes %in% female_only_up_protein        ~ "Female Only UP",
      Genes %in% female_only_down_protein      ~ "Female Only DOWN",
      Genes %in% male_only_up_protein          ~ "Male Only UP",
      Genes %in% male_only_down_protein        ~ "Male Only DOWN",
      TRUE                                     ~ "Not Significant"
    ),
    Main_Category = case_when(
      Genes %in% overlapping_proteins_up                               ~ "Concordant UP",
      Genes %in% overlapping_proteins_down                             ~ "Concordant DOWN",
      Genes %in% c(divergent_f_up_m_down_protein, divergent_m_up_f_down_protein) ~ "Sex-Divergent",
      Genes %in% c(female_only_up_protein, female_only_down_protein)   ~ "Female-Specific",
      Genes %in% c(male_only_up_protein, male_only_down_protein)       ~ "Male-Specific",
      TRUE                                                             ~ "Not Significant"
    ),
    Distance_from_diagonal = abs(Female_FC - Male_FC) / sqrt(2),
    On_diagonal            = abs(Female_FC - Male_FC) < 0.5
  )

top_divergent <- plot_data_divergence %>%
  filter(Main_Category != "Not Significant") %>%
  arrange(desc(Distance_from_diagonal)) %>%
  slice_head(n = 15)

divergence_plot <- ggplot(plot_data_divergence, aes(x = Male_FC, y = Female_FC)) +
  geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", color = "#dddddd", linewidth = 0.5, alpha = 0.6) +
  geom_hline(yintercept = c(-FC_THRESHOLD, FC_THRESHOLD),
             linetype = "dashed", color = "#dddddd", linewidth = 0.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_point(aes(color = Distance_from_diagonal,
                 size  = Distance_from_diagonal,
                 alpha = Main_Category)) +
  scale_color_gradientn(
    colors = c("#4575B4", "#91BFDB", "#FFFFBF", "#FDAE61", "#D73027"),
    values = scales::rescale(c(0, 1, 2, 4, 8)),
    name   = expression(bold(atop("Sex-difference", "|" * Delta * " FC|"))),
    guide  = guide_colorbar(barwidth = 1.5, barheight = 15)
  ) +
  scale_size_continuous(range = c(1, 6), guide = "none") +
  scale_alpha_manual(
    values = c(
      "Concordant UP"   = 0.9, "Concordant DOWN" = 0.9,
      "Sex-Divergent"   = 0.9, "Female-Specific" = 0.8,
      "Male-Specific"   = 0.8, "Not Significant" = 0.2
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  labs(title    = "Sex comparison of the proteomic response to PLX",
       subtitle = "Distance from diagonal = Degree of sex-biased magnitude",
       x = "Male log2 Fold Change",
       y = "Female log2 Fold Change") +
  theme_minimal(base_size = 14) +
  theme(
    text         = element_text(family = "Helvetica", color = "black"),
    plot.title   = element_text(size = 20, face = "bold",   hjust = 0.5),
    plot.subtitle= element_text(size = 18, hjust = 0.5, color = "gray30"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.position  = "right",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
  ) +
  annotate("text", x = 3,  y = 4.5,  label = "Similar response in both sexes\n(near diagonal)",
           size = 6, fontface = "italic", color = "gray40") +
  annotate("text", x = 3,  y = -4.5, label = "Sex-biased magnitude\n(far from diagonal)",
           size = 6, fontface = "italic", color = "gray40")

divergence_plot_labeled <- divergence_plot +
  geom_text_repel(data = top_divergent, aes(label = Genes),
                  size = 5.5, fontface = "bold", max.overlaps = 15,
                  box.padding = 0.5, point.padding = 0.3,
                  segment.color = "gray30", segment.size = 0.3)

print(divergence_plot_labeled)

#ggsave(file.path(DIR_PNG, "Divergence Plot Labeled.png"), divergence_plot_labeled, width = 14, height = 10, dpi = 300)
#ggsave(file.path(DIR_PDF, "Divergence Plot Labeled.pdf"), divergence_plot_labeled, width = 14, height = 10, dpi = 300)

# ---------------------- C1QB — PROTEIN & mRNA EXPRESSION PLOTS ----------------------

# Shared stat comparisons used by both C1qb plots
C1QB_COMPARISONS <- list(c("Control_M", "PLX_M"), c("Control_F", "PLX_F"))

# Perseus raw p per sex for C1qb
c1qb_int  <- perseus_proc_plx %>% filter(Genes == "C1qb")
c1qb_pers <- sapply(c("M","F"), function(sx) {
  negp <- c1qb_int[[sprintf("-Log Student's T-test p-value PLX_%s_Control_%s", sx, sx)]]
  10^(-as.numeric(negp))
})

# x positions: Control_M(1), PLX_M(2), Control_F(3), PLX_F(4)
c1qb_signif <- data.frame(
  xmin       = c(1, 3),
  xmax       = c(2, 4),
  y          = c(19, 19),
  annotation = c(sprintf("p = %.2g", c1qb_pers["M"]),
                 sprintf("p = %.2g", c1qb_pers["F"]))
)

### Protein
c1qb_protein <- perseus_proc_plx %>%
  filter(Genes == "C1qb") %>%
  select(all_of(ALL_SAMPLE_COLS)) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity") %>%
  mutate(
    Group     = ifelse(grepl("PLX", Sample), "PLX", "Control"),
    Sex       = ifelse(grepl("_F_", Sample), "F", "M"),
    Group_Sex = factor(paste(Group, Sex, sep = "_"),
                       levels = c("Control_M", "PLX_M", "Control_F", "PLX_F"))
  )

c1qb_protein_plot <- ggplot(c1qb_protein, aes(x = Group_Sex, y = Intensity, fill = Group)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.9), width = 0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "#dddddd", linewidth = 0.7) +
  geom_jitter(color = "black", position = position_jitter(width = 0.15),
              alpha = 0.6, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = TREATMENT_COLORS, name = "Treatment") +
  ggsignif::geom_signif(
    data = c1qb_signif,
    aes(xmin = xmin, xmax = xmax, annotations = annotation, y_position = y - 4),
    manual = TRUE, tip_length = 0.01, textsize = 3, inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = c(Control_M = "Control", PLX_M = "PLX",
                              Control_F = "Control", PLX_F = "PLX")) +
  ylim(0, 18) +
  labs(title = "C1QB Protein Expression", x = "", y = "Protein Intensity (log2)") +
  theme_minimal() +
  theme(text        = element_text(family = "Helvetica", size = 8, color = "black"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.title    = element_text(size = 10, face = "bold"),
        legend.text     = element_text(size = 8),
        legend.box.spacing = unit(0.1, "cm"),
        plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y    = element_text(size = 10, face = "bold"),
        axis.text.y     = element_text(size = 8, color = "black"),
        axis.text.x     = element_text(size = 1, color = "white"),
        axis.ticks.x    = element_blank()
  )

print(c1qb_protein_plot)

#ggsave(file.path(DIR_PNG, "C1qb protein plot.png"), c1qb_protein_plot, width = 4, height = 4)
#ggsave(file.path(DIR_PDF, "C1qb protein plot.pdf"), c1qb_protein_plot, width = 4, height = 4)

### mRNA
plot_data_c1qb <- read.csv(PATH_C1QB_CPM, row.names = 1) %>%
  mutate(
    Group     = factor(Group, levels = c("Control", "PLX")),
    Sex       = factor(Sex,   levels = c("M", "F")),
    Group_Sex = factor(Group_Sex,
                       levels = c("Control_M", "PLX_M", "Control_F", "PLX_F"))
  )

c1qb_plot <- ggplot(plot_data_c1qb, aes(x = Group_Sex, y = C1qb, fill = Group)) +
  geom_bar(stat = "summary", fun = "mean",
           position = position_dodge(width = 0.9), width = 0.6) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "#dddddd", linewidth = 0.7) +
  geom_jitter(color = "black", position = position_jitter(width = 0.15),
              alpha = 0.6, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = TREATMENT_COLORS, name = "Treatment",
                    breaks = c("Control", "PLX")) +
  ggsignif::geom_signif(
    comparisons = C1QB_COMPARISONS, test = "t.test", map_signif_level = function(p) sprintf("p = %.2g", p),
    tip_length = 0.01, textsize = 3) +
  scale_x_discrete(labels = c(Control_M = "Control", PLX_M = "PLX",
                              Control_F = "Control", PLX_F = "PLX")) +
  ylim(0, 130) +
  labs(
    title = expression(paste(bolditalic("C1qb"), bold(" mRNA Expression"))),
    x = "", y = "CPM (Counts Per Million)") +
  theme_minimal() +
  theme(text        = element_text(family = "Helvetica", size = 8, color = "black"),
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.box.spacing = unit(0.1, "cm"),
    plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.y    = element_text(size = 10, face = "bold"),
    axis.text.y     = element_text(size = 8, color = "black"),
    axis.text.x     = element_text(size = 1, color = "white"),
    axis.ticks.x    = element_blank()
  )


print(c1qb_plot)

#ggsave(file.path(DIR_PNG, "C1qb plot.png"), c1qb_plot, width = 4, height = 4)
#ggsave(file.path(DIR_PDF, "C1qb plot.pdf"), c1qb_plot, width = 4, height = 4)

# ---------------------- SPATIAL DEG DIVERGING BAR PLOTS ----------------------
deg_fill_scale_simple <- scale_fill_manual(
  values = c(
    Down_regulated = unname(DEG_SEGMENT_COLORS[["Down-regulated"]]),  # #538f8f
    Up_regulated   = unname(DEG_SEGMENT_COLORS[["Up-regulated"]])     # #D6616B
  ),
  labels = c("Down-regulated", "Up-regulated"),
  name   = "Legend"
)

### Unfiltered — cell types
deg_summary_ct <- map_dfr(names(spatial_deg_list), function(ct_sex) {
  entry <- spatial_deg_list[[ct_sex]]
  df    <- entry$df; col <- entry$col
  data.frame(
    Cell_Type       = ct_sex,
    Up_regulated    = sum(df[[col]] == "Up-regulated",    na.rm = TRUE),
    Down_regulated  = sum(df[[col]] == "Down-regulated",  na.rm = TRUE),
    Non_Significant = sum(df[[col]] == "Non-Significant", na.rm = TRUE),
    Total           = nrow(df)
  )
})

deg_long_unfiltered <- deg_summary_ct %>%
  select(Cell_Type, Up_regulated, Down_regulated) %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated),
               names_to = "Regulation", values_to = "Count") %>%
  mutate(
    Regulation   = factor(Regulation, levels = c("Down_regulated", "Up_regulated")),
    Count_signed = ifelse(Regulation == "Down_regulated", -Count, Count),
    Cell_Type    = factor(Cell_Type, levels = CELL_TYPE_ORDER)
  )

bar_totals_unfiltered <- deg_long_unfiltered %>%
  group_by(Cell_Type, Regulation) %>%
  summarise(bar_end = sum(Count_signed), .groups = "drop")

labels_up_unfiltered <- deg_long_unfiltered %>%
  filter(Regulation == "Up_regulated") %>%
  left_join(bar_totals_unfiltered, by = c("Cell_Type", "Regulation")) %>%
  mutate(label_y = bar_end + 15)

labels_down_unfiltered <- deg_long_unfiltered %>%
  filter(Regulation == "Down_regulated") %>%
  left_join(bar_totals_unfiltered, by = c("Cell_Type", "Regulation")) %>%
  mutate(label_y = bar_end - 24)

deg_diverging_unfiltered <- ggplot(
  deg_long_unfiltered,
  aes(x = Cell_Type, y = Count_signed, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0,   color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 6.5, color = "gray40", linewidth = 0.5) +
  add_diverging_labels(
    labels_up_unfiltered,
    labels_down_unfiltered,
    "Cell_Type",
    label_size = 6
  ) +
  scale_x_discrete(labels = gsub("_(F|M)$", "", CELL_TYPE_ORDER)) +
  deg_fill_scale_simple +
  scale_y_continuous(
    breaks = Y_BREAKS_CELL_TYPES,
    labels = abs(Y_BREAKS_CELL_TYPES),
    limits = c(-Y_LIM_CELL_TYPES, Y_LIM_CELL_TYPES)
  ) +
  coord_flip() +
  labs(
    title = "Sex-Specific DEG Patterns Across Cell Types",
    x     = "Cell Type",
    y     = "Number of DEGs"
  ) +
  diverging_bar_theme()

print(deg_diverging_unfiltered)

#ggsave(file.path(DIR_PNG, "DEG Summary Diverging Sex Specific.png"), deg_diverging_unfiltered, width = 7, height = 8, dpi = 300)
#ggsave(file.path(DIR_PDF, "DEG Summary Diverging Sex Specific.pdf"), deg_diverging_unfiltered, width = 7, height = 8, dpi = 300)


### Filtered (microglial marker + Gm/Rik) — cell types
deg_long_filtered <- map_dfr(names(spatial_deg_list), function(ct_sex) {
  entry <- spatial_deg_list[[ct_sex]]
  degs  <- entry$df %>%
    filter(.data[[entry$col]] %in% c("Up-regulated", "Down-regulated")) %>%
    select(real_gene_name, direction = .data[[entry$col]])
  
  map_dfr(c("Up-regulated", "Down-regulated"), function(dir) {
    d <- degs %>% filter(direction == dir) %>% pull(real_gene_name)
    if (length(d) == 0) return(NULL)
    counts <- classify_deg_counts(d, ALL_MICROGLIAL_MARKERS)
    data.frame(
      Cell_Type = ct_sex,
      Sex       = ifelse(grepl("_F$", ct_sex), "Female", "Male"),
      Direction = dir,
      as.data.frame(counts)
    )
  })
}) %>%
  pivot_longer(cols = c(Genuine_DEGs, GmRik_DEGs, Microglial_DEGs),
               names_to = "Segment_raw", values_to = "Count") %>%
  mutate(
    Segment = case_when(
      Segment_raw == "Microglial_DEGs" & Direction == "Down-regulated" ~ "Microglial marker (down)",
      Segment_raw == "Microglial_DEGs" & Direction == "Up-regulated"   ~ "Microglial marker (up)",
      Segment_raw == "GmRik_DEGs"      & Direction == "Down-regulated" ~ "Gm/Rik (down)",
      Segment_raw == "GmRik_DEGs"      & Direction == "Up-regulated"   ~ "Gm/Rik (up)",
      Direction == "Down-regulated"                                     ~ "Down-regulated",
      TRUE                                                              ~ "Up-regulated"
    ) %>% factor(levels = c("Down-regulated", "Gm/Rik (down)", "Microglial marker (down)",
                            "Up-regulated",  "Gm/Rik (up)",   "Microglial marker (up)")),
    Count_signed = ifelse(Direction == "Down-regulated", -Count, Count),
    Cell_Type    = factor(Cell_Type, levels = CELL_TYPE_ORDER)
  )

bar_totals_filtered <- deg_long_filtered %>%
  group_by(Cell_Type, Direction) %>%
  summarise(bar_end = sum(Count_signed), .groups = "drop")

genuine_labels_ct <- deg_long_filtered %>%
  filter(Segment_raw == "Genuine_DEGs", Count > 10) %>%
  left_join(bar_totals_filtered, by = c("Cell_Type", "Direction")) %>%
  mutate(label_y = ifelse(Direction == "Down-regulated", bar_end - 13, bar_end + 13))

deg_diverging_filtered <- ggplot(
  deg_long_filtered, aes(x = Cell_Type, y = Count_signed, fill = Segment)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0,   color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = 6.5, color = "gray40", linewidth = 0.5) +
  add_diverging_labels(
    genuine_labels_ct %>% filter(Direction == "Up-regulated"),
    genuine_labels_ct %>% filter(Direction == "Down-regulated"),
    "Cell_Type", label_size = 5
  ) +
  scale_x_discrete(labels = gsub("_(F|M)$", "", CELL_TYPE_ORDER)) +
  deg_fill_scale +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(breaks = Y_BREAKS_CELL_TYPES, labels = abs(Y_BREAKS_CELL_TYPES),
                     limits = c(-Y_LIM_CELL_TYPES, Y_LIM_CELL_TYPES)) +
  coord_flip() +
  labs(
    title    = "Sex-Specific DEG Patterns Across Cell Types",
    subtitle = paste("For clarity, numerical labels are displayed only for bona fide DEGs",
                     "\nand when more than 10 DEGs are annotated"),
    x = "Cell Type", y = "Number of DEGs"
  ) +
  diverging_bar_theme()

print(deg_diverging_filtered)

#ggsave(file.path(DIR_PNG, "DEG Diverging Both Filters.png"), deg_diverging_filtered, width = 8, height = 9, dpi = 300)
#ggsave(file.path(DIR_PDF, "DEG Diverging Both Filters.pdf"), deg_diverging_filtered, width = 8, height = 9)

# ---------------------- COMBINED VOLCANO — SPATIAL (UNFILTERED & FILTERED) ----------------------
create_combined_sex_volcano <- function(deg_df_f, deg_df_m,
                                        cell_type_name,
                                        fc_threshold   = FC_THRESHOLD,
                                        pval_threshold = PVAL_THRESHOLD,
                                        apply_filter   = FALSE,
                                        marker_list    = ALL_MICROGLIAL_MARKERS,
                                        show_y_label   = TRUE,
                                        show_x_label   = TRUE,
                                        n_label        = 5,
                                        label_only_filtered = TRUE,
                                        priority_genes = GENES_TO_PLOT,
                                        y_max          = 100) {
  
  prep <- function(df, sex_label) {
    df %>%
      mutate(
        Sex         = sex_label,
        sig_base    = padj < pval_threshold & abs(log2FoldChange) > fc_threshold,
        Significant = if (apply_filter) sig_base & !is_artifact(real_gene_name, marker_list)
        else sig_base
      ) %>%
      select(real_gene_name, log2FoldChange, padj, Sex, Significant)
  }
  
  plot_data <- bind_rows(prep(deg_df_f, "Female"), prep(deg_df_m, "Male")) %>%
    mutate(
      Group = factor(
        case_when(!Significant ~ "Non-significant",
                  Sex == "Male" ~ "Male",
                  TRUE          ~ "Female"),
        levels = c("Non-significant", "Female", "Male")
      )
    ) %>%
    arrange(Significant, Sex)
  
  do_label <- (apply_filter || !label_only_filtered) && n_label > 0
  
  label_data <- if (do_label) {
    sig_pool <- plot_data %>%
      filter(Significant, Group %in% c("Female", "Male"))
    
    top_by_padj <- sig_pool %>%
      group_by(Sex) %>%
      slice_min(padj, n = n_label, with_ties = FALSE) %>%
      ungroup()
    
    priority_hits <- sig_pool %>%
      filter(real_gene_name %in% priority_genes)
    
    bind_rows(top_by_padj, priority_hits) %>%
      distinct(real_gene_name, Sex, .keep_all = TRUE)
  } else {
    plot_data[0, ]
  }
  
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold),
               linetype = "dashed", color = "gray60") +
    geom_hline(yintercept = -log10(pval_threshold),
               linetype = "dashed", color = "gray60") +
    geom_point(aes(color = Group, alpha = Significant, size = Significant),
               position = position_jitter(width = 0.15, height = 0.05, seed = 42),
               shape = 16) +
    scale_alpha_manual(values = c("TRUE" = 0.5,  "FALSE" = 0.2)) +
    scale_size_manual( values = c("TRUE" = 1.5,  "FALSE" = 0.4)) +
    scale_color_manual(
      values = c("Non-significant" = "gray",
                 "Female"          = SEX_COLORS[["Female"]],
                 "Male"            = SEX_COLORS[["Male"]]),
      breaks = c("Male", "Female", "Non-significant"),
      name   = "Legend", drop = FALSE
    ) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, y_max)) +
    labs(title = cell_type_name,
         x = if (show_x_label) "log2(Fold Change)"        else "",
         y = if (show_y_label) "-log10(adjusted p-value)" else "") +
    theme_minimal() +
    theme(
      text          = element_text(family = "Helvetica", size = 12, color = "black"),
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x   = element_text(size = 14, color = "black"),
      axis.text.y   = if (show_y_label) element_text(size = 14, color = "black") else element_blank(),
      axis.title.x  = element_text(size = 14, face = "bold"),
      axis.title.y  = if (show_y_label) element_text(size = 14, face = "bold") else element_blank(),
      axis.ticks.y  = if (show_y_label) element_line() else element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), nrow = 1),
           alpha = "none", size = "none")
  
  if (do_label && nrow(label_data) > 0) {
    p <- p + geom_text_repel(
      data = label_data,
      aes(x = log2FoldChange, y = -log10(padj),
          label = real_gene_name, color = Group),
      size = 4, fontface = "italic", family = "Helvetica",
      max.overlaps = Inf, box.padding = 0.6, point.padding = 0.5, force = 3,
      segment.color = "gray50", segment.size = 0.3,
      min.segment.length = 0, show.legend = FALSE,
      inherit.aes = FALSE
    )
  }
  
  p
}

### Volcano panel: list(female_df, male_df, title, show_y_label, show_x_label)
VOLCANO_SPECS <- list(
  list(deg_Astrocyte_F,       deg_Astrocyte_M,       "Astrocyte",       TRUE,  FALSE),
  list(deg_GABAergic_F,       deg_GABAergic_M,       "GABAergic",       FALSE, FALSE),
  list(deg_Glutamatergic_F,   deg_Glutamatergic_M,   "Glutamatergic",   FALSE, FALSE),
  list(deg_Microglia_F,       deg_Microglia_M,       "Microglia",       FALSE, FALSE),
  list(deg_Oligodendrocyte_F, deg_Oligodendrocyte_M, "Oligodendrocyte", FALSE, FALSE),
  list(deg_Vascular_F,        deg_Vascular_M,        "Vascular",        FALSE, FALSE)
)

build_combined_volcano <- function(apply_filter, n_label = 5) {
  y_max <- if (apply_filter) 65 else 100
  
  plots <- lapply(VOLCANO_SPECS, function(s)
    create_combined_sex_volcano(s[[1]], s[[2]], s[[3]],
                                apply_filter = apply_filter,
                                show_y_label = s[[4]],
                                show_x_label = s[[5]],
                                n_label      = n_label,
                                y_max        = y_max))
  
  caption_text <- if (apply_filter)
    "Significant genes (p < 0.05, |log2FC| > 1) colored by sex. Non-significant genes, microglial markers, and uncharacterized genes shown in gray."
  else
    "Significant genes (p < 0.05, |log2FC| > 1) are colored by sex; non-significant genes shown in gray."
  
  combined <- Reduce(`|`, plots) + plot_layout(guides = "collect")
  combined <- combined &
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 16),
      legend.margin    = margin(t = 13),
      legend.key.size  = unit(0.5, "cm"),
      text             = element_text(family = "Helvetica", size = 12, color = "black"),
      axis.text.x      = element_text(size = 16),
      axis.title.x     = element_text(size = 16, face = "bold", margin = margin(t = 5))
    ) &
    labs(x = "log2(Fold Change)")
  
  combined + plot_annotation(
    title   = "Cell-Type-Specific Transcriptional Response by Sex",
    caption = caption_text,
    theme   = theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 13, face = "italic", hjust = 0.5, margin = margin(t = 5))
    )
  )
}

combined_volcano_unfiltered <- build_combined_volcano(apply_filter = FALSE)
combined_volcano_filtered   <- build_combined_volcano(apply_filter = TRUE)

print(combined_volcano_unfiltered)
print(combined_volcano_filtered)

#ggsave(file.path(DIR_PNG, "Volcano Combined Sex Comparison.png"),  combined_volcano_unfiltered, width = 18, height = 7, dpi = 300)
#ggsave(file.path(DIR_PNG, "Combined Volcano Both Filters.png"),    combined_volcano_filtered,   width = 20, height = 7, dpi = 300)

#ggsave(file.path(DIR_PDF, "Volcano Combined Sex Comparison.pdf"),  combined_volcano_unfiltered, width = 18, height = 7, dpi = 300)
#ggsave(file.path(DIR_PDF, "Combined Volcano Both Filters.pdf"),    combined_volcano_filtered,   width = 20, height = 7, dpi = 300)

# ---------------------- BRAIN REGION DEG PLOTS ----------------------
results_br_simple <- map_dfr(seq_along(REGION_FILE_NAMES), function(i) {
  map_dfr(c("F", "M"), function(sex) {
    file_path <- file.path(ifelse(sex == "F", DIR_REGION_F, DIR_REGION_M),
                           REGION_FILE_NAMES[i])
    if (!file.exists(file_path)) return(NULL)
    counts    <- count_degs_simple(file_path)
    sex_label <- ifelse(sex == "F", "Female", "Male")
    data.frame(
      Region    = REGION_NAMES[i],
      Sex       = sex_label,
      Direction = c("Down-regulated", "Up-regulated"),
      Count     = c(counts$down, counts$up)
    )
  })
})

### Region order by total DEG count
REGION_BASE_ORDER <- results_br_simple %>%
  group_by(Region) %>% summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>% pull(Region)

REGION_SEX_ORDER <- c(rev(paste0(REGION_BASE_ORDER, "_M")),
                      rev(paste0(REGION_BASE_ORDER, "_F")))

results_br_diverging <- results_br_simple %>%
  mutate(
    Region_Sex   = factor(paste0(Region, "_", substr(Sex, 1, 1)), levels = REGION_SEX_ORDER),
    Count_signed = ifelse(Direction == "Down-regulated", -Count, Count),
    Direction    = factor(Direction, levels = c("Down-regulated", "Up-regulated"))
  )

bar_totals_br_simple <- results_br_diverging %>%
  group_by(Region_Sex, Direction) %>% summarise(bar_end = sum(Count_signed), .groups = "drop")

labels_up_br <- results_br_diverging %>%
  filter(Direction == "Up-regulated",   Count > 0) %>%
  left_join(bar_totals_br_simple, by = c("Region_Sex", "Direction")) %>%
  mutate(label_y = bar_end + 15)

labels_down_br <- results_br_diverging %>%
  filter(Direction == "Down-regulated", Count > 0) %>%
  left_join(bar_totals_br_simple, by = c("Region_Sex", "Direction")) %>%
  mutate(label_y = bar_end - 18)

deg_brain_unfiltered <- ggplot(
  results_br_diverging, aes(x = Region_Sex, y = Count_signed, fill = Direction)
) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = length(REGION_BASE_ORDER) + 0.5, color = "gray40", linewidth = 0.5) +
  add_diverging_labels(labels_up_br, labels_down_br, "Region_Sex", label_size = 4.5) +
  scale_x_discrete(labels = gsub("_(F|M)$", "", levels(results_br_diverging$Region_Sex))) +
  scale_fill_manual(values = c("Down-regulated" = "#538f8f", "Up-regulated" = "#D6616B"),
                    name = "Legend") +
  scale_y_continuous(breaks = Y_BREAKS_BRAIN_REGION, labels = abs(Y_BREAKS_BRAIN_REGION),
                     limits = c(-Y_LIM_BRAIN_REGION, Y_LIM_BRAIN_REGION)) +
  coord_flip() +
  labs(title = "Sex-Specific DEG Patterns Across Brain Regions",
       x = "Brain Region", y = "Number of DEGs") +
  diverging_bar_theme(base_size = 13) +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 18))

print(deg_brain_unfiltered)

#ggsave(file.path(DIR_PDF, "Region Analysis/Brain region DEG comparison.pdf"), deg_brain_unfiltered, width = 10, height = 8)
#ggsave(file.path(DIR_PNG, "Region Analysis/Brain region DEG comparison.png"), deg_brain_unfiltered, width = 10, height = 8)


### Segmented brain region plot (genuine / Gm/Rik / microglial)
results_br_segmented <- map_dfr(seq_along(REGION_FILE_NAMES), function(i) {
  map_dfr(c("F", "M"), function(sex) {
    file_path <- file.path(ifelse(sex == "F", DIR_REGION_F, DIR_REGION_M),
                           REGION_FILE_NAMES[i])
    if (!file.exists(file_path)) return(NULL)
    segs      <- load_brain_deg_segments(file_path, ALL_MICROGLIAL_MARKERS)
    sex_label <- ifelse(sex == "F", "Female", "Male")
    bind_rows(
      data.frame(Region = REGION_NAMES[i], Sex = sex_label, Direction = "Up-regulated",
                 Segment = c("Genuine", "GmRik", "Microglial"),
                 Count   = c(segs$up_genuine, segs$up_gmrik, segs$up_micro)),
      data.frame(Region = REGION_NAMES[i], Sex = sex_label, Direction = "Down-regulated",
                 Segment = c("Genuine", "GmRik", "Microglial"),
                 Count   = c(segs$down_genuine, segs$down_gmrik, segs$down_micro))
    )
  })
})

results_br_seg_diverging <- results_br_segmented %>%
  mutate(
    Region_Sex   = factor(paste0(Region, "_", substr(Sex, 1, 1)), levels = REGION_SEX_ORDER),
    Count_signed = ifelse(Direction == "Down-regulated", -Count, Count),
    Segment_Dir  = case_when(
      Segment == "Microglial" & Direction == "Down-regulated" ~ "Microglial marker (down)",
      Segment == "Microglial" & Direction == "Up-regulated"   ~ "Microglial marker (up)",
      Segment == "GmRik"      & Direction == "Down-regulated" ~ "Gm/Rik (down)",
      Segment == "GmRik"      & Direction == "Up-regulated"   ~ "Gm/Rik (up)",
      Direction == "Down-regulated"                           ~ "Down-regulated",
      TRUE                                                    ~ "Up-regulated"
    ) %>% factor(levels = c("Down-regulated", "Gm/Rik (down)", "Microglial marker (down)",
                            "Up-regulated",  "Gm/Rik (up)",   "Microglial marker (up)"))
  )

bar_totals_br_seg <- results_br_seg_diverging %>%
  group_by(Region_Sex, Direction) %>% summarise(bar_end = sum(Count_signed), .groups = "drop")

genuine_labels_br <- results_br_seg_diverging %>%
  filter(Segment == "Genuine", Count > 10) %>%
  left_join(bar_totals_br_seg, by = c("Region_Sex", "Direction")) %>%
  mutate(label_y = ifelse(Direction == "Down-regulated", bar_end - 18, bar_end + 18))

deg_brain_filtered <- ggplot(
  results_br_seg_diverging, aes(x = Region_Sex, y = Count_signed, fill = Segment_Dir)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_vline(xintercept = length(REGION_BASE_ORDER) + 0.5, color = "gray40", linewidth = 0.5) +
  add_diverging_labels(
    genuine_labels_br %>% filter(Direction == "Up-regulated"),
    genuine_labels_br %>% filter(Direction == "Down-regulated"),
    "Region_Sex", label_size = 4.5
  ) +
  scale_x_discrete(labels = gsub("_(F|M)$", "", levels(results_br_seg_diverging$Region_Sex))) +
  deg_fill_scale +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(breaks = Y_BREAKS_BRAIN_REGION, labels = abs(Y_BREAKS_BRAIN_REGION),
                     limits = c(-Y_LIM_BRAIN_REGION, Y_LIM_BRAIN_REGION)) +
  coord_flip() +
  labs(
    title    = "Sex-Specific DEG Patterns Across Brain Regions",
    subtitle = paste("For clarity, numerical labels are displayed only for bona fide DEGs",
                     "\nand when more than 10 DEGs are annotated"),
    x = "Brain Region", y = "Number of DEGs"
  ) +
  diverging_bar_theme(base_size = 13) +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 18))

print(deg_brain_filtered)

#ggsave(file.path(DIR_PDF, "Region Analysis/Brain region DEG Both Filters.pdf"), deg_brain_filtered, width = 10, height = 8)
#ggsave(file.path(DIR_PNG, "Region Analysis/Brain region DEG Both Filters.png"), deg_brain_filtered, width = 10, height = 8)

# ---------------------- UPSET PLOTS ----------------------

make_upset_plot <- function(validated_list, direction_label,
                            bar_colors = c(Female = "#a42e30", Male = "#254e7a")) {
  df         <- build_upset_df(validated_list)
  active_cts <- CELL_TYPES_UPSET_ORDER[
    CELL_TYPES_UPSET_ORDER %in% names(df) &
      sapply(CELL_TYPES_UPSET_ORDER, function(ct) any(df[[ct]]))
  ]
  
  ComplexUpset::upset(
    df,
    intersect = active_cts,
    base_annotations = list(
      "# DEGs" = (
        ggplot(mapping = aes(fill = Sex)) +
          geom_bar(stat = "count", position = "dodge", width = 0.7) +
          geom_text(
            stat = "count", aes(label = after_stat(count), group = Sex),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.5, fontface = "bold"
          ) +
          scale_fill_manual(values = bar_colors) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
          labs(y = "# DEGs") +
          theme(
            legend.position = "right",
            text            = element_text(family = "Helvetica"),
            axis.title.y    = element_text(face = "bold", size = 14)
          )
      )
    ),
    set_sizes = (
      upset_set_size(
        geom = geom_bar(aes(fill = group, x = group), width = 0.6)
      ) +
        scale_fill_manual(values = CELL_TYPE_COLORS, guide = "none") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ),
    sort_sets             = "ascending",
    width_ratio           = 0.15,
    min_size              = 1,
    n_intersections       = 30,
    sort_intersections_by = "cardinality",
    themes = upset_modify_themes(list(
      "intersections_matrix" = theme(text = element_text(family = "Helvetica", size = 12)),
      "overall_sizes"        = theme(text = element_text(family = "Helvetica", size = 11))
    ))
  ) +
    plot_annotation(
      title = "Bona Fide DEG Intersections Across Cell Types",
      theme = theme(plot.title = element_text(size = 14, face = "bold",
                                              hjust = 0.5, family = "Helvetica"))
    )
}

p_up   <- make_upset_plot(validated_up,   "Up-regulated")
p_down <- make_upset_plot(validated_down, "Down-regulated")

#ggsave(file.path(DIR_PNG, "UpSet_ComplexUpset_Upregulated.png"),   p_up,   width = 12, height = 6, dpi = 300)
#ggsave(file.path(DIR_PNG, "UpSet_ComplexUpset_Downregulated.png"), p_down, width = 12, height = 6, dpi = 300)

### DEG intersection tables
build_intersection_tables <- function(validated_list) {
  lapply(c("Male", "Female"), function(sex_label) {
    sex_char  <- substr(sex_label, 1, 1)
    sex_paren <- paste0("(", sex_char, ")")
    
    sex_lists <- validated_list[grepl(sex_paren, names(validated_list), fixed = TRUE)]
    names(sex_lists) <- sub(paste0(" \\(", sex_char, "\\)$"), "", names(sex_lists))
    
    all_genes <- unique(unlist(sex_lists))
    if (length(all_genes) == 0)
      return(list(summary = data.frame(), detail = data.frame()))
    
    membership <- setNames(
      sapply(all_genes, function(g)
        paste(sort(names(sex_lists)[sapply(sex_lists, function(gl) g %in% gl)]),
              collapse = " & ")),
      all_genes
    )
    
    summary_df <- as.data.frame(table(intersection = membership),
                                stringsAsFactors = FALSE) %>%
      rename(n_genes = Freq) %>%
      arrange(desc(n_genes))
    
    # One row per intersection group — matches screenshot format
    detail_df <- data.frame(
      gene         = all_genes,
      intersection = membership[all_genes],
      stringsAsFactors = FALSE
    ) %>%
      group_by(intersection) %>%
      summarise(
        n_cell_types = length(strsplit(intersection[1], " & ")[[1]]),
        n_genes      = n(),
        genes        = paste(sort(gene), collapse = ", "),
        .groups      = "drop"
      ) %>%
      arrange(desc(n_cell_types), desc(n_genes)) %>%
      select(intersection, n_cell_types, n_genes, genes)
    
    list(summary = summary_df, detail = detail_df)
  }) %>%
    setNames(c("Male", "Female"))
}

tables_up   <- build_intersection_tables(validated_up)
tables_down <- build_intersection_tables(validated_down)

# ---------------------- CELL PROPORTIONS SPATIAL ----------------------
cell_prop_sex <- cell_number_spatial %>%
  group_by(Group, Sex, anno_0825) %>%
  summarise(
    mean_prop = mean(prop) * 100,
    se_prop   = sd(prop) * 100 / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Group_Sex = interaction(Group, Sex, sep = " "))

# Microglia-specific stats
microglia_data <- cell_number_spatial %>%
  filter(anno_0825 == "Microglia") %>%
  mutate(prop_percent = prop * 100)

stat_test_microglia <- microglia_data %>%
  group_by(Sex) %>%
  t_test(prop_percent ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# Stacked bar — all cell types
plot_stacked <- cell_prop_sex %>%
  ggplot(aes(x = Group_Sex, y = mean_prop, fill = anno_0825)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.3, width = 0.6) +
  geom_text(data = cell_prop_sex %>% filter(mean_prop > 5),
            aes(label = sprintf("%.1f%%", mean_prop)),
            position = position_stack(vjust = 0.7),
            size = 7, fontface = "bold", color = "white") +
  scale_fill_manual(values = CELL_TYPE_COLORS, name = "Cell Type") +
  scale_x_discrete(labels = c("Control F" = "Control", "Control M" = "Control",
                              "PLX F"     = "PLX",     "PLX M"     = "PLX")) +
  labs(x = "", y = "Cell Type Proportion (%)") +
  theme_classic(base_size = 14) +
  theme(
    text         = element_text(family = "Helvetica", size = 12, color = "black"),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 35)),
    axis.text.x  = element_text(size = 18, face = "bold", margin = margin(t = 15), color = "black"),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 17),
    legend.position = "bottom"
  )

print(plot_stacked)

#ggsave(file.path(DIR_PDF, "Cell Type Stacked Barplot All Groups.pdf"), plot_stacked, width = 12, height = 9, dpi = 300)
#ggsave(file.path(DIR_PNG, "Cell Type Stacked Barplot All Groups.png"), plot_stacked, width = 12, height = 9, dpi = 300)

### Microglia proportion bar plot
microglia_summary <- microglia_data %>%
  group_by(Group, Sex) %>%
  summarise(
    mean_prop = mean(prop_percent),
    se_prop   = sd(prop_percent) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Group_Sex = interaction(Group, Sex, sep = " "))

stat_test_micro_plot <- stat_test_microglia %>%
  mutate(group1_label = paste("Control", Sex, sep = " "),
         group2_label = paste("PLX", Sex, sep = " "))

plot_microglia <- ggplot(microglia_summary,
                         aes(x = Group_Sex, y = mean_prop, fill = Group)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                width = 0.2, linewidth = 0.5) +
  geom_jitter(data = microglia_data %>%
                mutate(Group_Sex = interaction(Group, Sex, sep = " ")),
              aes(x = Group_Sex, y = prop_percent),
              width = 0.1, size = 2, alpha = 0.6, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", mean_prop)),
            vjust = -3.5, size = 7, fontface = "bold", color = "black") +
  geom_segment(
    data = stat_test_micro_plot,
    aes(x    = match(group1_label, levels(factor(microglia_summary$Group_Sex))),
        xend = match(group2_label, levels(factor(microglia_summary$Group_Sex))),
        y    = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.35,
        yend = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.35),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    data = stat_test_micro_plot,
    aes(x    = match(group1_label, levels(factor(microglia_summary$Group_Sex))),
        xend = match(group1_label, levels(factor(microglia_summary$Group_Sex))),
        y    = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.32,
        yend = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.35),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    data = stat_test_micro_plot,
    aes(x    = match(group2_label, levels(factor(microglia_summary$Group_Sex))),
        xend = match(group2_label, levels(factor(microglia_summary$Group_Sex))),
        y    = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.32,
        yend = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.35),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_text(
    data = stat_test_micro_plot,
    aes(x = (match(group1_label, levels(factor(microglia_summary$Group_Sex))) +
               match(group2_label, levels(factor(microglia_summary$Group_Sex)))) / 2,
        y = max(microglia_summary$mean_prop + microglia_summary$se_prop) * 1.4,
        label = sprintf("p = %.2g", p.adj)),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(Control = "#C7AFD4", PLX = "#C7AFD4"), name = "Treatment") +
  scale_x_discrete(labels = c("Control F" = "Control", "PLX F" = "PLX",
                              "Control M" = "Control", "PLX M" = "PLX")) +
  ylim(0, 15) +
  labs(x = "", y = "Microglia Proportion (%)") +
  theme_classic(base_size = 14) +
  theme(
    text         = element_text(family = "Helvetica", size = 12, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", margin = margin(t = 15), color = "black"),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )

print(plot_microglia)

#ggsave(file.path(DIR_PDF, "Microglia Proportion Barplot.pdf"), plot_microglia, width = 6, height = 8, dpi = 300)
#ggsave(file.path(DIR_PNG, "Microglia Proportion Barplot.png"), plot_microglia, width = 6, height = 8, dpi = 300)


### Spatial cell type proportions (Supp. Table 1)
prop_summary_full <- cell_number_spatial %>%
  group_by(Cell_Type = anno_0825, Group, Sex) %>%
  summarise(
    N        = n(),
    Mean_prop = mean(prop * 100),
    SD_prop   = sd(prop * 100),
    SE_prop   = SD_prop / sqrt(N),
    .groups = "drop"
  )

# t-test for all cell types (PLX vs Control)
stat_all_cts <- cell_number_spatial %>%
  group_by(Cell_Type = anno_0825, Sex) %>%
  summarise(
    p_value   = tryCatch({
      ctrl <- prop[Group == "Control"] * 100
      plx  <- prop[Group == "PLX"]     * 100
      if (length(ctrl) >= 2 && length(plx) >= 2) t.test(ctrl, plx)$p.value else NA_real_
    }, error = function(e) NA_real_),
    mean_diff = mean(prop[Group == "PLX"] * 100) - mean(prop[Group == "Control"] * 100),
    .groups = "drop"
  ) %>%
  group_by(Cell_Type) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    sig_label = case_when(
      is.na(p_adj)  ~ "n.s.",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "n.s."
    )
  )

final_table <- prop_summary_full %>%
  pivot_wider(names_from = Group,
              values_from = c(N, Mean_prop, SD_prop, SE_prop)) %>%
  left_join(stat_all_cts, by = c("Cell_Type", "Sex")) %>%
  transmute(
    Cell_Type,
    Sex,
    N_Control         = N_Control,
    Mean_prop_Control = round(Mean_prop_Control, 3),
    SD_prop_Control   = round(SD_prop_Control,   3),
    N_PLX             = N_PLX,
    Mean_prop_PLX     = round(Mean_prop_PLX,     3),
    SD_prop_PLX       = round(SD_prop_PLX,       3),
    Mean_difference   = round(mean_diff,          3),
    p_value           = round(p_value,            4),
    p_adj_BH          = round(p_adj,              4),
    Significance      = sig_label
  ) %>%
  arrange(Cell_Type, Sex)

# ---------------------- PSEUDOBULK CORRELATION HEATMAPS — SPATIAL ----------------------
cell_types_spatial <- sort(unique(metadata$anno_0825))
cor_matrices_list  <- list()
heatmap_list       <- list()

for (ct in cell_types_spatial) {
  ct_metadata <- metadata %>% filter(anno_0825 == ct)
  if (nrow(ct_metadata) < 3) next
  
  ct_rows   <- ct_metadata$row_id
  ct_counts <- counts %>%
    filter(row_id %in% ct_rows) %>%
    left_join(metadata %>% select(row_id, Sample), by = "row_id") %>%
    select(-row_id) %>%
    select(Sample, everything()) %>%
    column_to_rownames("Sample") %>%
    t() %>% as.matrix()
  
  ct_log <- log2(ct_counts + 1)
  ct_cor  <- cor(ct_log, method = "spearman", use = "pairwise.complete.obs")
  cor_matrices_list[[ct]] <- ct_cor
  
  ct_samples     <- colnames(ct_cor)
  ct_annotation  <- data.frame(
    Sample    = ct_samples,
    Treatment = ifelse(ct_samples %in% PLX_SAMPLES, "PLX", "Control"),
    Sex       = case_when(
      ct_samples %in% paste0("PLX", 1:6)   ~ "Female",
      ct_samples %in% paste0("PLX", 7:12)  ~ "Male",
      ct_samples %in% paste0("PLX", 13:16) ~ "Female",
      ct_samples %in% paste0("PLX", 17:20) ~ "Male",
      TRUE                                  ~ "Unknown"
    ),
    row.names = ct_samples
  )
  
  col_ha_ct <- HeatmapAnnotation(
    Treatment = ct_annotation$Treatment,
    Sex       = ct_annotation$Sex,
    col       = list(Treatment = TREATMENT_COLORS,
                     Sex       = c(Male = "#384E77", Female = "#B6465F")),
    show_annotation_name = FALSE, show_legend = FALSE
  )
  row_ha_ct <- rowAnnotation(
    Treatment = ct_annotation$Treatment,
    Sex       = ct_annotation$Sex,
    col       = list(Treatment = TREATMENT_COLORS,
                     Sex       = c(Male = "#384E77", Female = "#B6465F")),
    show_annotation_name = FALSE, show_legend = FALSE
  )
  
  heatmap_list[[ct]] <- Heatmap(
    ct_cor,
    name                        = paste0("Correlation_", gsub("[^A-Za-z0-9]", "_", ct)),
    col                         = CORRELATION_COLORS,
    cluster_rows                = TRUE,
    cluster_columns             = TRUE,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows      = "complete",
    clustering_method_columns   = "complete",
    top_annotation              = col_ha_ct,
    left_annotation             = row_ha_ct,
    rect_gp                     = gpar(col = "white", lwd = 0.5),
    show_row_names              = FALSE,
    show_column_names           = FALSE,
    show_heatmap_legend         = FALSE,
    column_title                = ct,
    column_title_gp             = gpar(fontsize = 14, fontface = "bold")
  )
}

KEY_CELL_TYPES_HEATMAP <- c("Astrocyte", "GABAergic", "Glutamatergic",
                            "Microglia", "Oligodendyocyte", "Vascular")
available_key_celltypes <- KEY_CELL_TYPES_HEATMAP[KEY_CELL_TYPES_HEATMAP %in% names(heatmap_list)]

#pdf(file.path(DIR_PDF, "Correlations/Spatial_AllCellTypes_Combined.pdf"), width = 20, height = 14)
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                           heights = unit(c(5.5, 5.5), "inches"))))
for (i in seq_along(available_key_celltypes)) {
  ct      <- available_key_celltypes[i]
  row_idx <- ceiling(i / 3)
  col_idx <- ((i - 1) %% 3) + 1
  pushViewport(viewport(layout.pos.row = row_idx, layout.pos.col = col_idx))
  draw(heatmap_list[[ct]], newpage = FALSE)
  upViewport()
}
upViewport()

# Shared legends
pushViewport(viewport(x = 0.5, y = 0.05, width = 1, height = 0.1))
pd <- packLegend(
  Legend(title = "Spearman\nCorrelation", col_fun = CORRELATION_COLORS,
         at = c(0.85, 0.90, 0.95, 1.00),
         title_gp = gpar(fontsize = 16, fontface = "bold"),
         labels_gp = gpar(fontsize = 14), direction = "horizontal"),
  Legend(labels = c("Control", "PLX"), legend_gp = gpar(fill = TREATMENT_COLORS),
         title = "Treatment",
         title_gp = gpar(fontsize = 16, fontface = "bold"),
         labels_gp = gpar(fontsize = 14), direction = "horizontal"),
  Legend(labels = c("Male", "Female"),
         legend_gp = gpar(fill = c(SEX_COLORS["Male"], SEX_COLORS["Female"])),
         title = "Sex",
         title_gp = gpar(fontsize = 16, fontface = "bold"),
         labels_gp = gpar(fontsize = 14), direction = "horizontal"),
  direction = "horizontal", column_gap = unit(2, "cm")
)

draw(pd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "center")
upViewport()

#dev.off()

# ---------------------- miRNA TARGET ANALYSIS ----------------------

### Filter a gene vector of intersection hits
filter_mirna_hits <- function(hits, marker_list = ALL_MICROGLIAL_MARKERS) {
  hits[!is_artifact(hits, marker_list)]
}

### PROTEOMICS
perseus_proc_plx_cohen <- compute_cohens_d(perseus_proc_plx, PLX_COLS, CTRL_COLS)

genes_protein <- list(
  F = perseus_proc_plx_cohen %>% filter(cohens_d_F > COHENS_D_THRESHOLD) %>% pull(Genes) %>% unique(),
  M = perseus_proc_plx_cohen %>% filter(cohens_d_M > COHENS_D_THRESHOLD) %>% pull(Genes) %>% unique()
)

# Build long table storing both counts and gene lists (pre- and post-filter)
protein_overlaps <- expand.grid(
  miRNA = MIRNA_LABELS,
  Sex   = c("Female", "Male"),
  Type  = c("Validated", "Predicted"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    mir_key = MIRNA_KEYS[match(miRNA, MIRNA_LABELS)],
    sex_key = ifelse(Sex == "Female", "F", "M")
  ) %>%
  rowwise() %>%
  mutate(
    target_genes  = list({
      if (Type == "Validated") unique(validated_targets[[mir_key]]$`Target Gene`)
      else                     unique(predicted_targets[[mir_key]]$Gene.Symbol)
    }),
    hits_raw      = list(intersect(genes_protein[[sex_key]], target_genes)),
    hits_filtered = list(filter_mirna_hits(hits_raw)),
    n_raw         = length(hits_raw),
    n_filtered    = length(hits_filtered),
    genes_removed = list(setdiff(hits_raw, hits_filtered)),
    n_removed     = length(genes_removed)
  ) %>%
  ungroup() %>%
  select(-target_genes) %>%
  mutate(miRNA = factor(miRNA, levels = MIRNA_LABELS))

protein_overlaps %>%
  select(miRNA, Sex, Type, n_raw, n_filtered, n_removed) %>%
  print(n = Inf)

protein_plot <- protein_overlaps %>%
  select(miRNA, Sex, Type, Count = n_filtered) %>%
  ggplot(aes(x = Sex, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 1.0, width = 0.5) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold", color = "black") +
  facet_wrap(~ miRNA, nrow = 1, scales = "fixed") + 
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0))) +
  scale_fill_manual(values = MIRNA_TARGET_COLORS) +
  labs(
    title    = paste0("miRNA Target Overlap with Cohen's D > ",
                      COHENS_D_THRESHOLD, " (Filtered)"),
    x = "Sex", y = "Number of Overlapping Targets", fill = "Target Type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    text             = element_text(family = "Helvetica", color = "black"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle    = element_text(hjust = 0.5, size = 12, face = "italic"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(size = 12),
    axis.text.x      = element_text(size = 1, margin = margin(b = 40)),
    axis.text.y      = element_text(size = 12),
    axis.title       = element_text(size = 13, face = "bold")
  )

print(protein_plot)

#ggsave(file.path(DIR_PNG, "miRNA/Protein miRNA Overlap Filtered.png"), protein_plot, width = 12, height = 7, dpi = 300)
#ggsave(file.path(DIR_PDF, "miRNA/Protein miRNA Overlap Filtered.pdf"), protein_plot, width = 12, height = 7)

### SPATIAL
counts_fresh  <- read.csv(PATH_PSEUDOBULK_COUNTS, row.names = 1)
counts_matrix <- as.matrix(counts_fresh)

cell_types_mirna <- setdiff(unique(metadata$anno_0825), CELL_TYPES_EXCLUDE_SPATIAL)
results_list     <- list()
overlap_results_filtered <- list()

for (ct in cell_types_mirna) {
  for (sx in c("F", "M")) {
    
    ct_sx_idx  <- which(metadata$anno_0825 == ct & metadata$Sex == sx)
    ct_sx_meta <- metadata[ct_sx_idx, ]
    sample_ids  <- ct_sx_meta$row_id
    
    ct_sx_counts <- t(counts_matrix[sample_ids, , drop = FALSE])
    
    n_plx  <- sum(ct_sx_meta$Group == "PLX")
    n_ctrl <- sum(ct_sx_meta$Group == "Control")
    if (n_plx < 2 || n_ctrl < 2) next
    
    min_samples       <- ceiling(MIN_SAMPLE_FRAC * ncol(ct_sx_counts))
    keep              <- rowSums(ct_sx_counts >= MIN_CELLS_FILTER) >= min_samples
    ct_sx_counts_filt <- ct_sx_counts[keep, , drop = FALSE]
    
    dds      <- DESeqDataSetFromMatrix(countData = round(ct_sx_counts_filt),
                                       colData = ct_sx_meta, design = ~ Group)
    vst_data <- assay(vst(dds, blind = FALSE))
    
    plx_idx  <- which(ct_sx_meta$Group == "PLX")
    ctrl_idx <- which(ct_sx_meta$Group == "Control")
    
    cohens_d <- sapply(seq_len(nrow(vst_data)), function(i)
      calc_cohens_d(vst_data[i, plx_idx], vst_data[i, ctrl_idx])
    )
    names(cohens_d) <- rownames(vst_data)
    
    key <- paste(ct, sx, sep = "_")
    results_list[[key]] <- data.frame(
      gene         = names(cohens_d),
      cohens_d     = cohens_d,
      abs_cohens_d = abs(cohens_d),
      sex          = sx,
      cell_type    = ct
    )
  }
}

all_results_spatial <- bind_rows(results_list)
high_eff_spatial    <- all_results_spatial %>% filter(cohens_d > COHENS_D_THRESHOLD)

for (ct in cell_types_mirna) {
  genes_ct <- list(
    F = high_eff_spatial %>% filter(cell_type == ct, sex == "F") %>% pull(gene),
    M = high_eff_spatial %>% filter(cell_type == ct, sex == "M") %>% pull(gene)
  )
  if (length(genes_ct$F) == 0 && length(genes_ct$M) == 0) next
  
  ct_overlaps <- expand.grid(
    miRNA = MIRNA_LABELS,
    Sex   = c("Female", "Male"),
    Type  = c("Validated", "Predicted"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      cell_type = ct,
      mir_key   = MIRNA_KEYS[match(miRNA, MIRNA_LABELS)],
      sex_key   = ifelse(Sex == "Female", "F", "M")
    ) %>%
    
    rowwise() %>%
    mutate(
      target_genes  = list({
        if (Type == "Validated") unique(validated_targets[[mir_key]]$`Target Gene`)
        else                     unique(predicted_targets[[mir_key]]$Gene.Symbol)
      }),
      hits_raw      = list(intersect(genes_ct[[sex_key]], target_genes)),
      hits_filtered = list(filter_mirna_hits(hits_raw)),
      n_raw         = length(hits_raw),
      n_filtered    = length(hits_filtered),
      genes_removed = list(setdiff(hits_raw, hits_filtered)),
      n_removed     = length(genes_removed)
    ) %>%
    ungroup() %>%
    select(-target_genes)
  
  overlap_results_filtered[[ct]] <- ct_overlaps
}

all_overlaps_spatial_filtered <- bind_rows(overlap_results_filtered)

# Audit: any removals?
all_overlaps_spatial_filtered %>%
  filter(n_removed > 0) %>%
  select(cell_type, miRNA, Sex, Type, n_raw, n_filtered, n_removed, genes_removed) %>%
  print(width = 160)

# Plot — filtered spatial
spatial_plot_filtered <- ggplot(
  all_overlaps_spatial_filtered,
  aes(x = cell_type, y = n_filtered, fill = Type)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.5,
           color = "white", linewidth = 1.0) +
  geom_text(data = subset(all_overlaps_spatial_filtered, n_filtered > 0),
            aes(label = n_filtered), position = position_stack(vjust = 0.6),
            size = 4, fontface = "bold", color = "black") +
  facet_grid(Sex ~ miRNA, scales = "free_y") +
  scale_fill_manual(values = MIRNA_TARGET_COLORS) +
  ylim(0, 80) +
  labs(
    title    = paste0("miRNA Target Overlap with Cohen's D > ", COHENS_D_THRESHOLD,
                      " Across Cell Types (Filtered)"),
    x = "Cell Type", y = "Number of Overlapping Targets", fill = "Target Type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    text             = element_text(family = "Helvetica", color = "black"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text.x     = element_text(face = "bold", size = 16),
    strip.text.y     = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 16),
    axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_text(size = 14, face = "bold")
  )

print(spatial_plot_filtered)

#ggsave(file.path(DIR_PNG, "miRNA/Spatial miRNA Overlap Filtered.png"), spatial_plot_filtered, width = 14, height = 10, dpi = 300)
#ggsave(file.path(DIR_PDF, "miRNA/Spatial miRNA Overlap Filtered.pdf"), spatial_plot_filtered, width = 14, height = 10)


# ---------------------- SPATIAL QC PLOTS ----------------------
meta <- read.csv(PATH_SPATIAL_META, row.names = 1)

# Define sample ordering: Control_F → PLX_F → Control_M → PLX_M
ORDERED_SAMPLES <- c(paste0("PLX", 13:16), paste0("PLX", 1:6),
                     paste0("PLX", 17:20), paste0("PLX", 7:12))

plot_data_qc <- meta %>%
  select(Sample, Group, Sex, anno_0825,
         nFeature_Spatial, nCount_Spatial, percent.mito) %>%
  mutate(
    Group_Sex  = paste(Group, Sex, sep = "_"),
    Sample     = factor(Sample, levels = ORDERED_SAMPLES),
    Group      = factor(Group, levels = c("Control", "PLX")),
    Treatment_Sex = factor(paste(Group, Sex, sep = "_"),
                           levels = c("Control_M", "PLX_M", "Control_F", "PLX_F"))
  )

N_FEMALE_SAMPLES <- length(c(paste0("PLX", 13:16), paste0("PLX", 1:6)))
DIVIDER_X        <- N_FEMALE_SAMPLES + 0.5

# Common QC violin theme
QC_VIOLIN_THEME <- theme_classic(base_size = 12) +
  theme(
    text          = element_text(family = "Helvetica", color = "black"),
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.y   = element_text(size = 12),
    axis.title.y  = element_text(size = 13, face = "bold", margin = margin(r = 20)),
    legend.position  = "bottom",
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(size = 12)
  )

make_qc_violin <- function(y_var, y_label, title_label, hline_val) {
  ggplot(plot_data_qc, aes(x = Sample, y = .data[[y_var]], fill = Group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1, outlier.alpha = 0.3) +
    geom_vline(xintercept = DIVIDER_X, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_hline(yintercept = hline_val, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_fill_manual(values = TREATMENT_COLORS, name = "Treatment") +
    guides(fill = guide_legend(override.aes = list(color = NA))) +
    labs(title = title_label, x = "", y = y_label) +
    QC_VIOLIN_THEME
}

nFeature_Spatial <- make_qc_violin("nFeature_Spatial", "Number of features", "nFeature_Spatial", 200)
nCount_Spatial   <- make_qc_violin("nCount_Spatial",   "Total UMIs",         "nCount_Spatial",   300)
percent_mito     <- make_qc_violin("percent.mito",      "% mito",             "Mitochondrial Percentage", 15)

combined_qc_spatial <- (nFeature_Spatial / nCount_Spatial / percent_mito) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_qc_spatial)

#(file.path(DIR_PDF, "QC_All_Metrics_Combined.pdf"), combined_qc_spatial, width = 14, height = 10)
#ggsave(file.path(DIR_PNG, "QC_All_Metrics_Combined.png"), combined_qc_spatial, width = 14, height = 10, dpi = 300)

### QC summary statistics
summary_stats <- plot_data_qc %>%
  group_by(Sample, Group, Sex) %>%
  summarise(
    N_cells             = n(),
    Mean_nFeature       = mean(nFeature_Spatial, na.rm = TRUE),
    Median_nFeature     = median(nFeature_Spatial, na.rm = TRUE),
    SD_nFeature         = sd(nFeature_Spatial, na.rm = TRUE),
    Mean_nCount         = mean(nCount_Spatial, na.rm = TRUE),
    Median_nCount       = median(nCount_Spatial, na.rm = TRUE),
    SD_nCount           = sd(nCount_Spatial, na.rm = TRUE),
    Mean_percent_mito   = mean(percent.mito, na.rm = TRUE),
    Median_percent_mito = median(percent.mito, na.rm = TRUE),
    SD_percent_mito     = sd(percent.mito, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(Sample, ORDERED_SAMPLES))

print(summary_stats, n = Inf)

# Cell retention plot
RETENTION_DATA <- data.frame(
  Sex        = c("F", "F", "M", "M"),
  Group      = c("Control", "PLX", "Control", "PLX"),
  retention_rate = c(99.575177, 98.774915, 97.958761, 99.207598)
) %>%
  mutate(Group_Sex = factor(paste0(Group, "_", Sex),
                            levels = c("Control_F", "PLX_F", "Control_M", "PLX_M")))

retention_plot <- ggplot(RETENTION_DATA, aes(x = Group_Sex, y = retention_rate, fill = Group)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "#dddddd", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", retention_rate)),
            vjust = 1.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = TREATMENT_COLORS, name = "Treatment") +
  scale_x_discrete(labels = c(Control_M = "Control", PLX_M = "PLX",
                              Control_F = "Control", PLX_F = "PLX")) +
  labs(title = "Proportion of High-Quality Cells",
       x = "", y = "High quality cells (%)") +
  theme_minimal() +
  theme(
    text            = element_text(family = "Helvetica", size = 14, color = "black"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.margin   = margin(t = 45),
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.y    = element_text(size = 16, face = "bold"),
    axis.text.y     = element_text(size = 14, color = "black"),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  )

print(retention_plot)

#ggsave(file.path(DIR_PDF, "Spatial QC/Filtered cells QC.pdf"), retention_plot, width = 8, height = 6)
#ggsave(file.path(DIR_PNG, "Spatial QC/Filtered cells QC.png"), retention_plot, width = 8, height = 6, dpi = 300)

# Cell composition change after filtering
cell_composition_filtered <- read.csv(PATH_CELL_COMP) %>%
  arrange(diff) %>%
  mutate(simplified_celltype = factor(simplified_celltype, levels = simplified_celltype))

cell_composition_filtered_plot <- ggplot(
  cell_composition_filtered, aes(x = simplified_celltype, y = diff, fill = diff > 0)
) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f%%", diff),
                hjust = ifelse(diff > 0, -0.1, 1.1)), size = 5.5) +
  scale_fill_manual(values = c("TRUE" = "gray70", "FALSE" = "gray70"), guide = "none") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  coord_flip() +
  labs(x = "Cell Type", y = "Change in Proportion (%)",
       title = "Cell Type Composition Changes", subtitle = "After QC Filtering") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle= element_text(hjust = 0.5, size = 14, face = "italic"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dashed")
  )

print(cell_composition_filtered_plot)

#ggsave(file.path(DIR_PDF, "Spatial QC/Celltype_composition_change.pdf"), cell_composition_filtered_plot, width = 8, height = 6)
#ggsave(file.path(DIR_PNG, "Spatial QC/Celltype_composition_change.png"), cell_composition_filtered_plot, width = 8, height = 6, dpi = 300)

# ---------------------- SUPPLEMENTARY GENE / PROTEIN EXPRESSION PANELS ----------------------
GENE_YLIMS <- list(
  Ang     = c(0, 15),
  Cd163   = c(0, 30),
  Cd209f  = c(0, 15),
  Lzts3   = c(0, 90),
  Fgfbp1  = c(0, 40),
  Nlrp3   = c(0, 10),
  Shank3  = c(0, 50),
  Vinac1  = c(0, 7.5)
)

plot_gene_all_celltypes <- function(gene_name,
                                    cell_types_plot = CELL_TYPES,
                                    treatment_colors = TREATMENT_COLORS,
                                    ylim = c(0, 80)) {
  cts_raw  <- read.csv(PATH_PSEUDOBULK_COUNTS, check.names = FALSE)
  meta_raw <- read.csv(PATH_PSEUDOBULK_META,   check.names = FALSE)
  colnames(cts_raw)[1]  <- "row_id"
  colnames(meta_raw)[1] <- "row_id"
  
  # Fix typo in cell type annotation
  meta_raw$anno_0825 <- recode(meta_raw$anno_0825,
                               "Oligodendyocyte" = "Oligodendrocyte")
  
  meta_raw$Sex <- case_when(
    meta_raw$Sample %in% FEMALE_SAMPLES ~ "F",
    meta_raw$Sample %in% MALE_SAMPLES   ~ "M"
  )
  
  meta_filt <- meta_raw %>%
    filter(anno_0825 %in% cell_types_plot) %>%
    mutate(anno_0825 = factor(anno_0825, levels = cell_types_plot))
  
  if (!gene_name %in% colnames(cts_raw))
    stop(paste("Gene", gene_name, "not found in counts matrix"))
  
  merged <- meta_filt %>%
    left_join(cts_raw %>% select(row_id, all_of(gene_name)), by = "row_id") %>%
    left_join(
      cts_raw %>%
        rowwise() %>%
        mutate(lib_size = sum(c_across(-row_id), na.rm = TRUE)) %>%
        ungroup() %>%
        select(row_id, lib_size),
      by = "row_id"
    ) %>%
    mutate(CPM = (.data[[gene_name]] / lib_size) * 1e6)
  
  plot_data <- merged %>%
    select(CellType = anno_0825, Group, Sex, Sample, CPM) %>%
    mutate(
      Group     = factor(Group, levels = c("Control", "PLX")),
      Sex       = factor(Sex,   levels = c("F", "M")),
      Sex_label = ifelse(Sex == "F", "Female", "Male")
    ) %>%
    filter(!is.na(CPM) & is.finite(CPM)) %>%
    mutate(CellType = droplevels(CellType))
  
  if (nrow(plot_data) == 0) stop("No valid data points for this gene")
  
  sig_df <- expand.grid(
    CellType = levels(plot_data$CellType),
    Sex      = c("F", "M"),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      padj = {
        entry <- spatial_deg_list[[paste0(CellType, "_", Sex)]]
        if (is.null(entry)) NA_real_ else {
          v <- entry$df$padj[entry$df$real_gene_name == gene_name]
          if (length(v) == 0) NA_real_ else v[1]
        }
      }
    ) %>%
    ungroup() %>%
    mutate(
      CellType  = factor(CellType, levels = levels(plot_data$CellType)),
      Sex       = factor(Sex, levels = c("F", "M")),
      Sex_label = ifelse(Sex == "F", "Female", "Male"),
      sig_label = ifelse(is.na(padj) | padj >= 0.05, "", sprintf("p = %.2g", padj))
    ) %>%
    filter(sig_label != "")
  
  if (nrow(sig_df) > 0) {
    max_vals <- plot_data %>%
      group_by(CellType, Sex, Sex_label) %>%
      summarise(max_CPM = max(CPM, na.rm = TRUE), .groups = "drop")
    sig_df <- sig_df %>%
      left_join(max_vals, by = c("CellType", "Sex", "Sex_label")) %>%
      mutate(
        x_num = as.numeric(CellType),
        xmin  = x_num - 0.1875, xmax = x_num + 0.1875,
        y_bar  = ylim[2] * 0.88,
        y_star = ylim[2] * 0.92
      )
  }
  
  p <- ggplot(plot_data, aes(x = CellType, y = CPM, fill = Group)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.6,
                 outlier.shape = NA, color = "black", linewidth = 0.4, alpha = 0.9) +
    geom_point(aes(group = Group),
               position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
               size = 1.8, alpha = 0.6, color = "black", shape = 16, show.legend = FALSE) +
    facet_wrap(~ Sex_label, ncol = 2) +
    coord_cartesian(ylim = ylim) +
    scale_fill_manual(values = treatment_colors, name = "Treatment") +
    labs(title = bquote(bold(italic(.(gene_name)))),
         x = "Cell Type", y = "CPM (Counts Per Million)") +
    theme_classic(base_size = 14) +
    theme(
      text         = element_text(family = "Helvetica", color = "black"),
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black"),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 14)
    )
  
  if (nrow(sig_df) > 0) {
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin,
                       y = y_bar - ylim[2] * 0.015, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax,
                       y = y_bar - ylim[2] * 0.015, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax,
                       y = y_bar - max_CPM * 0.015, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_text(data = sig_df,
                aes(x = x_num, y = y_star, label = sig_label),
                inherit.aes = FALSE, size = 2.5, color = "black", fontface = "bold")
  }
  p
}

plot_protein <- function(protein_name,
                         plx_cols  = PLX_F_COLS,
                         ctrl_cols = CTRL_F_COLS,
                         treatment_colors = TREATMENT_COLORS) {
  protein_data <- perseus_proc_plx %>% filter(Genes == protein_name)
  if (nrow(protein_data) == 0) stop(paste("Protein", protein_name, "not found"))
  
  plot_data <- protein_data %>%
    select(all_of(ALL_SAMPLE_COLS)) %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity") %>%
    mutate(
      Group     = ifelse(grepl("PLX", Sample), "PLX", "Control"),
      Sex       = ifelse(grepl("_F_", Sample), "F", "M"),
      Group     = factor(Group, levels = c("Control", "PLX")),
      Sex       = factor(Sex,   levels = c("F", "M")),
      Sex_label = ifelse(Sex == "F", "Female", "Male")
    ) %>%
    filter(!is.na(Intensity) & is.finite(Intensity))
  
  if (nrow(plot_data) == 0) stop("No valid data points for this protein")
  
  sig_df <- purrr::map_dfr(c("F", "M"), function(sx) {
    col_p <- paste0("N: -Log Student's T-test p-value PLX_", sx, "_Control_", sx)
    row   <- t_test_comparisons_plx %>% filter(`T: Genes` == protein_name)
    p_val <- if (nrow(row) == 0) NA_real_ else 10^(-as.numeric(row[[col_p]][1]))
    tibble(Sex = sx, Sex_label = ifelse(sx == "F", "Female", "Male"), p_value = p_val)
  }) %>%
    mutate(sig_label = ifelse(is.na(p_value) | p_value >= 0.05, "",
                              sprintf("p = %.2g", p_value))) %>%
    filter(sig_label != "")
  
  if (nrow(sig_df) > 0) {
    max_vals <- plot_data %>%
      group_by(Sex, Sex_label) %>%
      summarise(max_Intensity = max(Intensity, na.rm = TRUE), .groups = "drop")
    sig_df <- sig_df %>%
      left_join(max_vals, by = c("Sex", "Sex_label")) %>%
      mutate(xmin = 1, xmax = 2,
             y_bar  = max_Intensity * 1.05,
             y_star = max_Intensity * 1.10)
  }
  
  p <- ggplot(plot_data, aes(x = Group, y = Intensity, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, color = "black",
                 linewidth = 0.4, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.1),
               size = 2.5, alpha = 0.6, color = "black", shape = 16,
               show.legend = FALSE) +
    facet_wrap(~ Sex_label, ncol = 2) +
    scale_fill_manual(values = treatment_colors, name = "Treatment") +
    ylim(8, 20) +
    labs(title = bquote(bold(.(protein_name))),
         x = "", y = "Protein Intensity (log2)") +
    theme_classic(base_size = 14) +
    theme(
      text             = element_text(family = "Helvetica", color = "black"),
      plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 14, face = "bold"),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(size = 12, color = "black"),
      axis.ticks.x     = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 16, face = "bold"),
      legend.text      = element_text(size = 14)
    )
  
  if (nrow(sig_df) > 0) {
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmax, y = y_bar, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin,
                       y = y_bar - max_Intensity * 0.015, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax,
                       y = y_bar - max_Intensity * 0.015, yend = y_bar),
                   inherit.aes = FALSE, linewidth = 0.4, color = "black") +
      geom_text(data = sig_df, aes(x = 1.5, y = y_star, label = sig_label),
                inherit.aes = FALSE, size = 5, color = "black", fontface = "bold")
  }
  p
}

### Generate supplementary panels
gene_plots <- lapply(GENES_TO_PLOT, function(gene) {
  cat(paste("\nProcessing", gene, "...\n"))
  tryCatch(
    plot_gene_all_celltypes(gene, ylim = GENE_YLIMS[[gene]] %||% c(0, 80)),
    error = function(e) { message(paste("Could not plot", gene, ":", e$message)); NULL }
  )
}) %>% .[!sapply(., is.null)]

genes_combined <- wrap_plots(gene_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#ggsave(file.path(DIR_PNG, "Gene expression supp p-value.png"), genes_combined, width = 16, height = 20, dpi = 300)
#ggsave(file.path(DIR_PDF, "Gene expression supp p-value.pdf"), genes_combined, width = 16, height = 20)

protein_plots <- lapply(PROTEINS_TO_PLOT, function(protein) {
  cat(paste("\nProcessing", protein, "...\n"))
  tryCatch(plot_protein(protein),
           error = function(e) { message(paste("Could not plot", protein, ":", e)); NULL })
}) %>% .[!sapply(., is.null)]

combined_proteins <- wrap_plots(protein_plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(DIR_PNG, "Protein expression supp p-value.png"), combined_proteins, width = 16, height = 18, dpi = 300)
ggsave(file.path(DIR_PDF, "Protein expression supp p-value.pdf"), combined_proteins, width = 16, height = 18)


# ---------------------- VENN DIAGRAMS (look for overlap between DEGs and DEPs for males and females)  ----------------------
### Show gene names if the region is small, else the count
venn_label <- function(v, max_names = 6) {
  if (length(v) == 0)          return("")
  if (length(v) > max_names)   return(as.character(length(v)))
  paste(sort(v), collapse = "\n")
}

# Reusable 2-set ggvenn builder (identical style for spatial + protein)
make_ggvenn_2set <- function(sets, fill_color,
                             title = "", subtitle = NULL, max_names = 6) {
  A <- sets[[1]]; B <- sets[[2]]
  region_labels <- data.frame(
    x   = c(-0.85, 0.85, 0.0),
    y   = c( 0.0,  0.0,  0.0),
    lab = c(
      venn_label(setdiff(A, B),    max_names),
      venn_label(setdiff(B, A),    max_names),
      venn_label(intersect(A, B),  max_names)
    ),
    stringsAsFactors = FALSE
  )
  
  ggvenn(
    sets,
    fill_color      = fill_color,
    stroke_size     = 0.3,
    set_name_size   = 5,
    text_size       = 0,
    show_percentage = FALSE
  ) +
    geom_text(data = region_labels, aes(x, y, label = lab),
              inherit.aes = FALSE, size = 6, fontface = "bold",
              family = "Helvetica", lineheight = 0.9, color = "black") +
    labs(title = title, subtitle = subtitle) +
    coord_fixed() +
    theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5,
                                   family = "Helvetica"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic",
                                   family = "Helvetica")
    )
}

SEX_VENN_FILL <- c(SEX_COLORS[["Female"]], SEX_COLORS[["Male"]])  # Female, Male

### SPATIAL: Female vs Male (Cohen's D > 1)
FILTER_BONA_FIDE <- TRUE

spatial_cohen_sets <- list(
  Female = high_eff_spatial %>%
    filter(sex == "F", cohens_d > COHENS_D_THRESHOLD) %>%
    { if (FILTER_BONA_FIDE) filter(., !is_artifact(gene, ALL_MICROGLIAL_MARKERS)) else . } %>%
    pull(gene) %>% unique(),
  Male = high_eff_spatial %>%
    filter(sex == "M", cohens_d > COHENS_D_THRESHOLD) %>%
    { if (FILTER_BONA_FIDE) filter(., !is_artifact(gene, ALL_MICROGLIAL_MARKERS)) else . } %>%
    pull(gene) %>% unique()
)

cat("Spatial — Female:", length(spatial_cohen_sets$Female),
    "| Male:", length(spatial_cohen_sets$Male),
    "| Shared:", length(intersect(spatial_cohen_sets$Female, spatial_cohen_sets$Male)), "\n")
cat("Shared spatial genes:",
    paste(sort(intersect(spatial_cohen_sets$Female, spatial_cohen_sets$Male)), collapse = ", "), "\n")

venn_spatial_cohen <- make_ggvenn_2set(
  spatial_cohen_sets,
  fill_color = SEX_VENN_FILL,
  title      = sprintf("Spatial genes with Cohen's D > %.1f", COHENS_D_THRESHOLD),
  subtitle   = if (FILTER_BONA_FIDE) "Bona fide genes, pooled across cell types"
  else "Pooled across cell types (full pool, = Table 7)"
)

print(venn_spatial_cohen)

#ggsave(file.path(DIR_PNG, "miRNA/Spatial CohensD Venn.png"), venn_spatial_cohen, width = 6, height = 6, dpi = 300)
#ggsave(file.path(DIR_PDF, "miRNA/Spatial CohensD Venn.pdf"), venn_spatial_cohen, width = 6, height = 6)

### PROTEIN: Female vs Male (restyled to match)
cohen_protein_sets <- list(
  Female = perseus_proc_plx_cohen %>%
    filter(cohens_d_F > COHENS_D_THRESHOLD) %>%
    pull(Genes) %>% unique() %>% na.omit() %>% as.character(),
  Male = perseus_proc_plx_cohen %>%
    filter(cohens_d_M > COHENS_D_THRESHOLD) %>%
    pull(Genes) %>% unique() %>% na.omit() %>% as.character()
)

cat("Protein — Female:", length(cohen_protein_sets$Female),
    "| Male:", length(cohen_protein_sets$Male),
    "| Shared:", length(intersect(cohen_protein_sets$Female, cohen_protein_sets$Male)), "\n")

venn_protein_cohen <- make_ggvenn_2set(
  cohen_protein_sets,
  fill_color = SEX_VENN_FILL,
  title      = sprintf("Proteins with Cohen's D > %.1f", COHENS_D_THRESHOLD),
  subtitle   = "PLX vs Control"
)

print(venn_protein_cohen)

#ggsave(file.path(DIR_PNG, "Protein QC/Protein CohensD Venn.png"), venn_protein_cohen, width = 6, height = 6, dpi = 300)
#ggsave(file.path(DIR_PDF, "Protein QC/Protein CohensD Venn.pdf"), venn_protein_cohen, width = 6, height = 6)

### PER-CELL-TYPE COHEN'S D VENNS (Female vs Male)
ct_levels <- sort(unique(high_eff_spatial$cell_type))
clean_ct  <- function(x) recode(x, "Oligodendyocyte" = "Oligodendrocyte")

cohen_sets_by_celltype <- lapply(ct_levels, function(ct) {
  get_ct <- function(sx) {
    g <- high_eff_spatial %>%
      filter(cell_type == ct, sex == sx) %>%   # already cohens_d > threshold
      { if (FILTER_BONA_FIDE) filter(., !is_artifact(gene, ALL_MICROGLIAL_MARKERS)) else . } %>%
      pull(gene)
    unique(g)
  }
  list(Female = get_ct("F"), Male = get_ct("M"))
}) %>% setNames(ct_levels)

for (ct in ct_levels) {
  s <- cohen_sets_by_celltype[[ct]]
  cat(sprintf("%-16s F:%4d  M:%4d  shared:%4d\n",
              clean_ct(ct), length(s$Female), length(s$Male),
              length(intersect(s$Female, s$Male))))
}

cohen_venn_plots <- lapply(ct_levels, function(ct) {
  s <- cohen_sets_by_celltype[[ct]]
  if (length(s$Female) == 0 && length(s$Male) == 0) {
    message("Skipping ", clean_ct(ct), " (no bona fide genes in either sex)")
    return(NULL)
  }
  make_ggvenn_2set(s, fill_color = SEX_VENN_FILL, title = clean_ct(ct))
}) %>% setNames(ct_levels) %>% .[!sapply(., is.null)]

cohen_venn_grid <- wrap_plots(cohen_venn_plots, ncol = 3) +
  plot_annotation(
    title = sprintf("Genes with Cohen's D > %.1f by Cell Type",
                    COHENS_D_THRESHOLD),
    theme = theme(plot.title = element_text(size = 16, face = "bold",
                                            hjust = 0.5, family = "Helvetica"))
  )

print(cohen_venn_grid)

#ggsave(file.path(DIR_PNG, "Spatial QC/Spatial CohensD Venn by CellType.png"), cohen_venn_grid, width = 14, height = 9, dpi = 300)
#ggsave(file.path(DIR_PDF, "Spatial QC/Spatial CohensD Venn by CellType.pdf"), cohen_venn_grid, width = 14, height = 9)


# ---------------------- EXCEL SUPPLEMENTARY TABLES ----------------------
dir.create(DIR_SUPP, recursive = TRUE, showWarnings = FALSE)

export_table3 <- function() {
  fname <- file.path(DIR_SUPP, "Supplementary_Table_3_CellType_Statistics.xlsx")
  wb    <- createWorkbook()
  addWorksheet(wb, "Cell_Type_Statistics")
  writeData(wb, "Cell_Type_Statistics", final_table)
  apply_sheet_styles(wb, "Cell_Type_Statistics", nrow(final_table), ncol(final_table))
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table4 <- function() {
  fname      <- file.path(DIR_SUPP, "Supplementary_Table_4_Cell_type_DEGs.xlsx")
  wb         <- createWorkbook()
  row_styles <- list(
    micro    = createStyle(fgFill = "#FFCCCC"),
    gmrik    = createStyle(fgFill = "#FFE5CC"),
    bonafide = createStyle(fgFill = "#CCFFCC")
  )
  
  ### Summary sheet (added first, so it sits at the front)
  summary_df <- map_dfr(names(spatial_deg_list), function(ct_sex) {
    entry <- spatial_deg_list[[ct_sex]]
    map_dfr(c("Up-regulated", "Down-regulated"), function(dir) {
      df <- entry$df %>%
        filter(.data[[entry$col]] == dir) %>%
        mutate(
          is_micro    = toupper(real_gene_name) %in% toupper(ALL_MICROGLIAL_MARKERS),
          is_gmrik    = !is_micro & is_gm_rik(real_gene_name),
          is_bonafide = !is_micro & !is_gmrik
        )
      if (nrow(df) == 0) return(NULL)
      data.frame(
        Cell_Type_Sex     = CLEAN_LABELS[ct_sex],
        Direction         = dir,
        Total_DEGs        = nrow(df),
        Bona_fide         = sum(df$is_bonafide),
        Microglial_marker = sum(df$is_micro),
        Gm_Rik            = sum(df$is_gmrik)
      )
    })
  })
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_df)
  apply_sheet_styles(wb, "Summary", nrow(summary_df), ncol(summary_df))
  
  ### Per-cell-type detail sheets
  for (ct_sex in names(spatial_deg_list)) {
    entry      <- spatial_deg_list[[ct_sex]]
    sheet_name <- substr(gsub("[()/ ]", "_", CLEAN_LABELS[ct_sex]), 1, 31)
    
    df <- entry$df %>%
      filter(.data[[entry$col]] %in% c("Up-regulated", "Down-regulated")) %>%
      select(gene = real_gene_name, log2FC = log2FoldChange, pvalue, padj,
             direction = .data[[entry$col]]) %>%
      arrange(log2FC) %>%
      mutate(
        is_microglial_marker = toupper(gene) %in% toupper(ALL_MICROGLIAL_MARKERS),
        is_gm_rik            = is_gm_rik(gene),
        is_bona_fide         = !is_microglial_marker & !is_gm_rik
      )
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
    apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df))
    
    for (style_name in names(row_styles)) {
      rows <- switch(style_name,
                     "micro"    = which(df$is_microglial_marker),
                     "gmrik"    = which(df$is_gm_rik),
                     "bonafide" = which(df$is_bona_fide)
      ) + 1
      if (length(rows) > 0)
        addStyle(wb, sheet_name, row_styles[[style_name]],
                 rows = rows, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
    }
  }
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table5 <- function() {
  fname <- file.path(DIR_SUPP, "Supplementary_Table_5_DEG_Intersections.xlsx")
  
  write_sheet <- function(wb, sheet_name, df) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
    apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df))
    gene_col  <- which(colnames(df) == "genes")
    inter_col <- which(colnames(df) == "intersection")
    if (length(gene_col)  > 0) setColWidths(wb, sheet_name, cols = gene_col,  widths = 80)
    if (length(inter_col) > 0) setColWidths(wb, sheet_name, cols = inter_col, widths = 60)
  }
  
  wb     <- createWorkbook()
  sheets <- list(
    list(name = "Male_Up",     df = tables_up$Male$detail),
    list(name = "Female_Up",   df = tables_up$Female$detail),
    list(name = "Male_Down",   df = tables_down$Male$detail),
    list(name = "Female_Down", df = tables_down$Female$detail)
  )
  for (s in sheets) write_sheet(wb, s$name, s$df)
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table6 <- function() {
  fname <- file.path(DIR_SUPP, "Supplementary_Table_6_Proteomics_DEPs.xlsx")
  wb    <- createWorkbook()
  
  pretty_cols <- c(
    Gene        = "Gene",
    log2FC      = "log2(fold change)",
    neg_log10P  = "-log10(raw p-value)"
  )
  
  protein_data <- list(
    Proteins_M_Up   = male_up_proteins_fc,
    Proteins_M_Down = male_down_proteins_fc,
    Proteins_F_Up   = female_up_proteins_fc,
    Proteins_F_Down = female_down_proteins_fc
  )
  
  for (sheet_name in names(protein_data)) {
    data <- protein_data[[sheet_name]]
    if (nrow(data) == 0) next
    
    data <- data %>% select(Gene, log2FC, neg_log10P)
    colnames(data) <- pretty_cols[colnames(data)]
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data)
    apply_sheet_styles(wb, sheet_name, nrow(data), ncol(data))
  }
  
  sheet_names <- names(wb)
  add_index_sheet(wb, sheet_names,
                  descriptions = c(
                    "Male Up-regulated proteins (sorted by -log10 raw p-value)",
                    "Male Down-regulated proteins (sorted by -log10 raw p-value)",
                    "Female Up-regulated proteins (sorted by -log10 raw p-value)",
                    "Female Down-regulated proteins (sorted by -log10 raw p-value)"
                  ),
                  counts = sapply(protein_data, nrow))
  worksheetOrder(wb) <- c(length(sheet_names) + 1, seq_along(sheet_names))
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table7 <- function() {
  fname <- file.path(DIR_SUPP,
                     "Supplementary_Table_7_CohensD_HighEffect_Pool.xlsx")
  wb    <- createWorkbook()
  
  ### Build the spatial pool (matches miRNA analysis exactly)
  spatial_pool <- all_results_spatial %>%
    filter(cohens_d > COHENS_D_THRESHOLD) %>%
    mutate(
      is_microglial_marker = toupper(gene) %in% toupper(ALL_MICROGLIAL_MARKERS),
      is_Gm_Rik            = is_gm_rik(gene),
      is_bona_fide         = !is_microglial_marker & !is_Gm_Rik
    )
  
  ### Build the proteomics pool
  protein_pool <- perseus_proc_plx_cohen %>%
    select(Gene = Genes, cohens_d_F, cohens_d_M) %>%
    pivot_longer(c(cohens_d_F, cohens_d_M), names_to = "Sex", values_to = "cohens_d") %>%
    mutate(Sex = ifelse(Sex == "cohens_d_F", "F", "M")) %>%
    filter(!is.na(cohens_d), cohens_d > COHENS_D_THRESHOLD) %>%
    arrange(desc(cohens_d)) %>%
    distinct(Gene, Sex, .keep_all = TRUE) %>%
    mutate(
      is_microglial_marker = toupper(Gene) %in% toupper(ALL_MICROGLIAL_MARKERS),
      is_Gm_Rik            = is_gm_rik(Gene),
      is_bona_fide         = !is_microglial_marker & !is_Gm_Rik
    )
  
  sheet_registry <- list()
  
  ### Helper: write a pool data frame to a single sheet
  write_pool_sheet <- function(df, sheet_name, desc) {
    if (nrow(df) == 0) return(invisible(NULL))
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
    apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df))
    
    excl_rows <- which(!df$Is_bona_fide) + 1
    if (length(excl_rows) > 0)
      addStyle(wb, sheet_name,
               createStyle(fgFill = "#FFCCCC"),
               rows = excl_rows, cols = 1:ncol(df),
               gridExpand = TRUE, stack = TRUE)
    
    sheet_registry[[sheet_name]] <<- list(desc = desc, count = nrow(df))
  }
  
  ### Per-cell-type spatial sheets (both sexes, sorted)
  for (ct in sort(unique(spatial_pool$cell_type))) {
    ct_df <- spatial_pool %>%
      filter(cell_type == ct) %>%
      transmute(
        Gene                 = gene,
        Sex                  = sex,
        Cohens_D             = round(cohens_d, 3),
        Is_microglial_marker = is_microglial_marker,
        Is_Gm_Rik            = is_Gm_Rik,
        Is_bona_fide         = is_bona_fide
      ) %>%
      arrange(Sex, desc(Cohens_D))
    
    sheet_name <- substr(paste0("Spatial_", gsub("[^A-Za-z0-9]", "", ct)), 1, 31)
    write_pool_sheet(ct_df, sheet_name,
                     paste("Spatial |", ct, "| genes with Cohen's D > 1"))
  }
  
  ### Proteomics sheets, split by sex
  for (sx in c("F", "M")) {
    sx_df <- protein_pool %>%
      filter(Sex == sx) %>%
      transmute(
        Gene                 = Gene,
        Cohens_D             = round(cohens_d, 3),
        Is_microglial_marker = is_microglial_marker,
        Is_Gm_Rik            = is_Gm_Rik,
        Is_bona_fide         = is_bona_fide
      ) %>%
      arrange(desc(Cohens_D))
    
    sx_label <- ifelse(sx == "F", "Female", "Male")
    write_pool_sheet(sx_df, paste0("Proteomics_", sx),
                     paste("Proteomics |", sx_label,
                           "| proteins with Cohen's D > 1"))
  }
  
  ### Index sheet
  all_names  <- names(sheet_registry)
  all_descs  <- sapply(sheet_registry, `[[`, "desc")
  all_counts <- sapply(sheet_registry, `[[`, "count")
  
  add_index_sheet(wb, all_names, all_descs, all_counts)
  
  n_sheets <- length(names(wb))
  worksheetOrder(wb) <- c(n_sheets, seq_len(n_sheets - 1))
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table8 <- function() {
  fname <- file.path(DIR_SUPP, "Supplementary_Table_8_miRNA_target_database_lists.xlsx")
  wb    <- createWorkbook()
  
  # Long-format predicted targets (TargetScan, percentile-filtered)
  predicted_targets_export <- map2_dfr(MIRNA_LABELS, MIRNA_KEYS, function(label, key) {
    df <- predicted_targets[[key]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% mutate(miRNA_label = label) %>%
      select(miRNA = miRNA_label, miRNA_ID = miRNA, everything())
  }) %>%
    arrange(miRNA, desc(weighted.context...score.percentile))
  
  # Long-format validated targets (miRTarBase, deduplicated by Target Gene)
  validated_targets_export <- map2_dfr(MIRNA_LABELS, MIRNA_KEYS, function(label, key) {
    df <- validated_targets[[key]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>% mutate(miRNA_label = label) %>%
      select(miRNA = miRNA_label, miRNA_ID = miRNA, everything())
  }) %>%
    arrange(miRNA)
  
  sheet_names  <- c("Predicted_targets_TargetScan", "Validated_targets_miRTarBase")
  sheet_counts <- c(nrow(predicted_targets_export),
                    nrow(validated_targets_export))
  
  # Add Index sheet first so it sits at position 1
  add_index_sheet(wb, sheet_names, descriptions = rep("", length(sheet_names)), counts = sheet_counts)
  
  addWorksheet(wb, "Predicted_targets_TargetScan", gridLines = TRUE)
  writeData(wb, "Predicted_targets_TargetScan", predicted_targets_export)
  apply_sheet_styles(wb, "Predicted_targets_TargetScan",
                     n_rows = nrow(predicted_targets_export),
                     n_cols = ncol(predicted_targets_export))
  
  addWorksheet(wb, "Validated_targets_miRTarBase", gridLines = TRUE)
  writeData(wb, "Validated_targets_miRTarBase", validated_targets_export)
  apply_sheet_styles(wb, "Validated_targets_miRTarBase",
                     n_rows = nrow(validated_targets_export),
                     n_cols = ncol(validated_targets_export))
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

export_table9 <- function() {
  fname <- file.path(DIR_SUPP, "Supplementary_Table_9_miRNA_target_overlap.xlsx")
  wb    <- createWorkbook()
  
  ## ---- helpers --------------------------------------------------------
  get_target_vec <- function(mir_key, type) {
    g <- if (type == "Validated") validated_targets[[mir_key]]$`Target Gene`
    else                     predicted_targets[[mir_key]]$Gene.Symbol
    unique(g[!is.na(g) & nzchar(g)])
  }
  
  classify <- function(genes) {
    is_micro <- toupper(genes) %in% toupper(ALL_MICROGLIAL_MARKERS)
    is_gmrik <- is_gm_rik(genes)
    tibble(
      Is_microglial_marker = is_micro,
      Is_Gm_Rik            = is_gmrik,
      Is_bona_fide         = !is_micro & !is_gmrik
    )
  }
  
  overlap_status <- function(detected, is_bona_fide) {
    dplyr::case_when(
      !detected     ~ "Not detected",
      is_bona_fide  ~ "Retained (bona fide)",
      TRUE          ~ "Excluded (marker/Gm-Rik)"
    )
  }
  
  write_target_sheet <- function(sheet_name, df) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
    apply_sheet_styles(wb, sheet_name, nrow(df), ncol(df))
    excl_rows <- which(df$Overlap_status == "Excluded (marker/Gm-Rik)") + 1
    if (length(excl_rows) > 0)
      addStyle(wb, sheet_name, createStyle(fgFill = "#FFCCCC"),
               rows = excl_rows, cols = 1:ncol(df),
               gridExpand = TRUE, stack = TRUE)
  }
  
  sheet_registry <- list()  # name -> list(desc, count)
  
  ## ---- PROTEOMICS: every target, detected or not ----------------------
  prot_rows <- map_dfr(seq_along(MIRNA_LABELS), function(i) {
    label <- MIRNA_LABELS[i]; key <- MIRNA_KEYS[i]
    map_dfr(c("Female", "Male"), function(sx) {
      sx_key <- ifelse(sx == "Female", "F", "M")
      map_dfr(c("Validated", "Predicted"), function(tt) {
        genes <- get_target_vec(key, tt)
        if (length(genes) == 0) return(NULL)
        detected <- genes %in% genes_protein[[sx_key]]
        cls <- classify(genes)
        tibble(Gene = genes, Layer = "Proteomics", miRNA = label,
               Sex = sx, Target_type = tt, Detected = detected) %>%
          bind_cols(cls) %>%
          mutate(Overlap_status = overlap_status(Detected, Is_bona_fide))
      })
    })
  }) %>%
    mutate(miRNA = factor(miRNA, levels = MIRNA_LABELS)) %>%
    select(Gene, Layer, miRNA, Sex, Target_type, Detected,
           Is_microglial_marker, Is_Gm_Rik, Is_bona_fide, Overlap_status) %>%
    arrange(miRNA, Sex, Target_type, desc(Detected), Overlap_status, Gene)
  
  write_target_sheet("Proteomics_targets", prot_rows)
  sheet_registry[["Proteomics_targets"]] <- list(
    desc  = "All proteomics targets (Table 8 universe) with detection + filter status",
    count = nrow(prot_rows))
  
  ## ---- SPATIAL: one sheet per cell type -------------------------------
  # Detected := target overlaps a Cohen's D > 1 gene in THIS cell type, that sex.
  spat_summary <- list()
  
  for (ct in sort(unique(high_eff_spatial$cell_type))) {
    
    ct_rows <- map_dfr(seq_along(MIRNA_LABELS), function(i) {
      label <- MIRNA_LABELS[i]; key <- MIRNA_KEYS[i]
      map_dfr(c("Female", "Male"), function(sx) {
        sx_key <- ifelse(sx == "Female", "F", "M")
        ct_hits <- high_eff_spatial %>%
          filter(cell_type == ct, sex == sx_key) %>% pull(gene) %>% unique()
        map_dfr(c("Validated", "Predicted"), function(tt) {
          genes <- get_target_vec(key, tt)
          if (length(genes) == 0) return(NULL)
          detected <- genes %in% ct_hits
          cls <- classify(genes)
          tibble(Gene = genes, Layer = "Spatial", Cell_type = ct,
                 miRNA = label, Sex = sx, Target_type = tt, Detected = detected) %>%
            bind_cols(cls) %>%
            mutate(Overlap_status = overlap_status(Detected, Is_bona_fide))
        })
      })
    })
    
    if (is.null(ct_rows) || nrow(ct_rows) == 0) next
    
    ct_rows <- ct_rows %>%
      mutate(miRNA = factor(miRNA, levels = MIRNA_LABELS)) %>%
      select(Gene, Layer, Cell_type, miRNA, Sex, Target_type, Detected,
             Is_microglial_marker, Is_Gm_Rik, Is_bona_fide, Overlap_status) %>%
      arrange(miRNA, Sex, Target_type, desc(Detected), Overlap_status, Gene)
    
    sheet_name <- substr(paste0("Sp_", gsub("[^A-Za-z0-9]", "", ct)), 1, 31)
    write_target_sheet(sheet_name, ct_rows)
    sheet_registry[[sheet_name]] <- list(
      desc  = paste("Spatial |", ct, "| all targets with detection + filter status"),
      count = nrow(ct_rows))
    
    spat_summary[[ct]] <- ct_rows %>%
      group_by(Cell_type, miRNA, Sex, Target_type) %>%
      summarise(N_targets  = n(),
                N_detected = sum(Detected),
                N_retained = sum(Detected & Is_bona_fide),
                N_excluded = sum(Detected & !Is_bona_fide),
                .groups = "drop")
  }
  
  ## ---- Summary --------------------------------------------------------
  summary_prot <- prot_rows %>%
    group_by(miRNA, Sex, Target_type) %>%
    summarise(N_targets = n(), N_detected = sum(Detected),
              N_retained = sum(Detected & Is_bona_fide),
              N_excluded = sum(Detected & !Is_bona_fide), .groups = "drop") %>%
    mutate(Layer = "Proteomics", Cell_type = NA_character_)
  
  summary_spat <- bind_rows(spat_summary) %>% mutate(Layer = "Spatial")
  
  summary_all <- bind_rows(summary_prot, summary_spat) %>%
    select(Layer, Cell_type, miRNA, Sex, Target_type,
           N_targets, N_detected, N_retained, N_excluded) %>%
    arrange(Layer, Cell_type, miRNA, Sex, Target_type)
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_all)
  apply_sheet_styles(wb, "Summary", nrow(summary_all), ncol(summary_all))
  
  # Move Summary to the front
  n_sheets <- length(names(wb))
  worksheetOrder(wb) <- c(n_sheets, seq_len(n_sheets - 1))
  
  saveWorkbook(wb, fname, overwrite = TRUE)
  cat("Saved:", fname, "\n")
}

# Run the functions
#export_table3()
#export_table4()
#export_table5()
#export_table6()
#export_table7()
#export_table8()
#export_table9()


# ---------------------- miRNA-target pathway convergence onto aging / microglia / sex ----------------------
# CONFIG
GO_QVAL   <- 0.05
GO_MIN_GS <- 10
GO_MAX_GS <- 500
MIN_GENES <- 5

# Established spillover / implausible-term exclusions
ENRICH_EXCLUDE <- c("cardiac", "muscle", "Schwann", "sarcolemma",
                    "cytokinesis", "mitotic", "myofibril", "contractile")

# Functional theme keywords (the enrichment layer — what GO can actually return)
THEME_KEYWORDS <- list(
  Synaptic          = c("synap", "neurotransmitter", "axon", "dendrit",
                        "vesicle", "postsynap", "presynap", "glutamat", "GABA"),
  Mito_OXPHOS       = c("mitochond", "oxidative phosphoryl", "respiratory chain",
                        "electron transport", "ATP synth", "aerobic", "NADH"),
  Immune_Microglial = c("immune", "inflammat", "cytokine", "interleukin",
                        "kappaB", "Toll", "defense response", "leukocyte",
                        "myeloid", "phagocyt", "complement"))

# Functional theme -> reviewer axis + evidence provenance
AXIS_CROSSWALK <- tibble::tribble(
  ~Theme,              ~Reviewer_axis,       ~Provenance,
  "Immune_Microglial", "Microglia function", "Enrichment (immune/inflammatory GO terms)",
  "Synaptic",          "Brain aging",        "Enrichment + Berchtold 2013 (direction-neutral)",
  "Mito_OXPHOS",       "Brain aging",        "Enrichment + OXPHOS-in-aging literature"
)

SEX_ANCHORS <- c("Igf1r", "Irak1")

### HELPERS
tag_theme <- function(descriptions) {
  vapply(descriptions, function(d) {
    hit <- names(THEME_KEYWORDS)[vapply(THEME_KEYWORDS, function(kw)
      any(grepl(paste(kw, collapse = "|"), d, ignore.case = TRUE)), logical(1))]
    if (length(hit) == 0) "Other" else paste(hit, collapse = "/")
  }, character(1))
}

# Returns list(result = df-or-NULL, status = one of the three below, n = set size)
run_go_bp <- function(genes, universe, label, exclude = ENRICH_EXCLUDE) {
  genes    <- unique(genes[!is.na(genes) & nzchar(genes)])
  universe <- unique(universe[!is.na(universe) & nzchar(universe)])
  
  if (length(genes) < MIN_GENES) {
    message(sprintf("  %-20s n=%d  -> NOT TESTABLE (< %d genes)", label, length(genes), MIN_GENES))
    return(list(result = NULL, status = "not_testable", n = length(genes)))
  }
  
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_QVAL,
    qvalueCutoff  = GO_QVAL,
    minGSSize     = GO_MIN_GS,
    maxGSSize     = GO_MAX_GS,
    readable      = FALSE
  )
  res <- as.data.frame(ego)
  
  if (nrow(res) == 0) {
    message(sprintf("  %-20s n=%d  -> tested, no enriched terms", label, length(genes)))
    return(list(result = NULL, status = "not_enriched", n = length(genes)))
  }
  
  res <- res[!grepl(paste(exclude, collapse = "|"), res$Description, ignore.case = TRUE), ]
  if (nrow(res) == 0) {
    message(sprintf("  %-20s n=%d  -> only spillover terms, all excluded", label, length(genes)))
    return(list(result = NULL, status = "not_enriched", n = length(genes)))
  }
  
  res <- res %>%
    dplyr::mutate(Set = label, Theme = tag_theme(Description)) %>%
    dplyr::mutate(dplyr::across(c(Set, Theme), as.character)) %>%
    dplyr::select(Set, Theme, ID, Description, GeneRatio, BgRatio,
                  p.adjust, qvalue, Count, geneID) %>%
    dplyr::arrange(p.adjust)
  
  message(sprintf("  %-20s n=%d  -> %d terms retained", label, length(genes), nrow(res)))
  list(result = res, status = "enriched", n = length(genes))
}

### BUILD PER-SEX TARGET-HIT SETS
## PROTEOMICS
protein_target_hits <- protein_overlaps %>%
  dplyr::select(Sex, hits_filtered) %>%
  tidyr::unnest(hits_filtered) %>%
  dplyr::distinct(Sex, gene = hits_filtered)

protein_universe <- unique(perseus_proc_plx_cohen$Genes)

protein_sets <- list(
  Female = protein_target_hits %>% dplyr::filter(Sex == "Female") %>% dplyr::pull(gene),
  Male   = protein_target_hits %>% dplyr::filter(Sex == "Male")   %>% dplyr::pull(gene)
)

## SPATIAL (pooled across cell types, per sex)
spatial_target_hits <- all_overlaps_spatial_filtered %>%
  dplyr::select(Sex, hits_filtered) %>%
  tidyr::unnest(hits_filtered) %>%
  dplyr::distinct(Sex, gene = hits_filtered)

spatial_universe <- unique(all_results_spatial$gene)

spatial_sets <- list(
  Female = spatial_target_hits %>% dplyr::filter(Sex == "Female") %>% dplyr::pull(gene),
  Male   = spatial_target_hits %>% dplyr::filter(Sex == "Male")   %>% dplyr::pull(gene)
)

cat("\n=== Target-hit set sizes (union across 5 miRNAs, filtered) ===\n")
cat(sprintf("Proteomics  F: %3d  M: %3d\n", length(protein_sets$Female), length(protein_sets$Male)))
cat(sprintf("Spatial     F: %3d  M: %3d\n", length(spatial_sets$Female), length(spatial_sets$Male)))

### RUN ENRICHMENT
cat("\n=== Running GO:BP enrichment (detected-feature background) ===\n")
runs <- list(
  Proteomics_Female = run_go_bp(protein_sets$Female, protein_universe, "Proteomics_Female"),
  Proteomics_Male   = run_go_bp(protein_sets$Male,   protein_universe, "Proteomics_Male"),
  Spatial_Female    = run_go_bp(spatial_sets$Female, spatial_universe, "Spatial_Female"),
  Spatial_Male      = run_go_bp(spatial_sets$Male,   spatial_universe, "Spatial_Male")
)

enrich_results <- dplyr::bind_rows(lapply(runs, `[[`, "result"))

# Status ledger — distinguishes "tested, nothing" from "couldn't test"
set_status <- tibble::tibble(
  Set    = names(runs),
  n      = vapply(runs, `[[`, integer(1),   "n"),
  status = vapply(runs, `[[`, character(1), "status")
) %>%
  tidyr::separate(Set, into = c("Layer", "Sex"), sep = "_", remove = FALSE)

cat("\n=== Enrichment status per set ===\n")
print(set_status)

### CONVERGENCE SUMMARY (functional themes)
if (nrow(enrich_results) > 0) {
  
  convergence_summary <- enrich_results %>%
    dplyr::filter(Theme != "Other") %>%
    dplyr::count(Set, Theme, name = "n_terms") %>%
    tidyr::pivot_wider(names_from = Theme, values_from = n_terms, values_fill = 0)
  
  cat("\n=== Pathway convergence by functional theme (enriched BP terms per set) ===\n")
  print(convergence_summary)
  
  axis_from_enrichment <- enrich_results %>%
    dplyr::filter(Theme != "Other") %>%
    tidyr::separate_rows(Theme, sep = "/") %>%
    dplyr::left_join(AXIS_CROSSWALK, by = "Theme") %>%
    dplyr::filter(!is.na(Reviewer_axis)) %>%
    tidyr::separate(Set, into = c("Layer", "Sex"), sep = "_", remove = FALSE)
  
  # Axes 1-2: microglia function, brain aging (from enrichment)
  axis_summary <- axis_from_enrichment %>%
    dplyr::count(Layer, Sex, Reviewer_axis, name = "n_terms") %>%
    tidyr::pivot_wider(names_from = Reviewer_axis, values_from = n_terms, values_fill = 0)
  
  print(axis_summary)
  
  # Axis 3: sex differences = male-vs-female contrast, joined to status so a
  sex_asymmetry <- axis_from_enrichment %>%
    dplyr::count(Layer, Reviewer_axis, Sex, name = "n_terms") %>%
    tidyr::pivot_wider(names_from = Sex, values_from = n_terms, values_fill = 0)
  
  # ensure both sex columns exist even if one sex enriched nothing anywhere
  for (col in c("Female", "Male")) if (!col %in% names(sex_asymmetry)) sex_asymmetry[[col]] <- 0L
  
  sex_asymmetry <- sex_asymmetry %>%
    dplyr::mutate(
      Male_minus_Female = Male - Female,
      Female_status = set_status$status[match(paste0(Layer, "_Female"), set_status$Set)],
      Male_status   = set_status$status[match(paste0(Layer, "_Male"),   set_status$Set)]
    )
  
  print(sex_asymmetry)
  
  # Top themed terms for the response letter / figure
  top_themed <- enrich_results %>%
    dplyr::filter(Theme != "Other") %>%
    dplyr::group_by(Set, Theme) %>%
    dplyr::slice_min(p.adjust, n = 5, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  cat("\n=== Top themed terms per set (for rebuttal / figure) ===\n")
  print(top_themed %>% dplyr::select(Set, Theme, Description, p.adjust, Count), n = Inf)
  
} else {
  cat("\n!! No enriched terms in any set — nothing to summarise.\n")
  cat("   Check set_status above: are the sets 'not_testable' (n<5) or 'not_enriched'?\n")
}


# ---------------------- DIFFERENTIAL ABUNDANCE ANALYSIS ----------------------
### Reshape to one row per sample
prop_wide <- cell_number_spatial %>%
  filter(!anno_0825 %in% CELL_TYPES_EXCLUDE_SPATIAL) %>%
  group_by(Sample, Group, Sex, anno_0825) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  pivot_wider(names_from  = anno_0825,
              values_from = mean_prop,
              values_fill = 0) %>%
  mutate(
    Group = factor(Group, levels = c("Control", "PLX")),
    Sex   = factor(Sex,   levels = c("F", "M"))
  )

cell_cols <- setdiff(colnames(prop_wide), c("Sample", "Group", "Sex"))
Y <- DR_data(prop_wide[, cell_cols])

# Glutamatergic = reference
base_idx <- which(cell_cols == "Glutamatergic")
fit_alt <- DirichReg(Y ~ Group + Sex | Group + Sex,
                     data  = prop_wide,
                     model = "alternative",
                     base  = base_idx)


### Extract GroupPLX mean-model effects + SEs robustly
beta <- coef(fit_alt)$beta                    
beta <- beta[!sapply(beta, is.null)]          
ses <- sqrt(diag(vcov(fit_alt)))          

da <- tibble(
  Cell_Type = names(beta),
  estimate  = vapply(beta, function(v) unname(v["GroupPLX"]), numeric(1)),
  se        = unname(ses[paste0(names(beta), ":GroupPLX")])
) %>%
  bind_rows(tibble(Cell_Type = "Glutamatergic", estimate = 0, se = NA_real_)) %>%
  mutate(
    lo   = estimate - 1.96 * se,
    hi   = estimate + 1.96 * se,
    p    = 2 * pnorm(-abs(estimate / se)),
    p_BH = p.adjust(p, method = "BH"),
    sig  = case_when(
      Cell_Type == "Glutamatergic" ~ "reference",
      !is.na(p_BH) & p_BH < 0.05    ~ "PLX effect",
      TRUE                          ~ "Non-significant"
    ),
    Cell_Type = fct_reorder(Cell_Type, estimate)
  )

print(da)

### Plot the differential abundance analysis
diff_abun_analysis <- ggplot(da, aes(estimate, Cell_Type, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c(
      "PLX effect" = "#E8836B",
      "Non-significant" = "grey50",
      "reference" = "grey75"),
    breaks = c("PLX effect", "Non-significant", "reference"),
    labels = c("PLX effect", "Non-significant", "Reference"),
    name = NULL) +
  labs(x = "PLX effect (log-ratio vs. glutamatergic reference)", y = NULL, title = "Compositional Differential Abundance",
       subtitle = "Treatment effect controlling for sex") +
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
    plot.subtitle = element_text(size = 15,color = "grey35",hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 12)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank())

print(diff_abun_analysis)


### Save
#ggsave("Figures and plots/PNGs/FigR1_dirichlet_DA.png", diff_abun_analysis, width = 10, height = 6, dpi = 300)
#ggsave("Figures and plots/PDFs/FigR1_dirichlet_DA.pdf", diff_abun_analysis, width = 10, height = 6)


# ---------------------- VARIANCE PARTITION — SPATIAL PSEUDOBULK ----------------------
# Expression matrix: genes x pseudobulk samples
# Filter metadata to the six main cell types used throughout
spatial_meta_vp <- metadata %>%
  filter(anno_0825 %in% CELL_TYPES) %>%
  mutate(
    Animal    = factor(Sample),
    CellType  = factor(recode(anno_0825, "Oligodendyocyte" = "Oligodendrocyte")),
    Treatment = factor(Group, levels = c("Control", "PLX")),
    Sex       = factor(Sex,   levels = c("F", "M"))
  )

# Build genes x samples matrix from pseudobulk counts
spatial_lcpm <- counts %>%
  filter(row_id %in% spatial_meta_vp$row_id) %>%
  column_to_rownames("row_id") %>%
  as.matrix() %>%
  t()

# log-CPM normalization
lib_sizes       <- colSums(spatial_lcpm)
spatial_lcpm    <- log2(t(t(spatial_lcpm) / lib_sizes * 1e6) + 1)

# Filter: keep genes with mean log-CPM > 1 and non-zero variance
expressed_spatial <- rowMeans(spatial_lcpm) > 1 &
  matrixStats::rowVars(spatial_lcpm) > 0
spatial_lcpm      <- spatial_lcpm[expressed_spatial, , drop = FALSE]

# Enforce sample order match between matrix and metadata
spatial_meta_vp        <- spatial_meta_vp[match(colnames(spatial_lcpm),
                                                spatial_meta_vp$row_id), ]
rownames(spatial_meta_vp) <- spatial_meta_vp$row_id
stopifnot(identical(colnames(spatial_lcpm), rownames(spatial_meta_vp)))

cat(sprintf("Spatial variance partition: %d genes x %d pseudobulk samples\n",
            nrow(spatial_lcpm), ncol(spatial_lcpm)))

# Audit: each animal should contribute ~6 cell-type pseudobulks
cat("\nAnimal x CellType table (should be ~1 per cell):\n")
print(table(spatial_meta_vp$Animal, spatial_meta_vp$CellType))

# Global partition
formula_spatial_all  <- ~ (1|CellType) + (1|Sex) + (1|Treatment) + (1|Animal)
vpart_spatial_all    <- fitExtractVarPartModel(spatial_lcpm,
                                               formula_spatial_all,
                                               spatial_meta_vp)
vpart_spatial_all    <- sortCols(vpart_spatial_all)

cat("\n=== Spatial global variance partition (median %) ===\n")
print(round(sort(colMeans(vpart_spatial_all) * 100, decreasing = TRUE), 2))

VPART_COLORS_SPATIAL <- c(
  CellType  = "#66b2ff",
  Animal    = "#538f8f",
  Treatment = "#8b3a7a",
  Sex       = "#a42e30",
  Residuals = "#d3d3d3"
)

VPART_THEME <- theme_minimal(base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 14, face = "bold",   hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.text        = element_text(size = 12, color = "black"),
    axis.title.y     = element_text(size = 13, face = "bold"),
    panel.border     = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

p_vpart_spatial_all <- plotVarPart(vpart_spatial_all,
                                   col = VPART_COLORS_SPATIAL[colnames(vpart_spatial_all)]) +
  labs(title    = "Spatial pseudobulk variance partition — all cell types",
       y        = "Variance explained (%)") +
  VPART_THEME

print(p_vpart_spatial_all)

### VARIANCE PARTITION — PROTEOMICS
# NOTE: In proteomics each animal IS a sample
# Expression matrix
prot_mat <- perseus_proc_plx %>%
  select(Genes, all_of(ALL_SAMPLE_COLS)) %>%
  filter(!is.na(Genes), nzchar(Genes)) %>%
  distinct(Genes, .keep_all = TRUE) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# Keep only proteins with full observations and non-zero variance
fully_observed <- rowSums(!is.na(prot_mat)) == ncol(prot_mat)
variable_prot  <- matrixStats::rowVars(prot_mat, na.rm = TRUE) > 0
prot_mat       <- prot_mat[fully_observed & variable_prot, , drop = FALSE]

cat(sprintf("\nProteomics variance partition: %d proteins x %d samples\n",
            nrow(prot_mat), ncol(prot_mat)))

# Sample metadata
prot_sample_info <- data.frame(
  Sample    = ALL_SAMPLE_COLS,
  Treatment = factor(ifelse(grepl("^PLX",     ALL_SAMPLE_COLS), "PLX", "Control"),
                     levels = c("Control", "PLX")),
  Sex       = factor(ifelse(grepl("_F_", ALL_SAMPLE_COLS), "F", "M"),
                     levels = c("F", "M")),
  row.names = ALL_SAMPLE_COLS
)
prot_sample_info <- prot_sample_info[colnames(prot_mat), ]
stopifnot(identical(rownames(prot_sample_info), colnames(prot_mat)))

cat("\nGroup sizes (proteomics):\n")
print(table(prot_sample_info$Sex, prot_sample_info$Treatment))

### Partition: Sex + Treatment + Residual
formula_prot  <- ~ (1|Sex) + (1|Treatment)
vpart_prot    <- fitExtractVarPartModel(prot_mat, formula_prot, prot_sample_info)
vpart_prot    <- sortCols(vpart_prot)

cat("\n=== Proteomics variance partition (median %) ===\n")
print(round(sort(colMeans(vpart_prot) * 100, decreasing = TRUE), 2))

VPART_COLORS_PROT <- c(
  Treatment = "#8b3a7a",
  Sex       = "#a42e30",
  Residuals = "#d3d3d3"
)

p_vpart_prot <- plotVarPart(vpart_prot,
                            col = VPART_COLORS_PROT[colnames(vpart_prot)]) +
  labs(title    = "Proteomics variance partition",
       y        = "Variance explained (%)") +
  VPART_THEME

print(p_vpart_prot)

### Save
#ggsave(file.path(DIR_PNG, "VarPart_Spatial_global.png"), p_vpart_spatial_all, width = 8, height = 6, dpi = 300)
#ggsave(file.path(DIR_PDF, "VarPart_Spatial_global.pdf"), p_vpart_spatial_all, width = 8, height = 6)

#ggsave(file.path(DIR_PNG, "VarPart_Proteomics.png"), p_vpart_prot, width = 8, height = 6, dpi = 300)
#ggsave(file.path(DIR_PDF, "VarPart_Proteomics.pdf"), p_vpart_prot, width = 8, height = 6)


