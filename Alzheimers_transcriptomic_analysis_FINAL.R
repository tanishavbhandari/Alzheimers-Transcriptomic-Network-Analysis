# ============================================================
# Alzheimer's Disease Transcriptomic and Network Analysis
# Dataset: GSE1297
# Comparison: Severe AD vs Control hippocampal samples
# ============================================================

# This script reconstructs the R-based part of the analysis
# from the original R history.
#
# Main steps:
# 1. Download GSE1297
# 2. Select Control and Severe AD samples
# 3. Log2-transform expression values
# 4. Differential expression analysis with limma
# 5. Create volcano plot
# 6. Select top 150 probes by nominal P-value
# 7. Annotate probes using hgu133a.db
# 8. Create the 143-gene exploratory list
# 9. Export the results for Enrichr and STRING
#
# GO/KEGG enrichment and STRING network analysis in the report
# were performed using the Enrichr and STRING web interfaces.
# Their exported result tables can be added later if desired.

# ============================================================
# 1. INSTALL / LOAD PACKAGES
# ============================================================

# Run these installation lines only if the packages are not
# already installed on your computer.

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery", update = FALSE, ask = FALSE)

if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)

if (!requireNamespace("hgu133a.db", quietly = TRUE))
    BiocManager::install("hgu133a.db", update = FALSE, ask = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

library(GEOquery)
library(limma)
library(hgu133a.db)
library(ggplot2)

# ============================================================
# 2. CREATE PROJECT FOLDERS
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# 3. DOWNLOAD GSE1297
# ============================================================

# GSE1297 contains hippocampal gene-expression data from
# control and Alzheimer's disease samples.

gse <- getGEO("GSE1297", GSEMatrix = TRUE)

# GSE1297 is expected to contain one expression-set object.
# Select the first object.

if (is.list(gse)) {
    gse <- gse[[1]]
}

# ============================================================
# 4. INSPECT SAMPLE INFORMATION
# ============================================================

pd <- pData(gse)

# The original analysis used the "group:ch1" field.
pd$group_simple <- pd$`group:ch1`

print(unique(pd$group_simple))

# ============================================================
# 5. SELECT CONTROL AND SEVERE AD SAMPLES
# ============================================================

keep_samples <- pd$group_simple %in% c("Control", "Severe")

gse_sub <- gse[, keep_samples]

# Expression matrix
exprs_data <- exprs(gse_sub)

# Define the experimental groups.
# Control is the reference group.
group <- factor(
    pData(gse_sub)$`group:ch1`,
    levels = c("Control", "Severe")
)

print(table(group))

# Expected:
# Control = 9
# Severe  = 7

# ============================================================
# 6. CHECK EXPRESSION SCALE
# ============================================================

print(range(exprs_data))

# Plot original expression distribution
png(
    "figures/Expression_distribution_before_log2.png",
    width = 2000,
    height = 1500,
    res = 250
)

boxplot(
    exprs_data,
    main = "Expression Value Distribution per Sample",
    las = 2,
    cex.axis = 0.7
)

dev.off()

# ============================================================
# 7. LOG2 TRANSFORMATION
# ============================================================

# The original analysis used:
# log2(expression + 1)

exprs_data_log <- log2(exprs_data + 1)

print(range(exprs_data_log))

# Plot log2-transformed expression distribution
png(
    "figures/Expression_distribution_log2.png",
    width = 2000,
    height = 1500,
    res = 250
)

boxplot(
    exprs_data_log,
    main = "Log2-Transformed Expression Distribution",
    las = 2,
    cex.axis = 0.7
)

dev.off()

# ============================================================
# 8. DESIGN MATRIX
# ============================================================

# Build the group vector directly from the selected samples and
# explicitly verify that it has one value per expression column.
selected_pd <- pData(gse_sub)
selected_group <- selected_pd$`group:ch1`

group <- factor(
    selected_group,
    levels = c("Control", "Severe")
)

if (length(group) != ncol(exprs_data_log)) {
    stop(
        paste0(
            "Group length (", length(group),
            ") does not match expression samples (",
            ncol(exprs_data_log), ")."
        )
    )
}

if (any(is.na(group))) {
    stop("Some selected samples have missing/unknown group labels.")
}

# Standard limma design: Control is the reference group.
design <- model.matrix(~ group)

if (nrow(design) != ncol(exprs_data_log)) {
    stop(
        paste0(
            "Design rows (", nrow(design),
            ") do not match expression columns (",
            ncol(exprs_data_log), ")."
        )
    )
}

colnames(design) <- c("Intercept", "Severe_vs_Control")

print(table(group))
print(dim(exprs_data_log))
print(dim(design))
print(design)

# ============================================================
# 9. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================

fit <- lmFit(exprs_data_log, design)

fit <- eBayes(fit)

results <- topTable(
    fit,
    coef = "Severe_vs_Control",
    number = Inf,
    sort.by = "P"
)

# Save complete differential-expression results
write.csv(
    results,
    "results/All_differential_expression_results.csv",
    row.names = TRUE
)

# ============================================================
# 10. BASIC DIFFERENTIAL EXPRESSION SUMMARY
# ============================================================

number_nominal <- sum(results$P.Value < 0.05)

number_fdr <- sum(results$adj.P.Val < 0.05)

number_fdr_up <- sum(
    results$adj.P.Val < 0.05 &
    results$logFC > 0
)

number_fdr_down <- sum(
    results$adj.P.Val < 0.05 &
    results$logFC < 0
)

cat("\n--------------------------------------------\n")
cat("Differential Expression Summary\n")
cat("--------------------------------------------\n")
cat("Total probes:", nrow(results), "\n")
cat("Nominal P < 0.05:", number_nominal, "\n")
cat("FDR-adjusted P < 0.05:", number_fdr, "\n")
cat("FDR-significant upregulated:", number_fdr_up, "\n")
cat("FDR-significant downregulated:", number_fdr_down, "\n")
cat("--------------------------------------------\n\n")

# ============================================================
# 11. VOLCANO PLOT
# ============================================================

results$significant <- ifelse(
    results$P.Value < 0.05,
    "Nominal P < 0.05",
    "Not nominally significant"
)

volcano <- ggplot(
    results,
    aes(
        x = logFC,
        y = -log10(P.Value),
        color = significant
    )
) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
        values = c(
            "Nominal P < 0.05" = "red",
            "Not nominally significant" = "grey"
        )
    ) +
    theme_classic(base_size = 13) +
    labs(
        title = "Volcano Plot: Severe AD vs Control",
        x = "Log2 Fold Change",
        y = "-log10(P-value)"
    )

print(volcano)

ggsave(
    "figures/Figure_1_Volcano_Plot.png",
    plot = volcano,
    width = 9,
    height = 7,
    dpi = 300
)

# ============================================================
# 12. TOP 150 PROBES
# ============================================================

top_genes <- results[
    order(results$P.Value),
]

top_genes <- head(top_genes, 150)

# Save the top 150 probes
write.csv(
    top_genes,
    "results/Top_150_probes.csv",
    row.names = TRUE
)

probe_ids <- rownames(top_genes)

cat("Top 150 probes selected:", length(probe_ids), "\n")

# ============================================================
# 13. PROBE ANNOTATION
# ============================================================

gene_map <- select(
    hgu133a.db,
    keys = probe_ids,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "PROBEID"
)

# Remove probes without a gene symbol
gene_map_clean <- gene_map[
    !is.na(gene_map$SYMBOL),
]

# Keep one entry per gene symbol
gene_map_clean <- gene_map_clean[
    !duplicated(gene_map_clean$SYMBOL),
]

# Final gene list
gene_list_final <- gene_map_clean$SYMBOL

cat(
    "Unique annotated genes:",
    length(gene_list_final),
    "\n"
)

# Save annotation table
write.csv(
    gene_map_clean,
    "results/Top_150_probe_annotation.csv",
    row.names = FALSE
)

# Save the gene list for Enrichr and STRING
writeLines(
    gene_list_final,
    "results/genes_for_Enrichr_and_STRING.txt"
)

# ============================================================
# 14. TOP 20 DIFFERENTIAL EXPRESSION TABLE
# ============================================================

table1 <- data.frame(
    Probe_ID = rownames(results)[1:20],
    logFC = results$logFC[1:20],
    AveExpr = results$AveExpr[1:20],
    t = results$t[1:20],
    P_value = results$P.Value[1:20],
    Adjusted_P_value = results$adj.P.Val[1:20]
)

write.csv(
    table1,
    "results/Table_1_Top_20_Differential_Expression.csv",
    row.names = FALSE
)

print(table1)

# ============================================================
# 15. SAVE R ANALYSIS OBJECTS
# ============================================================

save(
    gse_sub,
    pd,
    exprs_data_log,
    group,
    design,
    fit,
    results,
    gene_map,
    gene_map_clean,
    gene_list_final,
    file = "results/analysis_objects.RData"
)

# ============================================================
# 16. SOFTWARE VERSIONS
# ============================================================

version_info <- c(
    R = R.version.string,
    limma = as.character(packageVersion("limma")),
    GEOquery = as.character(packageVersion("GEOquery")),
    hgu133a.db = as.character(packageVersion("hgu133a.db")),
    ggplot2 = as.character(packageVersion("ggplot2"))
)

write.table(
    version_info,
    "results/software_versions.txt",
    quote = FALSE,
    col.names = FALSE
)

# ============================================================
# 17. FINAL MESSAGE
# ============================================================

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")
cat("Results saved in the 'results' folder.\n")
cat("Figures saved in the 'figures' folder.\n")
cat("The gene list is ready for Enrichr and STRING.\n")
cat("============================================\n")
