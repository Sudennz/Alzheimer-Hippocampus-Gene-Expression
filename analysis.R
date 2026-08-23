# =============================================================
# Alzheimer's Disease Gene Expression Analysis - Hippocampus
# Dataset: GSE5281 (NCBI GEO)
# =============================================================

# ---- 1. Load required packages ----
library(GEOquery)
library(hgu133plus2.db)
library(randomForest)

# ---- 2. Download the dataset ----
gse <- getGEO("GSE5281", GSEMatrix = TRUE)
gse <- gse[[1]]

# ---- 3. Clean phenotype (metadata) columns ----
pheno <- pData(gse)

pheno$disease_status <- ifelse(
  !is.na(pheno$"disease state:ch1"),
  pheno$"disease state:ch1",
  pheno$"Disease State:ch1"
)

pheno$organ_region <- as.character(ifelse(
  !is.na(as.character(pheno$"organ region:ch1")),
  as.character(pheno$"organ region:ch1"),
  as.character(pheno$"Organ Region:ch1")
))
pheno$organ_region <- trimws(pheno$organ_region)
pheno$organ_region[grepl("cingulate|singulate", pheno$organ_region, ignore.case = TRUE)] <- "Posterior Cingulate"

# ---- 4. Select hippocampus samples only ----
hippo_samples <- pheno$geo_accession[grepl("^hippocampus", pheno$organ_region, ignore.case = TRUE)]

# ---- 5. Extract and filter expression matrix ----
expr_matrix <- exprs(gse)
expr_hippo <- expr_matrix[, hippo_samples]

# ---- 6. Log2 transform ----
expr_hippo_log <- log2(expr_hippo + 1)

# ---- 7. Remove AFFX control probes ----
control_probes <- grepl("^AFFX", rownames(expr_hippo_log))
expr_hippo_clean <- expr_hippo_log[!control_probes, ]

# ---- 8. Group labels (AD vs Normal) ----
group <- ifelse(
  startsWith(as.character(pheno$disease_status[match(colnames(expr_hippo_clean), pheno$geo_accession)]), "A"),
  "AD", "Normal"
)

# ---- 9. PCA ----
pca_input <- t(expr_hippo_clean)
pca_result <- prcomp(pca_input, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$disease_status <- pheno$disease_status[match(rownames(pca_scores), pheno$geo_accession)]

png("pca_plot.png", width = 800, height = 600)
plot(pca_scores$PC1, pca_scores$PC2,
     col = ifelse(startsWith(as.character(pca_scores$disease_status), "A"), "red", "blue"),
     pch = 19,
     xlab = "PC1 (18.5%)",
     ylab = "PC2 (8.1%)",
     main = "Hippocampus PCA: AD vs Normal")
legend("topright", legend = c("Alzheimer's Disease", "Normal"),
       col = c("red", "blue"), pch = 19)
dev.off()

# ---- 10. Differential expression (t-test per gene) ----
t_test_clean <- apply(expr_hippo_clean, 1, function(gene_values) {
  ad_values <- gene_values[group == "AD"]
  normal_values <- gene_values[group == "Normal"]
  test <- t.test(ad_values, normal_values)
  c(mean_AD = mean(ad_values), mean_Normal = mean(normal_values),
    log2FC = mean(ad_values) - mean(normal_values), p_value = test$p.value)
})
t_test_clean <- as.data.frame(t(t_test_clean))
t_test_clean$p_adj <- p.adjust(t_test_clean$p_value, method = "BH")

png("volcano_plot.png", width = 800, height = 600)
plot(t_test_clean$log2FC, -log10(t_test_clean$p_value),
     col = ifelse(t_test_clean$log2FC > 1 & t_test_clean$p_adj < 0.05, "red",
                  ifelse(t_test_clean$log2FC < -1 & t_test_clean$p_adj < 0.05, "blue", "gray")),
     pch = 19, cex = 0.6,
     xlab = "log2 Fold Change (AD vs Normal)",
     ylab = "-log10(p-value)",
     main = "Volcano Plot: Hippocampus AD vs Normal")
abline(v = c(-1, 1), lty = 2, col = "black")
legend("topright", legend = c("Up in AD", "Down in AD", "Not Significant"),
       col = c("red", "blue", "gray"), pch = 19)
dev.off()

# ---- 11. Annotate top genes with real gene symbols ----
top_clean <- t_test_clean[order(t_test_clean$p_adj), ]
top_clean$PROBEID <- rownames(top_clean)

probe_to_symbol <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(top_clean),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

top_annotated <- merge(top_clean, probe_to_symbol, by = "PROBEID")
top_annotated <- top_annotated[order(top_annotated$p_adj), ]

# ---- 12. Classification models (LOOCV) ----
top10_probes <- head(top_annotated$PROBEID, 10)
model_data <- as.data.frame(t(expr_hippo_clean[top10_probes, ]))
model_data$status <- as.factor(group)

# Logistic Regression
predictions_glm <- character(nrow(model_data))
for (i in 1:nrow(model_data)) {
  train_data <- model_data[-i, ]
  test_data <- model_data[i, ]
  model <- glm(status ~ ., data = train_data, family = "binomial")
  prob <- predict(model, newdata = test_data, type = "response")
  predictions_glm[i] <- ifelse(prob > 0.5, levels(model_data$status)[2], levels(model_data$status)[1])
}
table(Predicted = predictions_glm, Actual = model_data$status)

# Random Forest
predictions_rf <- character(nrow(model_data))
set.seed(42)
for (i in 1:nrow(model_data)) {
  train_x <- model_data[-i, setdiff(colnames(model_data), "status")]
  train_y <- model_data$status[-i]
  test_x <- model_data[i, setdiff(colnames(model_data), "status"), drop = FALSE]
  rf_model <- randomForest(x = train_x, y = train_y, ntree = 500)
  pred <- predict(rf_model, newdata = test_x)
  predictions_rf[i] <- as.character(pred)
}
table(Predicted = predictions_rf, Actual = model_data$status)