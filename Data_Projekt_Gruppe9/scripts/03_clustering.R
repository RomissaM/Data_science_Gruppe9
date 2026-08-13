# ------------------------------------------------------------
# 1. Normalisierte Daten laden
# ------------------------------------------------------------

library(Biobase)

MLL.A_rma <- readRDS("../data/processed/MLL_A_rma.rds")
# Expressionsmatrix extrahieren
expr_matrix <- exprs(MLL.A_rma)

dim(expr_matrix)
# Für PCA Sample
expr_matrix_t <- t(expr_matrix)

dim(expr_matrix_t)
# ------------------------------------------------------------
# 3. PCA durchführen
# ------------------------------------------------------------

pca_result <- prcomp(
  expr_matrix_t,
  center = TRUE,
  scale. = FALSE
)
summary(pca_result)

# ------------------------------------------------------------
# 4. PCA-Plot
# ------------------------------------------------------------

plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  xlab = "PC1 (21.67 %)",
  ylab = "PC2 (14.81 %)",
  main = "PCA der normalisierten Expressionsdaten",
  pch = 19
)
text(
  pca_result$x[, 1],
  pca_result$x[, 2],
  labels = 1:20,
  pos = 3,
  cex = 0.7
)
# Die PCA zeigt Unterschiede in den globalen Expressionsprofilen
# der 20 Proben.
# PC1 erklärt 21.67 % und PC2 14.81 % der Gesamtvarianz.
# Zusammen erklären die ersten beiden Hauptkomponenten 36.48 %.
#
# Im PCA-Plot ist eine Strukturierung der Proben erkennbar.
# Einige Proben liegen deutlich vom zentralen Bereich entfernt.
# Eine biologische Interpretation dieser Gruppen ist anhand der
# PCA allein jedoch nicht möglich.

# ------------------------------------------------------------
# 5. Hierarchisches Clustering
# ------------------------------------------------------------

# Euklidische Distanz zwischen den Proben
sample_dist <- dist(expr_matrix_t)

# Hierarchisches Clustering
hc <- hclust(sample_dist, method = "complete")

# Dendrogramm
plot(
  hc,
  main = "Hierarchisches Clustering der Proben",
  xlab = "Proben",
  sub = "",
  cex = 0.7
)

# Das hierarchische Clustering zeigt eine erkennbare Struktur
# innerhalb der 20 Proben.
# Proben mit ähnlichen globalen Expressionsprofilen werden
# bei geringeren Distanzen zusammengeführt.
#
# Es sind größere sowie mehrere kleinere Gruppen erkennbar.
# Dies stimmt grundsätzlich mit der in der PCA beobachteten
# Heterogenität der Expressionsprofile überein.
#
# Eine biologische Zuordnung der Cluster ist anhand der
# Expressionsdaten allein jedoch nicht möglich.
# ------------------------------------------------------------
# 6. Heatmap der variabelsten Probe-Sets
# ------------------------------------------------------------

# Varianz jedes Probe-Sets berechnen
gene_var <- apply(expr_matrix, 1, var)

# 50 Probe-Sets mit der höchsten Varianz auswählen
top50 <- order(gene_var, decreasing = TRUE)[1:50]

# Expressionswerte dieser Probe-Sets
heatmap_data <- expr_matrix[top50, ]

dim(heatmap_data)

# Heatmap erstellen
heatmap(
  heatmap_data,
  scale = "row",
  Colv = TRUE,
  main = "Heatmap der 50 variabelsten Probe-Sets",
  xlab = "Proben",
  ylab = "Probe-Sets",
  cexCol = 0.7,
  cexRow = 0.4
)
# ------------------------------------------------------------
# 7. Biologische Gruppen der Proben untersuchen
# ------------------------------------------------------------

pheno <- pData(MLL.A_rma)

dim(pheno)
colnames(pheno)
pheno
# ============================================================
# Zusammenfassung der explorativen Analyse
# ============================================================

# PCA und hierarchisches Clustering zeigen eine heterogene Struktur
# der 20 Proben. Die Proben unterscheiden sich in ihren globalen
# Expressionsprofilen.

# Die Heatmap der 50 variabelsten Probe-Sets zeigt ebenfalls
# unterschiedliche Expressionsmuster zwischen den Proben.

# Eine biologische Zuordnung der beobachteten Cluster ist anhand
# der im ExpressionSet enthaltenen Metadaten derzeit nicht möglich.
# ============================================================