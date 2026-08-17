# ------------------------------------------------------------
# 1. Confusion Matrix der äußeren Cross-Validation
# ------------------------------------------------------------

# Confusion Matrix als Data Frame
cm_df <- as.data.frame(cv_cm)

cm_df

# ------------------------------------------------------------
# 2. Confusion Matrix visualisieren
# ------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

ggplot(
  cm_df,
  aes(
    x = Vorhergesagt,
    y = Tatsächlich,
    fill = Freq
  )
) +
  geom_tile() +
  geom_text(
    aes(label = Freq),
    size = 6
  ) +
  labs(
    title = "Confusion Matrix der äußeren 5-fachen Cross-Validation",
    x = "Vorhergesagte Klasse",
    y = "Tatsächliche Klasse"
  ) +
  theme_minimal()
# ------------------------------------------------------------
# 3. ROC-Kurve der äußeren Cross-Validation
# ------------------------------------------------------------

library(pROC)

target_num <- ifelse(target == "T", 1, 0)

roc_cv <- roc(
  response = target_num,
  predictor = cv_probabilities
)

auc_cv <- auc(roc_cv)

auc_cv
# ------------------------------------------------------------
# 4. ROC-Kurve visualisieren
# ------------------------------------------------------------

plot(
  roc_cv,
  main = "ROC-Kurve der äußeren 5-fachen Cross-Validation",
  legacy.axes = TRUE
)

# AUC in die Grafik schreiben
text(
  x = 0.6,
  y = 0.2,
  labels = paste("AUC =", round(as.numeric(auc_cv), 3))
)
# ------------------------------------------------------------
# 5. Ausgewählte LASSO-Gene visualisieren
# ------------------------------------------------------------

# Reihenfolge nach Koeffizient
selected_genes_annotated$SYMBOL <- factor(
  selected_genes_annotated$SYMBOL,
  levels = selected_genes_annotated$SYMBOL[
    order(selected_genes_annotated$Koeffizient)
  ]
)

ggplot(
  selected_genes_annotated,
  aes(
    x = SYMBOL,
    y = Koeffizient,
    fill = Koeffizient > 0
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Vom LASSO-Modell ausgewählte Gene",
    x = "Gen",
    y = "LASSO-Koeffizient"
  ) +
  theme_minimal() +
  guides(fill = "none")
# ------------------------------------------------------------
# 6. Expression ausgewählter Gene nach B/T-Gruppe
# ------------------------------------------------------------

# Probe-ID von CD3D
cd3d_probe <- "38319_at"

# Expressionswerte aus dem ALL-Datensatz
cd3d_expr <- exprs(ALL)[cd3d_probe, ]

# Data Frame für ggplot
cd3d_df <- data.frame(
  Expression = as.numeric(cd3d_expr),
  Klasse = target
)

# Boxplot
ggplot(
  cd3d_df,
  aes(
    x = Klasse,
    y = Expression,
    fill = Klasse
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.6
  ) +
  labs(
    title = "Expression von CD3D in B- und T-Zell-ALL",
    x = "Klasse",
    y = "Expressionswert"
  ) +
  theme_minimal() +
  guides(fill = "none")
# ------------------------------------------------------------
# 7. Expression von CD74
# ------------------------------------------------------------

cd74_probe <- "35016_at"

cd74_expr <- exprs(ALL)[cd74_probe, ]

cd74_df <- data.frame(
  Expression = as.numeric(cd74_expr),
  Klasse = target
)

ggplot(
  cd74_df,
  aes(
    x = Klasse,
    y = Expression,
    fill = Klasse
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.6
  ) +
  labs(
    title = "Expression von CD74 in B- und T-Zell-ALL",
    x = "Klasse",
    y = "Expressionswert"
  ) +
  theme_minimal() +
  guides(fill = "none")

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)
# ------------------------------------------------------------
# 8. Heatmap der ausgewählten LASSO-Gene
# ------------------------------------------------------------

# Probe-IDs der 10 vom LASSO ausgewählten Gene
selected_probes <- selected_genes_annotated$Gene

# Expressionswerte dieser Gene
heatmap_matrix <- exprs(ALL)[selected_probes, ]

# Gen-Symbole statt Probe-IDs als Zeilennamen
rownames(heatmap_matrix) <- selected_genes_annotated$SYMBOL

# Dimension kontrollieren
dim(heatmap_matrix)

# Klassenzugehörigkeit der Proben für die Heatmap
annotation_col <- data.frame(
  Klasse = target
)

# Probennamen müssen mit den Spalten der Expressionsmatrix übereinstimmen
rownames(annotation_col) <- colnames(heatmap_matrix)

# Kontrolle
dim(annotation_col)
table(annotation_col$Klasse)
# ------------------------------------------------------------
# 9. Heatmap der 10 LASSO-Gene
# ------------------------------------------------------------

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Expression der 10 vom LASSO ausgewählten Gene"
)
# ------------------------------------------------------------
# Zusammenfassung der Visualisierung
# ------------------------------------------------------------

# Die Visualisierungen zeigen eine klare Trennung zwischen
# B-Zell-ALL und T-Zell-ALL.
#
# Die Confusion Matrix der äußeren 5-fachen Cross-Validation
# zeigt keine Fehlklassifikationen (Accuracy = 1).
# Die ROC-Analyse bestätigt die hohe Trennleistung mit AUC = 1.
#
# Das LASSO-Modell wählte 10 Gene aus.
# CD3D, SH2D1A und TRDC besitzen positive Koeffizienten und sind
# mit der T-Zell-Klasse assoziiert, während unter anderem CD74,
# IGHM und BLNK negative Koeffizienten besitzen.
#
# Die Boxplots von CD3D und CD74 zeigen deutliche Unterschiede
# zwischen den beiden Klassen.
#
# Auch die Heatmap der 10 ausgewählten Gene zeigt unterschiedliche
# Expressionsmuster und eine klare Gruppierung der B- und T-Zell-ALL-Proben.
#
# Insgesamt stimmen die Ergebnisse der Visualisierung mit den
# Ergebnissen des LASSO-Klassifikationsmodells überein.