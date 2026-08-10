# ------------------------------------------------------------
# 1. Pakete laden
# ------------------------------------------------------------

library(ALLMLL)
library(affy)


# ------------------------------------------------------------
# 2. Rohdaten laden
# ------------------------------------------------------------

data("MLL.A")


# ------------------------------------------------------------
# 3. Quality Assessment der Rohdaten
# ------------------------------------------------------------

# Boxplot der Rohintensitäten
boxplot(
  MLL.A,
  main = "Rohdaten vor Normalisierung",
  las = 2,
  cex.axis = 0.6
)
# Dichteverteilungen der Rohintensitäten
hist(
  MLL.A,
  main = "Intensitätsverteilungen vor Normalisierung"
)
if (!require("affyPLM", quietly = TRUE))
  BiocManager::install("affyPLM")

library(affyPLM)

fit <- fitPLM(MLL.A)

RLE(fit)

NUSE(fit)

# ------------------------------------------------------------
# 4. RMA-Normalisierung
# ------------------------------------------------------------

MLL.A_rma <- rma(MLL.A)
# ------------------------------------------------------------
# 5. Kontrolle nach RMA-Normalisierung
# ------------------------------------------------------------

boxplot(
  exprs(MLL.A_rma),
  main = "Expressionswerte nach RMA-Normalisierung",
  las = 2,
  cex.axis = 0.6
)
# Nach dem RMA-Normalisierung weisen alle Proben nahezu identische Medianwerte und sehr aehnliche Verteilungen auf.
# Dies zeigt, dass techniche Unterschiede zwischen microarrays erfolfreich redutziert wurden und die expersionswert nin vergleichbar sind.
if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

library(limma)
plotDensities(
  exprs(MLL.A_rma),
  main = "Intensitätsverteilungen nach RMA-Normalisierung",
  legend = FALSE
)

#Qualitätskontrolle nach der RMA-Normalisierung: 
#Nach der RMA-Normalisierung zeigen alle 20 Microarray-Proben nahezu identische Intensitätsverteilungen. 
#Die Dichtekurven überlagern sich fast vollständig, was auf eine erfolgreiche Beseitigung systematischer technischer Unterschiede zwischen den Arrays hinweist. 
#Die Daten weisen somit eine hohe Vergleichbarkeit auf und sind für die anschließenden statistischen Analysen geeignet.

# ------------------------------------------------------------
# 6. Normalisierte Daten speichern
# ------------------------------------------------------------

saveRDS(
  MLL.A_rma,
  file = "../data/processed/MLL_A_rma.rds"
)

# ============================================================
# Zusammenfassung der Qualitätskontrolle
# ============================================================

# Die Rohdaten wurden erfolgreich auf ihre Qualität überprüft.
# Es wurden folgende Analysen durchgeführt:
# - Boxplot der Rohintensitäten
# - Dichteverteilungen der Rohdaten
# - RLE-Plot
# - NUSE-Plot
# - RMA-Normalisierung
# - Boxplot nach RMA-Normalisierung
# - Dichteplot nach RMA-Normalisierung
#
# Ergebnis:
# Die Proben zeigen nach der RMA-Normalisierung nahezu identische
# Intensitätsverteilungen. Es wurden keine auffälligen Ausreißer
# oder Qualitätsprobleme festgestellt.
#
# Die Daten sind damit für die explorative Datenanalyse
# (PCA, Clustering und Heatmaps) geeignet.
# ============================================================
