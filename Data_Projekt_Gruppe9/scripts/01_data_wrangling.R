# ============================================================
# Data Science Projekt - Gruppe 9
# Laden des Datensatzes und erste Datenexploration.
#
# Autor: Erfan Abdollahzadehgohari
# Datum: 09.08.2026
#
# Ziel:
# Installation der benötigten Pakete,
# Laden des Datensatzes und erste Datenexploration.
# ============================================================
# ------------------------------------------------------------
# Pakete installieren
# ------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# ------------------------------------------------------------
# Benötigte Bioconductor-Pakete installieren
# ------------------------------------------------------------

if (!require("ALLMLL", quietly = TRUE))
  BiocManager::install("ALLMLL")
# ------------------------------------------------------------
# Pakete laden
# ------------------------------------------------------------

library(ALLMLL)

# ------------------------------------------------------------
# Verfügbare Datensätze anzeigen
# ------------------------------------------------------------

data(package = "ALLMLL")

# ------------------------------------------------------------
# Datensatz laden
# ------------------------------------------------------------

data("MLL.A")
# Klasse des Datensatzes anzeigen
class(MLL.A)
# ------------------------------------------------------------
# Informationen über die Proben
# ------------------------------------------------------------

sampleNames(MLL.A)
# Informationen zu den Patienten (Metadaten)

pData
# ------------------------------------------------------------
# Struktur des Datensatzes untersuchen
# ------------------------------------------------------------

MLL.A
# ------------------------------------------------------------
# Informationen über das AffyBatch-Objekt
# ------------------------------------------------------------

annotation(MLL.A)
sampleNames(MLL.A)
protocolData(MLL.A)
experimentData(MLL.A)
# Anzahl der Proben
length(sampleNames(MLL.A))

# Verwendete Plattform
annotation(MLL.A)
# ------------------------------------------------------------
# Ende des Data-Wrangling-Schritts
# Die Daten wurden erfolgreich geladen und ihre Struktur geprüft.
# Die weitere Verarbeitung erfolgt in 02_QA_QC.R.
# ------------------------------------------------------------