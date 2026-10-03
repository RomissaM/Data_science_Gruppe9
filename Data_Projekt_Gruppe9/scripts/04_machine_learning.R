# ============================================================
# Data Science Projekt - Gruppe 9
# Datei: 04_machine_learning.R
#
# Ziel:
# Untersuchung, ob Genexpressionsdaten zur Vorhersage einer
# Fragestellung:
# Können Genexpressionsprofile zur Unterscheidung von
# B-Zell-ALL und T-Zell-ALL verwendet werden?

# Hypothese:
# B-Zell-ALL und T-Zell-ALL weisen unterschiedliche
# Genexpressionsprofile auf, die eine Klassifikation ermöglichen.
# ============================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ALL", quietly = TRUE)) {
  BiocManager::install("ALL")
}

library(ALL)
data(ALL)

ALL
dim(exprs(ALL))

colnames(pData(ALL))

table(ALL$BT)
table(ALL$mol.biol)
# ------------------------------------------------------------
# 2. Zielvariable B-Zell-ALL vs. T-Zell-ALL erstellen
# ------------------------------------------------------------

# Vorhandene B/T-Klassen ansehen
table(ALL$BT)

# Alle B-Untergruppen zu "B" und alle T-Untergruppen zu "T"
target <- ifelse(
  grepl("^B", ALL$BT),
  "B",
  "T"
)

target <- factor(target)

# Verteilung der Zielvariable kontrollieren
table(target)
# ------------------------------------------------------------

# 3.Erstellung von Eingabematrix für Machine Learning 
# ------------------------------------------------------------

# Expressionsmatrix extrahieren und transponieren
# 
X <- t(exprs(ALL))

# Für Machine Learning (kontroll):

dim(X)
nrow(X)
length(target)
# ------------------------------------------------------------
# 4. Trainings- und Testdaten 
# ------------------------------------------------------------

set.seed(123)

# Indizes für B- und T-Proben
idx_B <- which(target == "B")
idx_T <- which(target == "T")

# Auswahl 80 % jeder Klasse für Training 
train_B <- sample(idx_B, size = round(0.8 * length(idx_B)))
train_T <- sample(idx_T, size = round(0.8 * length(idx_T)))

train_idx <- c(train_B, train_T)

# Trainingsdaten
X_train <- X[train_idx, ]
y_train <- target[train_idx]

# Testdaten
X_test <- X[-train_idx, ]
y_test <- target[-train_idx]
dim(X_train)
dim(X_test)

table(y_train)
table(y_test)
# ------------------------------------------------------------
# 5. Feature-Auswahl
# ------------------------------------------------------------

# Berechnung der Varianz jedes Probe-Sets anhand der Trainingsdaten
gene_var <- apply(X_train, 2, var)

# Auswahl der 100 Probe-Sets mit der höchsten Varianz 
top_genes <- order(gene_var, decreasing = TRUE)[1:100]

# Reduzierung von Trainings- und Testdaten auf dieselben Probe-Sets 
X_train_selected <- X_train[, top_genes]
X_test_selected  <- X_test[, top_genes]

# Dimensionen kontrolle
dim(X_train_selected)
dim(X_test_selected)
# ------------------------------------------------------------
# 6. LASSO-logistische Regression
# ------------------------------------------------------------

# Paket installieren

if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet", repos = "https://cloud.r-project.org")
}

library(glmnet)

# Zielvariable numerisch codieren:
# B = 0, T = 1

y_train_num <- ifelse(y_train == "T", 1, 0)

table(y_train_num)

# ------------------------------------------------------------
# 7. LASSO-Modell mit Cross-Validation
# ------------------------------------------------------------

set.seed(123)

cv_model <- cv.glmnet(
  x = X_train_selected,
  y = y_train_num,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)

plot(cv_model)

# Optimale Lambda-Werte 

cv_model$lambda.min
cv_model$lambda.1se
coef_min <- coef(cv_model, s = "lambda.min")

# Anzahl der Koeffizienten ungleich 0
sum(coef_min != 0)

coef_min

# ------------------------------------------------------------
# 8. Vorhersage auf den Testdaten
# ------------------------------------------------------------

# Wahrscheinlichkeit für T-Zell-Klasse vorhersagen

prob_test <- predict(
  cv_model,
  newx = X_test_selected,
  s = "lambda.min",
  type = "response"
)

# Umwandlung In B bzw. T 

pred_test <- ifelse(prob_test >= 0.5, "T", "B")

length(pred_test)

# Vergleich des Vorhersage mit tatsächlicher Klasse 

table(
  Tatsächlich = y_test,
  Vorhergesagt = pred_test
)

# ------------------------------------------------------------
# 9. Modellgüte auf den Testdaten
# ------------------------------------------------------------

accuracy <- mean(pred_test == y_test)
accuracy

# ------------------------------------------------------------
# 10. Sensitivität und Spezifität
# ------------------------------------------------------------

# Confusion Matrix 

cm <- table(
  Tatsächlich = y_test,
  Vorhergesagt = pred_test
)

cm

# Werte aus der Confusion Matrix

TP <- cm["T", "T"]
TN <- cm["B", "B"]
FP <- cm["B", "T"]
FN <- cm["T", "B"]

# Sensitivität: Anteil korrekt erkannter T-Proben
sensitivity <- TP / (TP + FN)

# Spezifität: Anteil korrekt erkannter B-Proben
specificity <- TN / (TN + FP)

sensitivity
specificity
# ------------------------------------------------------------
# 11. Robustheit der Klassifikation
# ------------------------------------------------------------

# Für reproduzierbare Ergebnisse
set.seed(123)

# Anzahl der Folds
k <- 5

# Speicher für Vorhersagen aller 128 Proben
cv_predictions <- rep(NA, length(target))

length(cv_predictions)
# B- und T-Proben getrennt bestimmen
idx_B <- which(target == "B")
idx_T <- which(target == "T")

# Zufällig auf 5 Folds verteilen
fold_B <- sample(rep(1:k, length.out = length(idx_B)))
fold_T <- sample(rep(1:k, length.out = length(idx_T)))

# Fold-Zugehörigkeit speichern
fold_id <- rep(NA, length(target))

fold_id[idx_B] <- fold_B
fold_id[idx_T] <- fold_T

# Verteilung kontrollieren
table(Fold = fold_id, Klasse = target)
# ------------------------------------------------------------
# 11.3 Äußere 5-fache Cross-Validation
# ------------------------------------------------------------

cv_probabilities <- rep(NA, length(target))

# Alte Vorhersagen:
cv_predictions <- rep(NA, length(target))

for (i in 1:k) {
  
  # Trainings- und Testproben des aktuellen Folds
  train_i <- which(fold_id != i)
  test_i  <- which(fold_id == i)
  
  # Trainings- und Testdaten
  X_train_i <- X[train_i, ]
  X_test_i  <- X[test_i, ]
  
  y_train_i <- target[train_i]
  
  # Feature-Auswahl NUR anhand der Trainingsdaten
  gene_var_i <- apply(X_train_i, 2, var)
  top_genes_i <- order(gene_var_i, decreasing = TRUE)[1:100]
  
  X_train_i <- X_train_i[, top_genes_i]
  X_test_i  <- X_test_i[, top_genes_i]
  
  # B = 0, T = 1
  y_train_num_i <- ifelse(y_train_i == "T", 1, 0)
  
  # LASSO-Modell mit interner Cross-Validation
  model_i <- cv.glmnet(
    x = X_train_i,
    y = y_train_num_i,
    family = "binomial",
    alpha = 1,
    nfolds = 5
  )
  
  # Vorhersage für den jeweils ausgelassenen Fold
  prob_i <- predict(
    model_i,
    newx = X_test_i,
    s = "lambda.min",
    type = "response"
  )
cv_probabilities[test_i] <- as.numeric(prob_i)
  
  cv_predictions[test_i] <- ifelse(prob_i >= 0.5, "T", "B")
}
length(cv_probabilities)
sum(is.na(cv_probabilities))
# ------------------------------------------------------------
# 11.4 Ergebnisse der äußeren Cross-Validation
# ------------------------------------------------------------

# Kontrolle: liegt für alle Proben eine Vorhersage vor?
sum(is.na(cv_predictions))

# Confusion Matrix
cv_cm <- table(
  Tatsächlich = target,
  Vorhergesagt = cv_predictions
)

cv_cm

# Accuracy
cv_accuracy <- mean(cv_predictions == target)

cv_accuracy
dim(X)

dim(X_train)
dim(X_test)

dim(X_train_selected)
dim(X_test_selected)

length(y_train)
length(y_test)

length(pred_test)
dim(X_train)
dim(X_test)

dim(X_train_selected)
dim(X_test_selected)

# ------------------------------------------------------------
# Interpretation der äußeren Cross-Validation
# ------------------------------------------------------------

# Die äußere 5-fache Cross-Validation ergab eine Accuracy von 100 %.
# Alle 95 B-Zell- und alle 33 T-Zell-Proben wurden korrekt klassifiziert.
# Die Feature-Auswahl wurde in jedem Fold ausschließlich anhand der
# jeweiligen Trainingsdaten durchgeführt, wodurch Data Leakage vermieden wurde.
#
# Das Ergebnis deutet darauf hin, dass sich B- und T-Zell-ALL anhand
# der Genexpressionsprofile in diesem Datensatz sehr deutlich unterscheiden.
# Aufgrund der begrenzten Stichprobengröße sollte die Generalisierbarkeit
# auf unabhängige externe Datensätze dennoch vorsichtig interpretiert werden.


# ------------------------------------------------------------
# 12. ROC-Kurve und AUC
# ------------------------------------------------------------

# Paket install

if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC", repos = "https://cloud.r-project.org")
}

library(pROC)

# Numerische Zielvariable für die Testdaten:
# B = 0, T = 1
y_test_num <- ifelse(y_test == "T", 1, 0)

# ROC-Kurve aus den vorhergesagten Wahrscheinlichkeiten
roc_test <- roc(
  response = y_test_num,
  predictor = as.numeric(prob_test)
)

# AUC Berechnung
auc_test <- auc(roc_test)

auc_test

# ------------------------------------------------------------
# 13. Vom LASSO ausgewählte Probe-Sets
# ------------------------------------------------------------

# Koeffizienten des finalen Modells bei lambda.min

coef_lasso <- coef(cv_model, s = "lambda.min")

# Koeffizienten als Matrix umwandeln
coef_matrix <- as.matrix(coef_lasso)

# Probe-Sets mit Koeffizient ungleich 0 bestimmen
selected_idx <- which(coef_matrix[, 1] != 0)

# Intercept entfernen
selected_idx <- selected_idx[
  rownames(coef_matrix)[selected_idx] != "(Intercept)"
]

# Ausgewählte Probe-Sets und ihre Koeffizienten anzeigen
selected_genes <- data.frame(
  ProbeSet = rownames(coef_matrix)[selected_idx],
  Koeffizient = coef_matrix[selected_idx, 1]
)

selected_genes
nrow(selected_genes)
# ------------------------------------------------------------
# 13.1 Annotation der ausgewählten Gene
# ------------------------------------------------------------

# Passendes Annotationspaket für den hgu95av2-Chip
if (!requireNamespace("hgu95av2.db", quietly = TRUE)) {
  BiocManager::install("hgu95av2.db")
}

library(hgu95av2.db)
library(AnnotationDbi)

# Probe-Set-IDs annotieren
gene_annotation <- AnnotationDbi::select(
  hgu95av2.db,
  keys = selected_genes$ProbeSet,
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME")
)

gene_annotation
# ------------------------------------------------------------
# 13.2 Annotation mit LASSO-Koeffizienten verbinden
# ------------------------------------------------------------

selected_genes_annotated <- merge(
  selected_genes,
  gene_annotation,
  by.x = "ProbeSet",
  by.y = "PROBEID"
)

selected_genes_annotated

write.csv(
  selected_genes_annotated,
  "../data/processed/lasso_selected_genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Zusammenfassung der Machine-Learning-Analyse
# ------------------------------------------------------------

# Ziel:
# Untersuchung, ob B- und T-Zell-ALL anhand ihrer
# Genexpressionsprofile klassifiziert werden können.
#
# Datensatz:
# 128 Proben
# B-Zell-ALL: 95
# T-Zell-ALL: 33
#
# Methode:
#- Auswahl der 100 variabelsten Probe-Sets anhand der Trainingsdaten
# - LASSO-logistische Regression
# - interne Cross-Validation zur Wahl von lambda
# - äußere 5-fache Cross-Validation zur Modellbewertung
#
# Ergebnisse:
# Train/Test-Accuracy = 1.00
# Sensitivität = 1.00
# Spezifität = 1.00
# AUC = 1.00
#
# Äußere 5-fache Cross-Validation:
# 95/95 B-Zell-Proben korrekt klassifiziert
# 33/33 T-Zell-Proben korrekt klassifiziert
# Accuracy = 1.00
#
# Das finale LASSO-Modell verwendete 10 Probe-Sets.
# Unter den ausgewählten Genen befinden sich unter anderem
# CD3D, TRDC, SH2D1A, BLNK, CD74 und IGHM.
#
# Die ausgewählten Gene zeigen eine biologisch plausible
# Beziehung zur Unterscheidung von B- und T-Zell-Proben.
#
# Trotz der sehr hohen Klassifikationsleistung sollte die
# Übertragbarkeit auf unabhängige externe Datensätze
# vorsichtig interpretiert werden.