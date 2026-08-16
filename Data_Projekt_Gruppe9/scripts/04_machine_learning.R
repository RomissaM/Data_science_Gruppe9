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

BiocManager::install("ALL")
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
# 3. Eingabematrix für Machine Learning erstellen
# ------------------------------------------------------------

# Expressionsmatrix extrahieren
X <- exprs(ALL)


# Für Machine Learning:
# Samples müssen in den Zeilen und Gene in den Spalten stehen
X <- t(X)

dim(X)
# ------------------------------------------------------------
# 3. Eingabematrix für Machine Learning erstellen
# ------------------------------------------------------------

# Expressionsmatrix extrahieren und transponieren
# Patienten = Zeilen, Gene = Spalten
X <- t(exprs(ALL))

# Für Machine Learning:
# Samples müssen in den Zeilen und Gene in den Spalten stehen

dim(X)
nrow(X)
length(target)
# ------------------------------------------------------------
# 4. Trainings- und Testdaten erstellen
# ------------------------------------------------------------

set.seed(123)

# Indizes für B- und T-Proben
idx_B <- which(target == "B")
idx_T <- which(target == "T")

# 80 % jeder Klasse für Training auswählen
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

# Varianz jedes Gens nur anhand der Trainingsdaten berechnen
gene_var <- apply(X_train, 2, var)

# 100 Gene mit der höchsten Varianz auswählen
top_genes <- order(gene_var, decreasing = TRUE)[1:100]

# Trainings- und Testdaten auf dieselben Gene reduzieren
X_train_selected <- X_train[, top_genes]
X_test_selected  <- X_test[, top_genes]

# Dimensionen kontrollieren
dim(X_train_selected)
dim(X_test_selected)
# ------------------------------------------------------------
# 6. LASSO-logistische Regression
# ------------------------------------------------------------

# Paket installieren, falls noch nicht vorhanden
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
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
# Optimale Lambda-Werte anzeigen
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

# In B bzw. T umwandeln
pred_test <- ifelse(prob_test >= 0.5, "T", "B")

length(pred_test)
# Vorhersage mit tatsächlicher Klasse vergleichen
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

# Confusion Matrix speichern
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

# Alte Vorhersagen zurücksetzen
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
  
  cv_predictions[test_i] <- ifelse(prob_i >= 0.5, "T", "B")
}
# ------------------------------------------------------------
# 11.4 Ergebnisse der äußeren Cross-Validation
# ------------------------------------------------------------

# Prüfen, ob für alle Proben eine Vorhersage vorliegt
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
