project_root <- normalizePath(file.path(getwd()))
if (!file.exists(file.path(project_root, "data", "raw"))) {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_dir <- dirname(normalizePath(script_path))
    project_root <- normalizePath(file.path(script_dir, ".."))
  }
}

data_dir <- file.path(project_root, "data", "raw")
output_dir <- file.path(project_root, "outputs")
figure_dir <- file.path(output_dir, "figure_outcomes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

snake_case <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

safe_mean <- function(x) {
  mean(x, na.rm = TRUE)
}

pretty_term <- function(term) {
  pretty_map <- c(
    "smokingsmoker" = "smoker (vs non_smoker)",
    "genderM" = "male (vs female)",
    "smokingsmoker:genderM" = "smoker x male",
    "height_cm" = "height (cm)",
    "weight_kg" = "weight (kg)",
    "waist_cm" = "waist (cm)",
    "cholesterol" = "cholesterol",
    "triglyceride" = "triglyceride",
    "hemoglobin" = "hemoglobin",
    "age" = "age"
  )
  ifelse(term %in% names(pretty_map), pretty_map[term], term)
}

associations <- read.csv(file.path(data_dir, "genomic_associations.csv"), stringsAsFactors = FALSE)
smoking <- read.csv(file.path(data_dir, "smoking.csv"), stringsAsFactors = FALSE)

names(associations) <- snake_case(names(associations))
names(smoking) <- snake_case(names(smoking))

if ("smoking" %in% names(smoking)) {
  smoking$smoking <- factor(smoking$smoking, levels = c(0, 1), labels = c("non_smoker", "smoker"))
}
if ("gender" %in% names(smoking)) {
  smoking$gender <- factor(smoking$gender)
}
if ("oral" %in% names(smoking)) {
  smoking$oral <- factor(smoking$oral)
}
if ("tartar" %in% names(smoking)) {
  smoking$tartar <- factor(smoking$tartar)
}

dataset_sizes <- data.frame(
  dataset = c("genomic_associations", "smoking_clinical_phenotypes"),
  rows = c(nrow(associations), nrow(smoking))
)
write.csv(
  dataset_sizes,
  file.path(output_dir, "genomic_and_smoking_dataset_sizes.csv"),
  row.names = FALSE
)

association_summary <- aggregate(
  association_score ~ location + cpg + gwas + eqtl,
  data = associations,
  FUN = safe_mean
)
names(association_summary)[names(association_summary) == "association_score"] <- "mean_association_score"
write.csv(association_summary, file.path(output_dir, "genomic_association_score_summary.csv"), row.names = FALSE)

smoking_summary <- aggregate(
  cbind(eyesight_left, eyesight_right, height_cm, weight_kg) ~ smoking + gender,
  data = smoking,
  FUN = safe_mean
)
write.csv(smoking_summary, file.path(output_dir, "smoking_gender_clinical_means.csv"), row.names = FALSE)

left_eye_model <- lm(eyesight_left ~ smoking * gender + height_cm + weight_kg, data = smoking)
capture.output(summary(left_eye_model), file = file.path(output_dir, "linear_model_left_eyesight_smoking_gender.txt"))

smoking_glm <- glm(
  smoking ~ age + height_cm + weight_kg + waist_cm + cholesterol + triglyceride + hemoglobin + gender,
  data = smoking,
  family = binomial()
)
capture.output(summary(smoking_glm), file = file.path(output_dir, "logistic_model_smoking_status_clinical_predictors.txt"))

# Figure 1: dataset sizes
png(file.path(figure_dir, "01_dataset_sizes_barplot.png"), width = 1400, height = 900, res = 150)
par(mar = c(8, 5, 4, 2) + 0.1)
dataset_labels <- gsub("_", " ", dataset_sizes$dataset)
dataset_colors <- c("#4E79A7", "#F28E2B")
dataset_bar_positions <- barplot(
  height = dataset_sizes$rows,
  names.arg = dataset_labels,
  las = 2,
  col = dataset_colors,
  ylab = "Rows",
  main = "Dataset Sizes"
)
text(
  x = dataset_bar_positions,
  y = dataset_sizes$rows,
  labels = format(dataset_sizes$rows, big.mark = ","),
  pos = 3,
  cex = 0.9
)
dev.off()

# Figure 2: association score summary by annotation combination
association_plot <- association_summary[order(association_summary$mean_association_score, decreasing = TRUE), ]
association_labels <- paste0(
  "L", association_plot$location,
  " C", association_plot$cpg,
  " G", association_plot$gwas,
  " E", association_plot$eqtl
)
png(file.path(figure_dir, "02_genomic_association_summary_barplot.png"), width = 1600, height = 900, res = 150)
par(mar = c(9, 5, 4, 2) + 0.1)
association_bar_positions <- barplot(
  height = association_plot$mean_association_score,
  names.arg = association_labels,
  las = 2,
  col = "#59A14F",
  ylab = "Mean Association Score",
  main = "Mean Genomic Association Score by Annotation Pattern"
)
text(
  x = association_bar_positions,
  y = association_plot$mean_association_score,
  labels = sprintf("%.2f", association_plot$mean_association_score),
  pos = 3,
  cex = 0.85
)
mtext("Pattern key: L=location, C=CpG, G=GWAS, E=eQTL", side = 1, line = 7, cex = 0.85)
dev.off()

# Figure 3: smoking x gender phenotype means
smoking_plot <- smoking_summary[order(smoking_summary$smoking, smoking_summary$gender), ]
group_labels <- paste(smoking_plot$smoking, smoking_plot$gender, sep = " | ")
group_colors <- ifelse(smoking_plot$smoking == "smoker", "#E15759", "#76B7B2")

draw_panel <- function(values, panel_title, y_label) {
  panel_positions <- barplot(
    height = values,
    names.arg = group_labels,
    las = 2,
    col = group_colors,
    ylab = y_label,
    main = panel_title
  )
  text(
    x = panel_positions,
    y = values,
    labels = sprintf("%.2f", values),
    pos = 3,
    cex = 0.8
  )
}

png(file.path(figure_dir, "03_smoking_gender_clinical_means_panels.png"), width = 1700, height = 1200, res = 150)
par(mfrow = c(2, 2), mar = c(8, 5, 4, 2) + 0.1)
draw_panel(smoking_plot$eyesight_left, "Eyesight Left Mean", "Mean")
legend("topleft", legend = c("non_smoker", "smoker"), fill = c("#76B7B2", "#E15759"), bty = "n", cex = 0.9)
draw_panel(smoking_plot$eyesight_right, "Eyesight Right Mean", "Mean")
draw_panel(smoking_plot$height_cm, "Height Mean", "Centimeters")
draw_panel(smoking_plot$weight_kg, "Weight Mean", "Kilograms")
dev.off()

# Figure 4: linear model coefficient estimates with 95% CI
left_eye_coef <- summary(left_eye_model)$coefficients
left_eye_coef <- data.frame(
  term = rownames(left_eye_coef),
  estimate = left_eye_coef[, "Estimate"],
  std_error = left_eye_coef[, "Std. Error"],
  stringsAsFactors = FALSE
)
left_eye_coef <- left_eye_coef[left_eye_coef$term != "(Intercept)", ]
left_eye_coef$ci_low <- left_eye_coef$estimate - 1.96 * left_eye_coef$std_error
left_eye_coef$ci_high <- left_eye_coef$estimate + 1.96 * left_eye_coef$std_error
left_eye_coef$label <- pretty_term(left_eye_coef$term)
left_eye_coef <- left_eye_coef[order(left_eye_coef$estimate), ]

png(file.path(figure_dir, "04_linear_model_coefficients.png"), width = 1500, height = 900, res = 150)
par(mar = c(6, 14, 4, 2) + 0.1)
y_linear <- seq_len(nrow(left_eye_coef))
plot(
  left_eye_coef$estimate,
  y_linear,
  pch = 19,
  col = "#4E79A7",
  xlab = "Coefficient Estimate (95% CI)",
  ylab = "",
  yaxt = "n",
  main = "Linear Model: Left Eyesight Predictors",
  xlim = range(c(left_eye_coef$ci_low, left_eye_coef$ci_high))
)
segments(
  x0 = left_eye_coef$ci_low,
  y0 = y_linear,
  x1 = left_eye_coef$ci_high,
  y1 = y_linear,
  lwd = 2,
  col = "#4E79A7"
)
abline(v = 0, lty = 2, col = "gray40")
axis(2, at = y_linear, labels = left_eye_coef$label, las = 2)
dev.off()

# Figure 5: logistic model odds ratios with 95% CI
smoking_coef <- summary(smoking_glm)$coefficients
smoking_coef <- data.frame(
  term = rownames(smoking_coef),
  estimate = smoking_coef[, "Estimate"],
  std_error = smoking_coef[, "Std. Error"],
  stringsAsFactors = FALSE
)
smoking_coef <- smoking_coef[smoking_coef$term != "(Intercept)", ]
smoking_coef$odds_ratio <- exp(smoking_coef$estimate)
smoking_coef$ci_low <- exp(smoking_coef$estimate - 1.96 * smoking_coef$std_error)
smoking_coef$ci_high <- exp(smoking_coef$estimate + 1.96 * smoking_coef$std_error)
smoking_coef$label <- pretty_term(smoking_coef$term)
smoking_coef <- smoking_coef[order(smoking_coef$odds_ratio), ]

png(file.path(figure_dir, "05_logistic_model_odds_ratios.png"), width = 1500, height = 1000, res = 150)
par(mar = c(6, 14, 4, 2) + 0.1)
y_logistic <- seq_len(nrow(smoking_coef))
plot(
  smoking_coef$odds_ratio,
  y_logistic,
  pch = 19,
  col = "#F28E2B",
  log = "x",
  xlab = "Odds Ratio (log scale, 95% CI)",
  ylab = "",
  yaxt = "n",
  main = "Logistic Model: Smoking Status Predictors",
  xlim = range(c(smoking_coef$ci_low, smoking_coef$ci_high))
)
segments(
  x0 = smoking_coef$ci_low,
  y0 = y_logistic,
  x1 = smoking_coef$ci_high,
  y1 = y_logistic,
  lwd = 2,
  col = "#F28E2B"
)
abline(v = 1, lty = 2, col = "gray40")
axis(2, at = y_logistic, labels = smoking_coef$label, las = 2)
dev.off()

message("Genomic association and smoking phenotype analysis complete. Outputs written to: ", output_dir)
message("Outcome figures written to: ", figure_dir)
