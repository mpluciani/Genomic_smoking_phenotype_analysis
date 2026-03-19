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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snake_case <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
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

write.csv(
  data.frame(dataset = c("genomic_associations", "smoking_clinical_phenotypes"), rows = c(nrow(associations), nrow(smoking))),
  file.path(output_dir, "genomic_and_smoking_dataset_sizes.csv"),
  row.names = FALSE
)

association_summary <- aggregate(
  association_score ~ location + cpg + gwas + eqtl,
  data = associations,
  FUN = mean
)
names(association_summary)[names(association_summary) == "association_score"] <- "mean_association_score"
write.csv(association_summary, file.path(output_dir, "genomic_association_score_summary.csv"), row.names = FALSE)

smoking_summary <- aggregate(
  cbind(eyesight_left, eyesight_right, height_cm, weight_kg) ~ smoking + gender,
  data = smoking,
  FUN = mean
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

message("Genomic association and smoking phenotype analysis complete. Outputs written to: ", output_dir)
