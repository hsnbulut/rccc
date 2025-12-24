cat("=== Running rccc reproducibility scripts ===\n\n")

# 1) Print session info
source("scripts/session_info.R")

# 2) Load all function files under R/
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
if (length(r_files) == 0) stop("No .R files found under R/")

cat("Loading functions from R/ ...\n")
for (f in r_files) {
  cat(" - ", f, "\n", sep = "")
  source(f)
}
cat("\n")

# 3) Run all real-data scripts under real_data/
example_files <- list.files("real_data", pattern = "\\.R$", full.names = TRUE)
if (length(example_files) == 0) stop("No .R files found under real_data/")

cat("Running real-data examples from real_data/ ...\n")
if (!dir.exists("results")) dir.create("results")
for (f in example_files) {
  cat(" - ", f, "\n", sep = "")
  source(f)
}
cat("\n=== DONE ===\n")
