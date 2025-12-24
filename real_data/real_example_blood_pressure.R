############################################################
## Real Data Application 2: Sistolic Blood Pressure
############################################################

cat("\n--- Running: real_example_blood_pressure.R ---\n")


library(cccrm)
library(dplyr)
library(tidyr)
library(rrcov)
library(MASS)
library(L1pack)

#========================================================
# 1) Prepare wide data: manual vs automatic, using SIS
#========================================================
bp_wide <- bpres %>%
  filter(METODE %in% c(1, 2)) %>%
  mutate(METHOD = factor(METODE, levels = c(1, 2),
                         labels = c("manual", "automatic"))) %>%
  group_by(ID, METHOD) %>%
  summarise(SIS = median(SIS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = METHOD, values_from = SIS) %>%
  filter(!is.na(manual), !is.na(automatic))

x <- bp_wide$manual
y <- bp_wide$automatic

#========================================================
# 2) Clean data: compute coefficients
#========================================================
res_clean <- c(
  Pearson      = cor(x, y),
  CCC          = ccc(x, y),
  rCCC         = rccc(x, y),
  `L1-CCC`     = ccc_L1(x, y),
  `Bayes-CCC`  = ccc_bayes_mvt(x, y),
  Absolute     = robust_ccc(x, y, method = "absolute"),
  `Absolute 1.5` = robust_ccc(x, y, method = "absolute1.5"),
  Winsor       = robust_ccc(x, y, method = "winsor"),
  Huber        = robust_ccc(x, y, method = "huber")
)

#========================================================
# 3) Contaminated data: add one artificial outlier
#    (You can tweak these values if you want a stronger disruption)
#========================================================
set.seed(123)

k_prop <- 0.10              # contamination proportion (e.g., 0.05, 0.10, 0.20)
m <- ceiling(k_prop * length(x))

# robust-ish scale for extreme construction
sx <- sd(x, na.rm = TRUE)
sy <- sd(y, na.rm = TRUE)

# choose "strength" of outliers (increase to make it harsher)
K <- 10   # e.g., 10, 15, 20

# Create anti-concordant outliers:
# x very high, y very low (breaks both correlation and concordance)
x_bad <- rep(median(x, na.rm = TRUE) + K * sx, m)
y_bad <- rep(median(y, na.rm = TRUE) - K * sy, m)

# Optionally add slight jitter so they are not identical points
x_bad <- x_bad + rnorm(m, 0, 0.1 * sx)
y_bad <- y_bad + rnorm(m, 0, 0.1 * sy)

x_out <- c(x, x_bad)
y_out <- c(y, y_bad)


res_cont <- c(
  Pearson      = cor(x_out, y_out),
  CCC          = ccc(x_out, y_out),
  rCCC         = rccc(x_out, y_out),
  `L1-CCC`     = ccc_L1(x_out, y_out),
  `Bayes-CCC`  = ccc_bayes_mvt(x_out, y_out),
  Absolute     = robust_ccc(x_out, y_out, method = "absolute"),
  `Absolute 1.5` = robust_ccc(x_out, y_out, method = "absolute1.5"),
  Winsor       = robust_ccc(x_out, y_out, method = "winsor"),
  Huber        = robust_ccc(x_out, y_out, method = "huber")
)

#========================================================
# 4) Final table
#========================================================
result_bpres_SIS <- cbind(Clean = res_clean, Contaminated = res_cont)

write.csv(
  round(result_bpres_SIS, 3),
  file = file.path("results", "table_blood_pressure_clean_vs_contaminated.csv")
)

cat("\n--- Blood Pressure results ---\n")
print(round(result_bpres_SIS, 3))
