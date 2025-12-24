############################################################
## Real Data Application 1: Glucose
############################################################

cat("\n--- Running: real_example_glucose.R ---\n")

library(MethComp)
library(dplyr)
library(tidyr)
library(rrcov)
data(glucose)

# -----------------------------
# Starting, plasma vs serum
# -----------------------------
g0 <- glucose %>%
  filter(time == 0, type %in% c("plasma","serum")) %>%
  transmute(id = item, meth = type, y)

wide0 <- g0 %>%
  group_by(id, meth) %>%
  summarise(y = median(y, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = id, names_from = meth, values_from = y) %>%
  filter(!is.na(plasma), !is.na(serum))

x <- wide0$plasma
y <- wide0$serum

# -----------------------------
# Clean data
# -----------------------------
pearson.clean <- cor(x, y)
ccc.clean     <- ccc(x, y)
rccc.clean    <- rccc(x, y)
l1.clean      <- ccc_L1(x, y)
bayes.clean   <- ccc_bayes_mvt(x, y)

abs.clean     <- robust_ccc(x, y, method = "absolute")
abs1.5.clean  <- robust_ccc(x, y, method = "absolute1.5")
winsor.clean     <- robust_ccc(x, y, method = "winsor")
huber.clean      <- robust_ccc(x, y, method = "huber")

# -----------------------------
# Contaminated data
# -----------------------------
x_out <- c(x, 15)
y_out <- c(y,  2)

pearson.contaminated <- cor(x_out, y_out)
ccc.contaminated     <- ccc(x_out, y_out)
rccc.contaminated    <- rccc(x_out, y_out)
l1.contaminated      <- ccc_L1(x_out, y_out)
bayes.contaminated   <- ccc_bayes_mvt(x_out, y_out)

abs.contaminated     <- robust_ccc(x_out, y_out, method = "absolute")
abs1.5.contaminated  <- robust_ccc(x_out, y_out, method = "absolute1.5")
winsor.contaminated     <- robust_ccc(x_out, y_out, method = "winsor")
huber.contaminated      <- robust_ccc(x_out, y_out, method = "huber")

# -----------------------------
# Result Table
# -----------------------------
clean.data <- c(
  pearson.clean, ccc.clean, rccc.clean,
  l1.clean, bayes.clean,
  abs.clean, abs1.5.clean,
  winsor.clean, huber.clean
)

contaminated.data <- c(
  pearson.contaminated, ccc.contaminated, rccc.contaminated,
  l1.contaminated, bayes.contaminated,
  abs.contaminated, abs1.5.contaminated,
  winsor.contaminated, huber.contaminated
)

result.glucose <- cbind(Clean = clean.data, Contaminated= contaminated.data)

rownames(result.glucose) <- c(
  "Pearson", "CCC", "rCCC",
  "L1-CCC", "Bayes-CCC",
  "Absolute", "Absolute 1.5",
  "Winsor", "Huber"
)

if (!dir.exists("results")) dir.create("results")

write.csv(
  round(result.glucose, 3),
  file = file.path("results", "table_glucose_clean_vs_contaminated.csv")
)


cat("\n--- Glucose results ---\n")
print(round(result.glucose,3))

