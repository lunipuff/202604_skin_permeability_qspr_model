library(dplyr)
source("01_logL.R")
source("00loading.R")

raw <- read.csv("00cleaned_table.csv", row.names=1)
cv <- read.csv("03cvbestmodel.csv", row.names=1)
df <- raw
df$Compound <- trimws(df$Compound)

############################################################
# Predict logkpl for candidate in vitro compounds using lm()
############################################################

## 1. Current model formula

model <- "logkpl ~ MWa + logKowb + I(logKowb^2) +
          Mptc + I(Mptc^2) +
          LogSaqd + I(LogSaqd^2) +
          Texpi + I(Texpi^2)"


## 2. Candidate compound predictor table

new_compounds <- data.frame(
  Compound = c(
    "Lucifer Yellow CH dipotassium",
    "FITC-dextran 4 kDa",
    "Caffeine",
    "5-Fluorouracil",
    "Lidocaine",
    "Estradiol",
    "Methotrexate",
    "Melphalan"
  ),
  MWa = c(
    521.6,
    4000,
    194.19,
    130.08,
    234.34,
    272.38,
    454.44,
    305.20
  ),
  logKowb = c(
    NA,
    NA,
    -0.07,
    -0.58,
    2.30,
    3.57,
    -2.20,
    0.40
  ),
  Mptc = c(
    NA,
    NA,
    508.15,
    557.15,
    340.65,
    446.15,
    468.15,
    455.65
  ),
  LogSaqd = c(
    NA,
    NA,
    -3.952,
    -4.028,
    -5.597,
    -7.879,
    -6.427,
    -5.928
  ),
  Texpi = rep(298.15, 8)
)

new_compounds <- new_compounds[which(!new_compounds$Compound %in% c("Methotrexate",
    "Melphalan")),]
    
## 3. Fit lm() model on training data

fit <- lm(as.formula(model), data = df)

summary(fit)


## 4. Keep only compounds with complete predictor values

predictor_cols <- c("MWa", "logKowb", "Mptc", "LogSaqd", "Texpi")

predict_compounds <- new_compounds[
  complete.cases(new_compounds[, predictor_cols]),
]

excluded_compounds <- new_compounds[
  !complete.cases(new_compounds[, predictor_cols]),
]

cat("Compounds excluded from prediction due to missing predictors:\n")
print(excluded_compounds[, c("Compound", predictor_cols)])


## 5. Predict logkpl

pred <- predict(
  fit,
  newdata = predict_compounds,
  interval = "prediction",
  level = 0.95
)

pred <- as.data.frame(pred)

names(pred) <- c(
  "predicted_logkpl",
  "prediction_lower95",
  "prediction_upper95"
)


## 6. Combine predictions with compound information

prediction_results <- cbind(
  predict_compounds,
  pred
)


## 7. Add model-domain flags based on training-data ranges

domain_summary <- data.frame()

for (p in predictor_cols) {
  x <- df[[p]]

  domain_summary <- rbind(
    domain_summary,
    data.frame(
      predictor = p,
      lower_2.5 = quantile(x, 0.025, na.rm = TRUE),
      upper_97.5 = quantile(x, 0.975, na.rm = TRUE),
      minimum = min(x, na.rm = TRUE),
      maximum = max(x, na.rm = TRUE)
    )
  )
}

prediction_results$inside_2.5_97.5_domain <- TRUE
prediction_results$outside_predictors <- ""

for (i in seq_len(nrow(prediction_results))) {
  outside <- c()

  for (p in predictor_cols) {
    lower <- domain_summary$lower_2.5[domain_summary$predictor == p]
    upper <- domain_summary$upper_97.5[domain_summary$predictor == p]
    value <- prediction_results[i, p]

    if (value < lower | value > upper) {
      outside <- c(outside, p)
    }
  }

  if (length(outside) > 0) {
    prediction_results$inside_2.5_97.5_domain[i] <- FALSE
    prediction_results$outside_predictors[i] <- paste(outside, collapse = ", ")
  }
}


## 8. Print final result

prediction_results <- prediction_results[order(prediction_results$predicted_logkpl), ]

print(prediction_results)


## 9. Save outputs

write.csv(
  prediction_results,
  "candidate_compound_logkpl_predictions.csv",
  row.names = FALSE
)

write.csv(
  excluded_compounds,
  "candidate_compounds_excluded_from_prediction.csv",
  row.names = FALSE
)

write.csv(
  domain_summary,
  "training_predictor_domain_summary.csv",
  row.names = FALSE
)

############################################################
# Plot predicted logkpl as vertical compound markers
############################################################

## Order compounds by predicted logkpl
plot_df <- prediction_results[order(prediction_results$predicted_logkpl), ]

## Basic plot limits
x_min <- min(plot_df$prediction_lower95, plot_df$predicted_logkpl, na.rm = TRUE)
x_max <- max(plot_df$prediction_upper95, plot_df$predicted_logkpl, na.rm = TRUE)

## Add a little spacing
x_pad <- 0.1 * (x_max - x_min)
xlim_use <- c(x_min - x_pad, x_max + x_pad)

## Open PDF
pdf("candidate_compound_predicted_logkpl_axis.pdf",
    width = 11,
    height = 4)

par(mar = c(5, 2, 3, 2))

## Empty plot
plot(
  NA,
  xlim = xlim_use,
  ylim = c(0, 1),
  xlab = "Predicted logKp",
  ylab = "",
  yaxt = "n",
  main = "Predicted logKp of candidate in vitro compounds"
)

## Horizontal reference axis
abline(h = 0.5, lwd = 2)

## Add vertical lines and labels
for (i in seq_len(nrow(plot_df))) {
  x <- plot_df$predicted_logkpl[i]
  label <- plot_df$Compound[i]

  abline(v = x, lty = 2, lwd = 1)

  points(x, 0.5, pch = 16)

  text(
    x = x,
    y = 0.62 + 0.12 * (i %% 2),
    labels = label,
    srt = 45,
    adj = 0,
    cex = 0.8
  )
}

dev.off()