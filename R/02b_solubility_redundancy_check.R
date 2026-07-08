############################################################
# 02b_solubility_redundancy_check.R
# Check redundancy among logKow, melting point, and solubility descriptors
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset
############################################################

df <- read.csv(path_cleaned_dataset)

required_cols <- c(
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce"
)

missing_cols <- required_cols[!(required_cols %in% names(df))]

if (length(missing_cols) > 0) {
	stop(
		paste(
			"Missing required columns:",
		paste(missing_cols, collapse = ", ")
		)
	)
}

############################################################
# Convert relevant columns to numeric
############################################################

for (col in required_cols) {
	df[[col]] <- as.numeric(df[[col]])
}

############################################################
# General Solubility Equation estimate
############################################################

# The General Solubility Equation commonly uses melting point in Celsius:
# logS = 0.5 - 0.01 * (MP - 25) - logKow
#
# Mptc in this dataset is recorded in Kelvin, so convert to Celsius.

df$Mpt_C <- df$Mptc - 273.15

df$GSE_logS <- 0.5 - 0.01 * (df$Mpt_C - 25) - df$logKowb

############################################################
# Correlation matrix
############################################################

cor_vars <- c(
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"GSE_logS"
)

cor_mat <- cor(
	df[, cor_vars],
	use = "complete.obs"
)

cor_mat_round <- round(cor_mat, 3)

write.csv(
	cor_mat_round,
	"results/predictor_screening/solubility_descriptor_correlation_matrix.csv",
	row.names = TRUE
)

############################################################
# Linear relationship checks
############################################################

fit_saq_from_logkow_mpt <- lm(
	LogSaqd ~ logKowb + Mptc,
	data = df
)

fit_soc_from_logkow_mpt <- lm(
	LogSoce ~ logKowb + Mptc,
	data = df
)

fit_saq_from_gse <- lm(
	LogSaqd ~ GSE_logS,
	data = df
)

############################################################
# Base-R VIF function
############################################################

vif_one <- function(var, predictors, data) {
	others <- predictors[predictors != var]

	f <- as.formula(
		paste(var, "~", paste(others, collapse = " + "))
	)

	fit <- lm(f, data = data)
	r2 <- summary(fit)$r.squared

	1 / (1 - r2)
}

vif_predictors <- c(
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce"
)

vif_values <- sapply(
	vif_predictors,
	vif_one,
	predictors = vif_predictors,
	data = df
)

############################################################
# Summary table
############################################################

solubility_redundancy_table <- data.frame(
	Analysis = c(
		"Correlation between LogSaqd and GSE_logS",
		"Regression: LogSaqd ~ logKowb + Mptc",
		"Regression: LogSoce ~ logKowb + Mptc",
		"Regression: LogSaqd ~ GSE_logS",
		"VIF: logKowb",
		"VIF: Mptc",
		"VIF: LogSaqd",
		"VIF: LogSoce"
	),
	Statistic = c(
		"Pearson r",
		"R-squared",
		"R-squared",
		"R-squared",
		"VIF",
		"VIF",
		"VIF",
		"VIF"
	),
	Value = c(
		cor(df$LogSaqd, df$GSE_logS, use = "complete.obs"),
		summary(fit_saq_from_logkow_mpt)$r.squared,
		summary(fit_soc_from_logkow_mpt)$r.squared,
		summary(fit_saq_from_gse)$r.squared,
		vif_values["logKowb"],
		vif_values["Mptc"],
		vif_values["LogSaqd"],
		vif_values["LogSoce"]
	),
	Interpretation = c(
		"Agreement between dataset aqueous solubility and General Solubility Equation estimate",
		"Proportion of aqueous-solubility variation explained by logKow and melting point",
		"Proportion of organic-solubility variation explained by logKow and melting point",
		"Proportion of aqueous-solubility variation explained by GSE estimate",
		"Multicollinearity diagnostic",
		"Multicollinearity diagnostic",
		"Multicollinearity diagnostic",
		"Multicollinearity diagnostic"
	)
)

solubility_redundancy_table$value <- signif(
	solubility_redundancy_table$value,
	4
)

write.csv(
	solubility_redundancy_table,
	path_solubility_redundancy_table,
	row.names = FALSE
)

write.csv(
	solubility_redundancy_table,
	path_table_solubility_redundancy,
	row.names = FALSE
)

############################################################
# Figure S: GSE-predicted logS vs dataset LogSaqd
############################################################

plot_df <- df[
	complete.cases(df[, c("GSE_logS", "LogSaqd")]),
	c("GSE_logS", "LogSaqd")
]

fit_line <- lm(LogSaqd ~ GSE_logS, data = plot_df)

make_gse_logsaqd_plot <- function() {
	par(mar = c(5, 5, 3, 2))

	plot(
		plot_df$GSE_logS,
		plot_df$LogSaqd,
		pch = 16,
		xlab = "GSE-predicted log aqueous solubility",
		ylab = "Dataset aqueous solubility (LogSaqd)",
		main = "Relationship between GSE-predicted solubility and dataset LogSaqd"
	)

	abline(
		fit_line,
		lwd = 2
	)

	abline(
		0,
		1,
		lty = 2
	)

	legend(
		"topleft",
		legend = paste0(
			"r = ",
			round(cor(plot_df$LogSaqd, plot_df$GSE_logS), 3),
			"; R² = ",
			round(summary(fit_line)$r.squared, 3)
		),
		bty = "n"
	)
}

pdf(
	path_fig_gse_logs_vs_logsaqd_pdf,
	width = 7,
	height = 5
)

make_gse_logsaqd_plot()

dev.off()

png(
	path_fig_gse_logs_vs_logsaqd_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_gse_logsaqd_plot()

dev.off()

############################################################
# Figure S: descriptor correlation matrix as base-R heatmap
############################################################

pdf(
	path_fig_solubility_correlation,
	width = 6,
	height = 6
)

par(mar = c(7, 7, 3, 2))

image(
	1:ncol(cor_mat),
	1:nrow(cor_mat),
	cor_mat[nrow(cor_mat):1, ],
	axes = FALSE,
	xlab = "",
	ylab = "",
	main = "Correlation among solubility-related descriptors"
)

axis(
	1,
	at = 1:ncol(cor_mat),
	labels = colnames(cor_mat),
	las = 2
)

axis(
	2,
	at = 1:nrow(cor_mat),
	labels = rev(rownames(cor_mat)),
	las = 2
)

box()

for (i in seq_len(nrow(cor_mat))) {
	for (j in seq_len(ncol(cor_mat))) {
		text(
			j,
			nrow(cor_mat) - i + 1,
			labels = round(cor_mat[i, j], 2)
		)
	}
}

dev.off()

############################################################
# Save model summaries as text
############################################################

sink("results/predictor_screening/solubility_redundancy_model_summaries.txt")

cat("Correlation matrix\n")
print(cor_mat_round)

cat("\n\nModel: LogSaqd ~ logKowb + Mptc\n")
print(summary(fit_saq_from_logkow_mpt))

cat("\n\nModel: LogSoce ~ logKowb + Mptc\n")
print(summary(fit_soc_from_logkow_mpt))

cat("\n\nModel: LogSaqd ~ GSE_logS\n")
print(summary(fit_saq_from_gse))

cat("\n\nVIF values\n")
print(vif_values)

sink()

############################################################
# Console summary
############################################################

cat("Solubility redundancy check complete.\n")
cat("Correlation LogSaqd vs GSE_logS:",
	solubility_redundancy_table$value[
		solubility_redundancy_table$analysis ==
		"Correlation between LogSaqd and GSE_logS"
	],
	"\n"
)
cat("R2 LogSaqd ~ logKowb + Mptc:",
	solubility_redundancy_table$value[
		solubility_redundancy_table$analysis ==
		"Regression: LogSaqd ~ logKowb + Mptc"
	],
	"\n"
)