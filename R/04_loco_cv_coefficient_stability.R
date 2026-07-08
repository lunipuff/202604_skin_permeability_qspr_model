############################################################
# 04_loco_cv_coefficient_stability.R
# Assess coefficient stability across LOCO-CV training folds
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset and LOCO-CV model results
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

cv_results <- read.csv(
	path_loco_cv_final_results,
	stringsAsFactors = FALSE
)

selected_model <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model <- selected_model[1]

############################################################
# Check required columns
############################################################

required_df_cols <- c(
	"compound_id",
	"CAS.No",
	"logkpl"
)

missing_df_cols <- required_df_cols[
	!(required_df_cols %in% names(df))
]

if (length(missing_df_cols) > 0) {
	stop(
		paste(
			"Missing required columns in cleaned dataset:",
			paste(missing_df_cols, collapse = ", ")
		)
	)
}

if (!("formula" %in% names(cv_results))) {
	stop("The LOCO-CV results file must contain a column named 'formula'.")
}

############################################################
# Select models for coefficient-stability analysis
############################################################

# By default, analyze all final candidate models from the previous step.
# If this becomes slow, replace this with:
# model_formulas <- cv_results$formula[1:min(50, nrow(cv_results))]

model_formulas <- cv_results$formula

model_formulas <- model_formulas[
	!is.na(model_formulas) &
		model_formulas != ""
]

model_formulas <- unique(model_formulas)

############################################################
# Helper function: extract coefficient names from one model
############################################################

get_model_coefficient_names <- function(model_formula, data) {
	fit <- lm(
		as.formula(model_formula),
		data = data
	)

	names(coef(fit))
}

############################################################
# Helper function: coefficient stability for one model
############################################################

coefficient_stability_one_model <- function(model_formula, data) {

	ids <- unique(data$compound_id)

	full_fit <- lm(
		as.formula(model_formula),
		data = data
	)

	coef_names <- names(coef(full_fit))

	coef_by_fold <- as.data.frame(
		matrix(
			NA_real_,
			nrow = length(ids),
			ncol = length(coef_names)
		)
	)

	colnames(coef_by_fold) <- coef_names

	coef_by_fold$left_out_compound_id <- ids

	for (i in seq_along(ids)) {

		training <- data[
			data$compound_id != ids[i],
		]

		fit <- lm(
			as.formula(model_formula),
			data = training
		)

		fold_coef <- coef(fit)

		coef_by_fold[
			i,
			names(fold_coef)
		] <- fold_coef
	}

	coef_matrix <- coef_by_fold[
		,
		coef_names,
		drop = FALSE
	]

	coef_mean <- apply(
		coef_matrix,
		2,
		mean,
		na.rm = TRUE
	)

	coef_sd <- apply(
		coef_matrix,
		2,
		sd,
		na.rm = TRUE
	)

	coef_rsd <- coef_sd / abs(coef_mean) * 100

	coef_rsd[
		is.infinite(coef_rsd)
	] <- NA_real_

	coef_summary <- data.frame(
		formula = model_formula,
		term = coef_names,
		mean = as.numeric(coef_mean),
		sd = as.numeric(coef_sd),
		rsd_percent = as.numeric(coef_rsd),
		n_folds = length(ids),
		n_nonmissing = apply(
			coef_matrix,
			2,
			function(x) sum(!is.na(x))
		),
		stringsAsFactors = FALSE
	)

	coef_long <- data.frame()

	for (term in coef_names) {
		temp <- data.frame(
			formula = model_formula,
			left_out_compound_id = coef_by_fold$left_out_compound_id,
			term = term,
			coefficient = coef_by_fold[[term]],
			stringsAsFactors = FALSE
		)

		coef_long <- rbind(
			coef_long,
			temp
		)
	}

	list(
		summary = coef_summary,
		long = coef_long
	)
}

############################################################
# Run coefficient-stability analysis
############################################################

coefficient_summary_all <- data.frame()
coefficient_long_all <- data.frame()

for (m in seq_along(model_formulas)) {

	cat(
		"Coefficient stability",
		m,
		"/",
		length(model_formulas),
		"\n"
	)

	model_result <- coefficient_stability_one_model(
		model_formula = model_formulas[m],
		data = df
	)

	coefficient_summary_all <- rbind(
		coefficient_summary_all,
		model_result$summary
	)

	coefficient_long_all <- rbind(
		coefficient_long_all,
		model_result$long
	)
}

############################################################
# Create wide summary tables
############################################################

all_terms <- sort(
	unique(coefficient_summary_all$term)
)

make_wide_table <- function(summary_df, value_col) {

	wide_df <- data.frame(
		formula = model_formulas,
		stringsAsFactors = FALSE
	)

	for (term in all_terms) {
		wide_df[[term]] <- NA_real_
	}

	for (i in seq_len(nrow(summary_df))) {
		model_formula <- summary_df$formula[i]
		term <- summary_df$term[i]

		wide_df[
			wide_df$formula == model_formula,
			term
		] <- summary_df[[value_col]][i]
	}

	wide_df
}

coefficient_mean_wide <- make_wide_table(
	coefficient_summary_all,
	"mean"
)

coefficient_sd_wide <- make_wide_table(
	coefficient_summary_all,
	"sd"
)

coefficient_rsd_wide <- make_wide_table(
	coefficient_summary_all,
	"rsd_percent"
)

############################################################
# Save coefficient-stability outputs
############################################################

write.csv(
	coefficient_mean_wide,
	path_loco_cv_coefficient_mean,
	row.names = FALSE
)

write.csv(
	coefficient_sd_wide,
	path_loco_cv_coefficient_sd,
	row.names = FALSE
)

write.csv(
	coefficient_rsd_wide,
	path_loco_cv_coefficient_rsd,
	row.names = FALSE
)

write.csv(
	coefficient_long_all,
	path_loco_cv_coefficient_long,
	row.names = FALSE
)

############################################################
# Extract selected-model coefficient summary
############################################################

selected_model_coefficient_summary <- coefficient_summary_all[
	coefficient_summary_all$formula == selected_model,
]

selected_model_coefficient_summary <- selected_model_coefficient_summary[
	order(
		selected_model_coefficient_summary$term
	),
]

write.csv(
	selected_model_coefficient_summary,
	path_loco_cv_selected_model_coefficient_summary,
	row.names = FALSE
)

############################################################
# Figure: selected-model coefficient RSD
############################################################

make_selected_model_rsd_plot <- function() {

	plot_df <- selected_model_coefficient_summary

	plot_df <- plot_df[
		plot_df$term != "(Intercept)",
	]

	plot_df <- plot_df[
		order(plot_df$rsd_percent, decreasing = TRUE),
	]

	par(mar = c(8, 5, 3, 2))

	barplot(
		plot_df$rsd_percent,
		names.arg = plot_df$term,
		las = 2,
		ylab = "Relative standard deviation (%)",
		main = "Selected model coefficient stability"
	)

	abline(
		h = 100,
		lty = 2
	)
}

pdf(
	path_fig_selected_model_coefficient_rsd_pdf,
	width = 7,
	height = 5
)

make_selected_model_rsd_plot()

dev.off()

png(
	path_fig_selected_model_coefficient_rsd_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_selected_model_rsd_plot()

dev.off()

############################################################
# Figure: selected-model coefficient distributions
############################################################

make_selected_model_coefficient_boxplot <- function() {

	plot_df <- coefficient_long_all[
		coefficient_long_all$formula == selected_model &
			coefficient_long_all$term != "(Intercept)",
	]

	term_order <- selected_model_coefficient_summary$term[
		selected_model_coefficient_summary$term != "(Intercept)"
	]

	plot_df$term <- factor(
		plot_df$term,
		levels = term_order
	)

	par(mar = c(8, 5, 3, 2))

	boxplot(
		coefficient ~ term,
		data = plot_df,
		las = 2,
		xlab = "",
		ylab = "Coefficient estimate",
		main = "Selected model coefficient estimates across LOCO-CV folds"
	)

	abline(
		h = 0,
		lty = 2
	)
}

pdf(
	path_fig_selected_model_coefficient_boxplot_pdf,
	width = 7,
	height = 5
)

make_selected_model_coefficient_boxplot()

dev.off()

png(
	path_fig_selected_model_coefficient_boxplot_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_selected_model_coefficient_boxplot()

dev.off()

############################################################
# Console summary
############################################################

cat("LOCO-CV coefficient-stability analysis complete.\n")
cat("Models analyzed:", length(model_formulas), "\n")
cat("Selected model:\n")
cat(selected_model, "\n")
cat("Selected-model coefficient summary saved to:\n")
cat(path_loco_cv_selected_model_coefficient_summary, "\n")