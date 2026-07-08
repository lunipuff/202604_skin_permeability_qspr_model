############################################################
# 05_loco_cv_error_diagnostics.R
# Diagnose prediction errors from selected LOCO-CV model
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset, selected model, and LOCO-CV predictions
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

est <- read.csv(
	path_loco_cv_final_predictions,
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
	"Compound",
	"logkpl",
	all_predictors
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

required_est_cols <- c(
	"df_row",
	"compound_id",
	"CAS.No",
	"logkpl",
	"mu"
)

missing_est_cols <- required_est_cols[
	!(required_est_cols %in% names(est))
]

if (length(missing_est_cols) > 0) {
	stop(
		paste(
			"Missing required columns in LOCO-CV prediction file:",
			paste(missing_est_cols, collapse = ", ")
		)
	)
}

############################################################
# Prepare diagnostic dataset
############################################################

df$Compound <- trimws(df$Compound)

diagnostic_df <- cbind(
	df,
	loco_mu = est$mu
)

diagnostic_df$loco_residual <- diagnostic_df$logkpl - diagnostic_df$loco_mu
diagnostic_df$loco_abs_error <- abs(diagnostic_df$loco_residual)

############################################################
# Calculate error summary
############################################################

error_summary <- data.frame(
	item = c(
		"selected_model",
		"n_observations",
		"n_compounds",
		"RMSE",
		"MAE",
		"R2_pred",
		"mean_error",
		"median_absolute_error",
		"n_abs_error_gt_1",
		"proportion_abs_error_gt_1"
	),
	value = c(
		selected_model,
		nrow(diagnostic_df),
		length(unique(diagnostic_df$compound_id)),
		rmse(diagnostic_df$logkpl, diagnostic_df$loco_mu),
		mae(diagnostic_df$logkpl, diagnostic_df$loco_mu),
		r2_pred(diagnostic_df$logkpl, diagnostic_df$loco_mu),
		mean(diagnostic_df$loco_residual, na.rm = TRUE),
		median(diagnostic_df$loco_abs_error, na.rm = TRUE),
		sum(diagnostic_df$loco_abs_error > 1, na.rm = TRUE),
		mean(diagnostic_df$loco_abs_error > 1, na.rm = TRUE)
	)
)

write.csv(
	error_summary,
	path_loco_cv_error_summary,
	row.names = FALSE
)

write.csv(
	diagnostic_df,
	path_loco_cv_error_diagnostic_dataset,
	row.names = FALSE
)

############################################################
# Figure: observed vs LOCO-CV predicted logKp
############################################################

make_observed_predicted_error_plot <- function() {

	plot_range <- range(
		c(diagnostic_df$logkpl, diagnostic_df$loco_mu),
		na.rm = TRUE
	)

	par(
		pty = "s",
		mar = c(5, 5, 3, 2)
	)

	plot(
		diagnostic_df$logkpl,
		diagnostic_df$loco_mu,
		pch = 16,
		xlab = "Observed logKp",
		ylab = "LOCO-CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "Observed vs predicted logKp"
	)

	abline(
		0,
		1,
		lty = 2,
		lwd = 2
	)

	legend(
		"topleft",
		legend = paste0(
			"RMSE = ",
			round(rmse(diagnostic_df$logkpl, diagnostic_df$loco_mu), 3),
			"; MAE = ",
			round(mae(diagnostic_df$logkpl, diagnostic_df$loco_mu), 3),
			"; R²pred = ",
			round(r2_pred(diagnostic_df$logkpl, diagnostic_df$loco_mu), 3)
		),
		bty = "n"
	)
}

pdf(
	path_fig_loco_cv_error_observed_predicted_pdf,
	width = 6,
	height = 6
)

make_observed_predicted_error_plot()

dev.off()

png(
	path_fig_loco_cv_error_observed_predicted_png,
	width = 1500,
	height = 1500,
	res = 250
)

make_observed_predicted_error_plot()

dev.off()

############################################################
# Figure: residual error by descriptor
############################################################

diagnostic_predictors <- all_predictors[
	all_predictors %in% names(diagnostic_df)
]

make_residual_by_descriptor_plot <- function() {

	n_plots <- length(diagnostic_predictors)
	n_row <- ceiling(n_plots / 4)

	par(
		mfrow = c(n_row, 4),
		mar = c(4, 5, 2, 1)
	)

	for (p in diagnostic_predictors) {

		x <- diagnostic_df[[p]]
		y <- diagnostic_df$loco_residual

		plot(
			x,
			y,
			pch = 16,
			xlab = p,
			ylab = "Prediction error"
		)

		abline(
			h = 0,
			lty = 2
		)

		if (sum(complete.cases(x, y)) > 5) {
			lines(
				lowess(x, y),
				lwd = 2
			)
		}
	}
}

pdf(
	path_fig_loco_cv_residual_by_descriptor_pdf,
	width = 12,
	height = 9
)

make_residual_by_descriptor_plot()

dev.off()

png(
	path_fig_loco_cv_residual_by_descriptor_png,
	width = 2400,
	height = 1800,
	res = 250
)

make_residual_by_descriptor_plot()

dev.off()

############################################################
# Figure: absolute error by descriptor
############################################################

make_absolute_error_by_descriptor_plot <- function() {

	n_plots <- length(diagnostic_predictors)
	n_row <- ceiling(n_plots / 4)

	par(
		mfrow = c(n_row, 4),
		mar = c(4, 5, 2, 1)
	)

	for (p in diagnostic_predictors) {

		x <- diagnostic_df[[p]]
		y <- diagnostic_df$loco_abs_error

		plot(
			x,
			y,
			pch = 16,
			xlab = p,
			ylab = "Absolute prediction error"
		)

		abline(
			h = 0,
			lty = 2
		)

		if (sum(complete.cases(x, y)) > 5) {
			lines(
				lowess(x, y),
				lwd = 2
			)
		}
	}
}

pdf(
	path_fig_loco_cv_absolute_error_by_descriptor_pdf,
	width = 12,
	height = 9
)

make_absolute_error_by_descriptor_plot()

dev.off()

png(
	path_fig_loco_cv_absolute_error_by_descriptor_png,
	width = 2400,
	height = 1800,
	res = 250
)

make_absolute_error_by_descriptor_plot()

dev.off()

############################################################
# Figure: observed prediction-error density
############################################################

make_error_density_plot <- function() {

	error <- diagnostic_df$loco_residual
	error <- error[!is.na(error)]

	error_sd <- sd(error, na.rm = TRUE)

	d_obs <- density(error)

	x <- seq(
		min(
			error,
			-4 * error_sd,
			na.rm = TRUE
		),
		max(
			error,
			4 * error_sd,
			na.rm = TRUE
		),
		length.out = 1000
	)

	y_expected <- dnorm(
		x,
		mean = 0,
		sd = error_sd
	)

	par(
		mar = c(5, 5, 3, 2)
	)

	plot(
		x,
		y_expected,
		type = "l",
		lwd = 2,
		ylim = range(
			c(y_expected, d_obs$y),
			na.rm = TRUE
		),
		xlab = "Prediction error",
		ylab = "Density",
		main = "Observed LOCO-CV prediction-error density"
	)

	lines(
		d_obs$x,
		d_obs$y,
		lwd = 2,
		lty = 2
	)

	abline(
		v = 0,
		lty = 3
	)

	legend(
		"topright",
		legend = c(
			"Normal density using observed error SD",
			"Observed error density"
		),
		lwd = c(2, 2),
		lty = c(1, 2),
		bty = "n"
	)
}

pdf(
	path_fig_loco_cv_error_density_pdf,
	width = 7,
	height = 5
)

make_error_density_plot()

dev.off()

png(
	path_fig_loco_cv_error_density_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_error_density_plot()

dev.off()

############################################################
# High-error and watchlist compounds
############################################################

watchlist <- c(
	"water",
	"sucrose",
	"mannitol",
	"tetraethylammonium bromide",
	"methotrexate",
	"hydrocortisone",
	"corticosterone",
	"estradiol"
)

diagnostic_df$point_class <- "Other compounds"

diagnostic_df$point_class[
	tolower(diagnostic_df$Compound) %in% watchlist
] <- "Watchlist compounds"

diagnostic_df$point_class[
	diagnostic_df$loco_abs_error > 1
] <- "High-error compounds"

point_col <- ifelse(
	diagnostic_df$point_class == "High-error compounds",
	"red",
	ifelse(
		diagnostic_df$point_class == "Watchlist compounds",
		"gold",
		"gray70"
	)
)

make_high_error_compound_plot <- function() {

	plot_range <- range(
		c(diagnostic_df$logkpl, diagnostic_df$loco_mu),
		na.rm = TRUE
	)

	par(
		xpd = TRUE,
		mar = c(5, 5, 4, 8),
		pty = "s"
	)

	plot(
		diagnostic_df$logkpl,
		diagnostic_df$loco_mu,
		pch = 16,
		col = point_col,
		xlab = "Observed logKp",
		ylab = "LOCO-CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "High-error and watchlist compounds"
	)

	abline(
		0,
		1,
		lty = 2,
		lwd = 2
	)

	legend(
		"bottomright",
		inset = c(-0.35, 0),
		legend = c(
			"Other compounds",
			"Watchlist compounds",
			"High-error compounds"
		),
		pch = 16,
		col = c(
			"gray70",
			"gold",
			"red"
		),
		bty = "n"
	)

	label_idx <- diagnostic_df$point_class != "Other compounds"

	text(
		diagnostic_df$logkpl[label_idx],
		diagnostic_df$loco_mu[label_idx],
		labels = diagnostic_df$Compound[label_idx],
		pos = 4,
		cex = 0.65
	)
}

pdf(
	path_fig_loco_cv_high_error_compounds_pdf,
	width = 8,
	height = 6
)

make_high_error_compound_plot()

dev.off()

png(
	path_fig_loco_cv_high_error_compounds_png,
	width = 2000,
	height = 1500,
	res = 250
)

make_high_error_compound_plot()

dev.off()

############################################################
# Save high-error compound table
############################################################

high_error_df <- diagnostic_df[
	diagnostic_df$loco_abs_error > 1,
]

high_error_df <- high_error_df[
	order(
		high_error_df$Compound,
		-high_error_df$loco_abs_error
	),
]

high_error_df <- high_error_df[
	,
	c(
		"Compound",
		"CAS.No",
		"MWa",
		"logKowb",
		"Mptc",
		"LogSaqd",
		"LogSoce",
		"Hdf",
		"Hag",
		"MVh",
		"Texpi",
		"Skin.thicknessj",
		"logkpl",
		"loco_mu",
		"loco_residual",
		"loco_abs_error"
	)
]

write.csv(
	high_error_df,
	path_loco_cv_high_error_compounds,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("LOCO-CV error diagnostics complete.\n")
cat("Selected model:\n")
cat(selected_model, "\n")
cat("RMSE:", rmse(diagnostic_df$logkpl, diagnostic_df$loco_mu), "\n")
cat("MAE:", mae(diagnostic_df$logkpl, diagnostic_df$loco_mu), "\n")
cat("R2_pred:", r2_pred(diagnostic_df$logkpl, diagnostic_df$loco_mu), "\n")
cat("High-error observations |error| > 1:", nrow(high_error_df), "\n")