############################################################
# 08b_rdkit_descriptor_benchmark.R
# Benchmark RDKit descriptor models under LOCO-CV
############################################################

source("R/00_config.R")

############################################################
# Load datasets
############################################################

skin_df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

rdkit_as_reported <- read.csv(
	path_rdkit_descriptors_as_reported,
	stringsAsFactors = FALSE
)

rdkit_parent <- read.csv(
	path_rdkit_descriptors_parent,
	stringsAsFactors = FALSE
)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"

rdkit_molecular_predictors <- c(
	"rdkit_MolWt",
	"rdkit_MolLogP",
	"rdkit_MolMR",
	"rdkit_TPSA",
	"rdkit_NumHDonors",
	"rdkit_NumHAcceptors",
	"rdkit_NumRotatableBonds",
	"rdkit_RingCount",
	"rdkit_HeavyAtomCount",
	"rdkit_FractionCSP3",
	"rdkit_NumAromaticRings",
	"rdkit_FormalCharge"
)

experimental_predictors <- c(
	"Texpi",
	"Skin.thicknessj"
)

############################################################
# Check required columns
############################################################

required_skin_cols <- c(
	id_col,
	"CAS.No",
	"Compound",
	outcome_col,
	experimental_predictors
)

missing_skin_cols <- required_skin_cols[
	!(required_skin_cols %in% names(skin_df))
]

if (length(missing_skin_cols) > 0) {
	stop(
		paste(
			"Missing skin dataset columns:",
			paste(missing_skin_cols, collapse = ", ")
		)
	)
}

required_rdkit_cols <- c(
	id_col,
	"canonical_smiles",
	"representation",
	rdkit_molecular_predictors
)

missing_as_reported_cols <- required_rdkit_cols[
	!(required_rdkit_cols %in% names(rdkit_as_reported))
]

if (length(missing_as_reported_cols) > 0) {
	stop(
		paste(
			"Missing as-reported RDKit columns:",
			paste(missing_as_reported_cols, collapse = ", ")
		)
	)
}

missing_parent_cols <- required_rdkit_cols[
	!(required_rdkit_cols %in% names(rdkit_parent))
]

if (length(missing_parent_cols) > 0) {
	stop(
		paste(
			"Missing parent RDKit columns:",
			paste(missing_parent_cols, collapse = ", ")
		)
	)
}

############################################################
# Helper functions
############################################################

safe_rbind <- function(x) {
	x <- x[!sapply(x, is.null)]

	if (length(x) == 0) {
		return(data.frame())
	}

	all_cols <- unique(
		unlist(
			lapply(x, names)
		)
	)

	x2 <- lapply(
		x,
		function(d) {
			missing_cols <- setdiff(
				all_cols,
				names(d)
			)

			if (length(missing_cols) > 0) {
				for (m in missing_cols) {
					d[[m]] <- NA
				}
			}

			d[
				,
				all_cols,
				drop = FALSE
			]
		}
	)

	do.call(
		rbind,
		x2
	)
}

make_formula <- function(predictors) {
	as.formula(
		paste(
			outcome_col,
			"~",
			paste(predictors, collapse = " + ")
		)
	)
}

summarize_predictions <- function(pred_df,
								  model_name,
								  formula_text,
								  representation,
								  descriptor_set) {

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	data.frame(
		model = model_name,
		representation = representation,
		descriptor_set = descriptor_set,
		formula = formula_text,
		n_observations = nrow(pred_df),
		n_compounds = length(unique(pred_df$compound_id)),
		RMSE = rmse(pred_df$observed, pred_df$mu),
		MAE = mae(pred_df$observed, pred_df$mu),
		R2_pred = r2_pred(pred_df$observed, pred_df$mu),
		R = cor(
			pred_df$observed,
			pred_df$mu,
			use = "complete.obs"
		),
		mean_error = mean(pred_df$residual, na.rm = TRUE),
		median_abs_error = median(pred_df$abs_error, na.rm = TRUE),
		n_abs_error_gt_1 = sum(pred_df$abs_error > 1, na.rm = TRUE),
		proportion_abs_error_gt_1 = mean(pred_df$abs_error > 1, na.rm = TRUE),
		stringsAsFactors = FALSE
	)
}

############################################################
# Prepare merged dataset
############################################################

prepare_rdkit_model_dataset <- function(skin_df,
										rdkit_df,
										representation_name) {

	rdkit_keep_cols <- c(
		id_col,
		"canonical_smiles",
		"representation",
		rdkit_molecular_predictors
	)

	rdkit_df <- rdkit_df[
		,
		rdkit_keep_cols,
		drop = FALSE
	]

	model_df <- merge(
		skin_df,
		rdkit_df,
		by = id_col,
		all.x = FALSE,
		all.y = FALSE
	)

	model_df$rdkit_representation <- representation_name

	model_df
}

############################################################
# LOCO-CV linear model
############################################################

run_loco_cv_lm <- function(data,
						   model_name,
						   predictors,
						   representation,
						   descriptor_set,
						   id_col = "compound_id",
						   outcome_col = "logkpl") {

	model_formula <- make_formula(predictors)

	required_cols <- unique(c(
		id_col,
		"CAS.No",
		"Compound",
		"canonical_smiles",
		"rdkit_representation",
		outcome_col,
		predictors
	))

	model_data <- data[
		complete.cases(data[, required_cols, drop = FALSE]),
		,
		drop = FALSE
	]

	ids <- unique(model_data[[id_col]])

	pred_list <- vector(
		"list",
		length(ids)
	)

	for (i in seq_along(ids)) {
		id <- ids[i]

		training <- model_data[
			model_data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- model_data[
			model_data[[id_col]] == id,
			,
			drop = FALSE
		]

		fit <- lm(
			model_formula,
			data = training
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		pred_list[[i]] <- data.frame(
			model = model_name,
			representation = representation,
			descriptor_set = descriptor_set,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			canonical_smiles = test$canonical_smiles,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	summary_df <- summarize_predictions(
		pred_df = pred_df,
		model_name = model_name,
		formula_text = paste(deparse(model_formula), collapse = " "),
		representation = representation,
		descriptor_set = descriptor_set
	)

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# LOCO-CV random forest model
############################################################

run_loco_cv_rf <- function(data,
						   model_name,
						   predictors,
						   representation,
						   descriptor_set,
						   id_col = "compound_id",
						   outcome_col = "logkpl",
						   ntree = 1000,
						   nodesize = 5) {

	if (!requireNamespace("randomForest", quietly = TRUE)) {
		warning(
			"Package 'randomForest' is not installed. ",
			"Skipping ",
			model_name,
			"."
		)

		return(NULL)
	}

	model_formula <- make_formula(predictors)

	required_cols <- unique(c(
		id_col,
		"CAS.No",
		"Compound",
		"canonical_smiles",
		"rdkit_representation",
		outcome_col,
		predictors
	))

	model_data <- data[
		complete.cases(data[, required_cols, drop = FALSE]),
		,
		drop = FALSE
	]

	ids <- unique(model_data[[id_col]])

	pred_list <- vector(
		"list",
		length(ids)
	)

	mtry_value <- floor(
		sqrt(length(predictors))
	)

	if (mtry_value < 1) {
		mtry_value <- 1
	}

	for (i in seq_along(ids)) {
		id <- ids[i]

		cat(
			"RDKit RF LOCO-CV fold",
			i,
			"/",
			length(ids),
			"for",
			model_name,
			"\n"
		)

		training <- model_data[
			model_data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- model_data[
			model_data[[id_col]] == id,
			,
			drop = FALSE
		]

		set.seed(202604 + i)

		fit <- randomForest::randomForest(
			model_formula,
			data = training,
			ntree = ntree,
			mtry = mtry_value,
			nodesize = nodesize,
			importance = TRUE
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		pred_list[[i]] <- data.frame(
			model = model_name,
			representation = representation,
			descriptor_set = descriptor_set,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			canonical_smiles = test$canonical_smiles,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	summary_df <- summarize_predictions(
		pred_df = pred_df,
		model_name = model_name,
		formula_text = paste(deparse(model_formula), collapse = " "),
		representation = representation,
		descriptor_set = descriptor_set
	)

	summary_df$ntree <- ntree
	summary_df$mtry <- mtry_value
	summary_df$nodesize <- nodesize

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# Run one RDKit representation
############################################################

run_rdkit_representation_benchmarks <- function(skin_df,
												rdkit_df,
												representation_name) {

	model_df <- prepare_rdkit_model_dataset(
		skin_df = skin_df,
		rdkit_df = rdkit_df,
		representation_name = representation_name
	)

	cat("\n")
	cat("Running RDKit benchmarks for:", representation_name, "\n")
	cat("Merged rows:", nrow(model_df), "\n")
	cat("Merged compounds:", length(unique(model_df$compound_id)), "\n")

	prediction_list <- list()
	summary_list <- list()

	molecular_predictors <- rdkit_molecular_predictors

	molecular_plus_experimental_predictors <- c(
		rdkit_molecular_predictors,
		experimental_predictors
	)

	lm_molecular <- run_loco_cv_lm(
		data = model_df,
		model_name = paste0("rdkit_lm_molecular_only_", representation_name),
		predictors = molecular_predictors,
		representation = representation_name,
		descriptor_set = "molecular_only",
		id_col = id_col,
		outcome_col = outcome_col
	)

	prediction_list[[lm_molecular$summary$model]] <- lm_molecular$prediction
	summary_list[[lm_molecular$summary$model]] <- lm_molecular$summary

	lm_experimental <- run_loco_cv_lm(
		data = model_df,
		model_name = paste0("rdkit_lm_with_experimental_", representation_name),
		predictors = molecular_plus_experimental_predictors,
		representation = representation_name,
		descriptor_set = "molecular_plus_experimental",
		id_col = id_col,
		outcome_col = outcome_col
	)

	prediction_list[[lm_experimental$summary$model]] <- lm_experimental$prediction
	summary_list[[lm_experimental$summary$model]] <- lm_experimental$summary

	rf_molecular <- run_loco_cv_rf(
		data = model_df,
		model_name = paste0("rdkit_rf_molecular_only_", representation_name),
		predictors = molecular_predictors,
		representation = representation_name,
		descriptor_set = "molecular_only",
		id_col = id_col,
		outcome_col = outcome_col,
		ntree = 1000,
		nodesize = 5
	)

	if (!is.null(rf_molecular)) {
		prediction_list[[rf_molecular$summary$model]] <- rf_molecular$prediction
		summary_list[[rf_molecular$summary$model]] <- rf_molecular$summary
	}

	rf_experimental <- run_loco_cv_rf(
		data = model_df,
		model_name = paste0("rdkit_rf_with_experimental_", representation_name),
		predictors = molecular_plus_experimental_predictors,
		representation = representation_name,
		descriptor_set = "molecular_plus_experimental",
		id_col = id_col,
		outcome_col = outcome_col,
		ntree = 1000,
		nodesize = 5
	)

	if (!is.null(rf_experimental)) {
		prediction_list[[rf_experimental$summary$model]] <- rf_experimental$prediction
		summary_list[[rf_experimental$summary$model]] <- rf_experimental$summary
	}

	list(
		predictions = safe_rbind(prediction_list),
		summary = safe_rbind(summary_list)
	)
}

############################################################
# Run all RDKit benchmarks
############################################################

as_reported_result <- run_rdkit_representation_benchmarks(
	skin_df = skin_df,
	rdkit_df = rdkit_as_reported,
	representation_name = "as_reported"
)

parent_result <- run_rdkit_representation_benchmarks(
	skin_df = skin_df,
	rdkit_df = rdkit_parent,
	representation_name = "parent"
)

rdkit_predictions <- safe_rbind(
	list(
		as_reported_result$predictions,
		parent_result$predictions
	)
)

rdkit_summary <- safe_rbind(
	list(
		as_reported_result$summary,
		parent_result$summary
	)
)

rdkit_summary <- rdkit_summary[
	order(
		rdkit_summary$RMSE,
		rdkit_summary$MAE
	),
]

rdkit_summary$model_rank <- seq_len(
	nrow(rdkit_summary)
)

############################################################
# Save outputs
############################################################

write.csv(
	rdkit_predictions,
	path_rdkit_benchmark_predictions,
	row.names = FALSE
)

write.csv(
	rdkit_summary,
	path_rdkit_benchmark_summary,
	row.names = FALSE
)

write.csv(
	rdkit_summary,
	path_table_rdkit_benchmark_summary,
	row.names = FALSE
)

############################################################
# Figure: RDKit benchmark RMSE
############################################################

make_rdkit_benchmark_rmse_plot <- function() {

	plot_df <- rdkit_summary[
		order(
			rdkit_summary$RMSE,
			decreasing = TRUE
		),
	]

	par(
		mar = c(5, 18, 3, 2)
	)

	barplot(
		plot_df$RMSE,
		names.arg = plot_df$model,
		horiz = TRUE,
		las = 1,
		xlab = "LOCO-CV RMSE",
		main = "RDKit descriptor benchmark"
	)
}

pdf(
	path_fig_rdkit_benchmark_rmse_pdf,
	width = 10,
	height = 6
)

make_rdkit_benchmark_rmse_plot()

dev.off()

png(
	path_fig_rdkit_benchmark_rmse_png,
	width = 2500,
	height = 1500,
	res = 250
)

make_rdkit_benchmark_rmse_plot()

dev.off()

############################################################
# Figure: observed vs predicted for top RDKit models
############################################################

make_rdkit_observed_predicted_plot <- function() {

	top_models <- rdkit_summary$model[
		seq_len(
			min(4, nrow(rdkit_summary))
		)
	]

	plot_df <- rdkit_predictions[
		rdkit_predictions$model %in% top_models,
	]

	plot_range <- range(
		c(plot_df$observed, plot_df$mu),
		na.rm = TRUE
	)

	op <- par(
		mfrow = c(2, 2),
		mar = c(5, 5, 3, 2)
	)

	on.exit(
		par(op)
	)

	for (model_name in top_models) {

		model_df <- plot_df[
			plot_df$model == model_name,
		]

		plot(
			model_df$observed,
			model_df$mu,
			pch = 16,
			xlab = "Observed logKp",
			ylab = "LOCO-CV predicted logKp",
			xlim = plot_range,
			ylim = plot_range,
			main = model_name
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
				round(
					rmse(model_df$observed, model_df$mu),
					3
				)
			),
			bty = "n"
		)
	}
}

pdf(
	path_fig_rdkit_observed_predicted_pdf,
	width = 9,
	height = 9
)

make_rdkit_observed_predicted_plot()

dev.off()

png(
	path_fig_rdkit_observed_predicted_png,
	width = 2250,
	height = 2250,
	res = 250
)

make_rdkit_observed_predicted_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("\nRDKit descriptor benchmark complete.\n")
cat("Total prediction rows:", nrow(rdkit_predictions), "\n")
cat("Models evaluated:", nrow(rdkit_summary), "\n")
cat("\nModel summary:\n")

print(
	rdkit_summary[
		,
		c(
			"model",
			"representation",
			"descriptor_set",
			"n_observations",
			"n_compounds",
			"RMSE",
			"MAE",
			"R2_pred",
			"model_rank"
		)
	]
)