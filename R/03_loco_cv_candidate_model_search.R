############################################################
# 03_loco_cv_candidate_model_search.R
# Search candidate interpretable models using leave-one-compound-out CV
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset and initial screening results
############################################################

df <- read.csv(path_cleaned_dataset)

linear_screening_models <- read.csv(
	path_predictor_screening_selected,
	stringsAsFactors = FALSE
)

############################################################
# Check required columns
############################################################

required_cols <- c(
	"compound_id",
	"CAS.No",
	"logkpl",
	all_predictors
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

if (!("formula" %in% names(linear_screening_models))) {
	stop("The initial predictor-screening file must contain a column named 'formula'.")
}

############################################################
# Helper function: extract model terms
############################################################

extract_terms <- function(model_formula) {
	rhs <- gsub(
		"^logkpl ~ ",
		"",
		model_formula
	)

	terms <- unlist(
		strsplit(rhs, " \\+ ")
	)

	terms
}

############################################################
# Helper function: get base predictor names from model terms
############################################################

get_base_terms <- function(model_formula) {
	terms <- extract_terms(model_formula)

	terms <- gsub("^I\\(", "", terms)
	terms <- gsub("\\^2\\)$", "", terms)
	terms <- gsub("^log\\(", "", terms)
	terms <- gsub("\\)$", "", terms)

	unique(terms)
}

############################################################
# Helper function: remove duplicate model formulas
############################################################

remove_duplicate_models <- function(models) {
	models <- models[!is.na(models)]
	models <- models[models != ""]

	if (length(models) == 0) {
		return(character(0))
	}

	model_keys <- character(length(models))

	for (i in seq_along(models)) {
		terms <- extract_terms(models[i])
		model_keys[i] <- paste(sort(terms), collapse = " + ")
	}

	models[!duplicated(model_keys)]
}

############################################################
# Helper function: count quadratic terms
############################################################

count_quadratic_terms <- function(model_formula) {
	sum(grepl("^I\\(.+\\^2\\)$", extract_terms(model_formula)))
}

############################################################
# Helper function: count unique base predictors
############################################################

count_base_predictors <- function(model_formula) {
	length(get_base_terms(model_formula))
}

############################################################
# Helper function: LOCO-CV for one model
############################################################

loco_cv_one_model <- function(model_formula, data) {

	est <- data.frame(
		df_row = seq_len(nrow(data)),
		compound_id = data$compound_id,
		CAS.No = data$CAS.No,
		logkpl = data$logkpl,
		mu = NA_real_
	)

	for (id in unique(data$compound_id)) {

		training <- data[data$compound_id != id, ]
		test <- data[data$compound_id == id, ]

		fit <- lm(
			as.formula(model_formula),
			data = training
		)

		est$mu[est$compound_id == id] <- predict(
			fit,
			newdata = test
		)
	}

	observed <- est$logkpl
	predicted <- est$mu

	data.frame(
		formula = model_formula,
		RMSE = rmse(observed, predicted),
		MAE = mae(observed, predicted),
		R2_pred = r2_pred(observed, predicted),
		R = cor(observed, predicted, use = "complete.obs"),
		n_base_predictors = count_base_predictors(model_formula),
		n_quadratic_terms = count_quadratic_terms(model_formula),
		n_terms = length(extract_terms(model_formula)),
		stringsAsFactors = FALSE
	)
}

############################################################
# Helper function: LOCO-CV for multiple models
############################################################

loco_cv_models <- function(models, data, label = "models") {

	models <- remove_duplicate_models(models)

	if (length(models) == 0) {
		return(data.frame())
	}

	cv_results <- data.frame()

	for (m in seq_along(models)) {

		cat(
			"LOCO-CV",
			label,
			m,
			"/",
			length(models),
			"\n"
		)

		model_result <- loco_cv_one_model(
			model_formula = models[m],
			data = data
		)

		cv_results <- rbind(
			cv_results,
			model_result
		)
	}

	cv_results <- cv_results[
		order(
			cv_results$RMSE,
			-cv_results$R2_pred,
			cv_results$n_terms
		),
	]

	cv_results$model_rank <- seq_len(nrow(cv_results))

	cv_results
}

############################################################
# Helper function: add quadratic terms
############################################################

add_quadratic_terms <- function(model_formula) {

	terms <- extract_terms(model_formula)
	base_terms <- get_base_terms(model_formula)

	candidate_terms <- base_terms[
		!(paste0("I(", base_terms, "^2)") %in% terms)
	]

	if (length(candidate_terms) == 0) {
		return(character(0))
	}

	new_models <- character(0)

	for (n_add in seq_along(candidate_terms)) {

		term_combinations <- combn(
			candidate_terms,
			n_add,
			simplify = FALSE
		)

		for (term_set in term_combinations) {

			new_terms <- terms

			for (term in term_set) {
				new_terms <- c(
					new_terms,
					paste0("I(", term, "^2)")
				)
			}

			new_model <- paste(
				"logkpl ~",
				paste(new_terms, collapse = " + ")
			)

			new_models <- c(new_models, new_model)
		}
	}

	remove_duplicate_models(new_models)
}

############################################################
# Helper function: add log-transformed terms
############################################################

add_log_terms <- function(model_formula) {

	terms <- extract_terms(model_formula)
	base_terms <- get_base_terms(model_formula)

	excluded_log_terms <- c(
		"logKowb",
		"LogSaqd",
		"LogSoce",
		"Hdf",
		"Hag"
	)

	candidate_terms <- base_terms[
		!(base_terms %in% excluded_log_terms)
	]

	candidate_terms <- candidate_terms[
		!(paste0("log(", candidate_terms, ")") %in% terms)
	]

	if (length(candidate_terms) == 0) {
		return(character(0))
	}

	new_models <- character(0)

	for (n_add in seq_along(candidate_terms)) {

		term_combinations <- combn(
			candidate_terms,
			n_add,
			simplify = FALSE
		)

		for (term_set in term_combinations) {

			new_terms <- terms

			for (term in term_set) {
				new_terms[new_terms == term] <- paste0("log(", term, ")")
			}

			new_model <- paste(
				"logkpl ~",
				paste(new_terms, collapse = " + ")
			)

			new_models <- c(new_models, new_model)
		}
	}

	remove_duplicate_models(new_models)
}

############################################################
# Helper function: add one omitted linear predictor
############################################################

add_one_linear_predictor <- function(model_formula, candidate_predictors) {

	current_predictors <- get_base_terms(model_formula)

	omitted_predictors <- candidate_predictors[
		!(candidate_predictors %in% current_predictors)
	]

	if (length(omitted_predictors) == 0) {
		return(character(0))
	}

	new_models <- paste(
		model_formula,
		omitted_predictors,
		sep = " + "
	)

	remove_duplicate_models(new_models)
}

############################################################
# Helper function: add one interaction between linear terms
############################################################

add_one_interaction <- function(model_formula) {

	terms <- extract_terms(model_formula)

	linear_terms <- terms[
		!grepl("^I\\(", terms) &
		!grepl("^log\\(", terms) &
		!grepl(":", terms)
	]

	if (length(linear_terms) < 2) {
		return(character(0))
	}

	interaction_pairs <- combn(
		linear_terms,
		2,
		simplify = FALSE
	)

	new_models <- character(0)

	for (pair in interaction_pairs) {

		interaction_term <- paste(
			pair,
			collapse = ":"
		)

		new_model <- paste(
			model_formula,
			interaction_term,
			sep = " + "
		)

		new_models <- c(new_models, new_model)
	}

	remove_duplicate_models(new_models)
}

############################################################
# Helper function: select competitive models
############################################################

select_competitive_models <- function(cv_results, rmse_tolerance = 0.01) {

	if (nrow(cv_results) == 0) {
		return(character(0))
	}

	best_rmse <- min(cv_results$RMSE, na.rm = TRUE)

	competitive <- cv_results[
		cv_results$RMSE <= best_rmse + rmse_tolerance,
	]

	competitive <- competitive[
		order(
			competitive$n_terms,
			competitive$RMSE
		),
	]

	competitive$formula
}

############################################################
# Stage 1: LOCO-CV of best apparent linear models from screening
############################################################

linear_models <- linear_screening_models$formula

cv_linear <- loco_cv_models(
	models = linear_models,
	data = df,
	label = "linear models"
)

cv_linear$model_stage <- "linear"

write.csv(
	cv_linear,
	path_loco_cv_linear,
	row.names = FALSE
)

############################################################
# Stage 2: Add quadratic terms
############################################################

quadratic_models <- character(0)

for (model in cv_linear$formula) {
	quadratic_models <- c(
		quadratic_models,
		add_quadratic_terms(model)
	)
}

quadratic_models <- remove_duplicate_models(quadratic_models)

cv_quadratic <- loco_cv_models(
	models = quadratic_models,
	data = df,
	label = "quadratic models"
)

cv_quadratic$model_stage <- "quadratic"

write.csv(
	cv_quadratic,
	path_loco_cv_quadratic,
	row.names = FALSE
)

############################################################
# Stage 3: Select linear/quadratic candidates for transformation search
############################################################

cv_linear_quadratic <- rbind(
	cv_linear,
	cv_quadratic
)

cv_linear_quadratic <- cv_linear_quadratic[
	order(
		cv_linear_quadratic$RMSE,
		cv_linear_quadratic$n_terms
	),
]

competitive_linear_quadratic_models <- select_competitive_models(
	cv_linear_quadratic,
	rmse_tolerance = 0.01
)

############################################################
# Stage 4: Add log-transformed terms
############################################################

log_models <- character(0)

for (model in competitive_linear_quadratic_models) {
	log_models <- c(
		log_models,
		add_log_terms(model)
	)
}

log_models <- remove_duplicate_models(log_models)

cv_log <- loco_cv_models(
	models = log_models,
	data = df,
	label = "log-transformed models"
)

cv_log$model_stage <- "log-transformed"

cv_with_log <- rbind(
	cv_linear_quadratic,
	cv_log
)

cv_with_log <- cv_with_log[
	order(
		cv_with_log$RMSE,
		cv_with_log$n_terms
	),
]

write.csv(
	cv_with_log,
	path_loco_cv_log,
	row.names = FALSE
)

############################################################
# Stage 5: Add one omitted linear predictor to competitive models
############################################################

competitive_log_models <- select_competitive_models(
	cv_with_log,
	rmse_tolerance = 0.01
)

additive_models <- character(0)

for (model in competitive_log_models) {
	additive_models <- c(
		additive_models,
		add_one_linear_predictor(
			model_formula = model,
			candidate_predictors = all_predictors
		)
	)
}

additive_models <- remove_duplicate_models(additive_models)

cv_additive <- loco_cv_models(
	models = additive_models,
	data = df,
	label = "additive models"
)

cv_additive$model_stage <- "additive"

cv_with_additive <- rbind(
	cv_with_log,
	cv_additive
)

cv_with_additive <- cv_with_additive[
	order(
		cv_with_additive$RMSE,
		cv_with_additive$n_terms
	),
]

write.csv(
	cv_with_additive,
	path_loco_cv_additive,
	row.names = FALSE
)

############################################################
# Stage 6: Add one interaction term to competitive models
############################################################

competitive_additive_models <- select_competitive_models(
	cv_with_additive,
	rmse_tolerance = 0.01
)

interaction_models <- character(0)

for (model in competitive_additive_models) {
	interaction_models <- c(
		interaction_models,
		add_one_interaction(model)
	)
}

interaction_models <- remove_duplicate_models(interaction_models)

cv_interaction <- loco_cv_models(
	models = interaction_models,
	data = df,
	label = "interaction models"
)

cv_interaction$model_stage <- "interaction"

cv_final <- rbind(
	cv_with_additive,
	cv_interaction
)

cv_final <- cv_final[
	order(
		cv_final$RMSE,
		cv_final$n_terms
	),
]

cv_final$model_rank <- seq_len(nrow(cv_final))

write.csv(
	cv_final,
	path_loco_cv_final_results,
	row.names = FALSE
)

############################################################
# Select final candidate model
############################################################

best_rmse <- min(cv_final$RMSE, na.rm = TRUE)

competitive_final <- cv_final[
	cv_final$RMSE <= best_rmse + 0.01,
]

competitive_final <- competitive_final[
	order(
		competitive_final$n_terms,
		competitive_final$RMSE
	),
]

selected_model <- competitive_final$formula[1]

writeLines(
	selected_model,
	path_loco_cv_selected_model
)

############################################################
# Generate LOCO-CV predictions for selected model
############################################################

est <- data.frame(
	df_row = seq_len(nrow(df)),
	compound_id = df$compound_id,
	CAS.No = df$CAS.No,
	logkpl = df$logkpl,
	mu = NA_real_
)

for (id in unique(df$compound_id)) {

	training <- df[df$compound_id != id, ]
	test <- df[df$compound_id == id, ]

	fit <- lm(
		as.formula(selected_model),
		data = training
	)

	est$mu[est$compound_id == id] <- predict(
		fit,
		newdata = test
	)
}

est$residual <- est$logkpl - est$mu
est$abs_error <- abs(est$residual)

write.csv(
	est,
	path_loco_cv_final_predictions,
	row.names = FALSE
)

############################################################
# Make observed-vs-predicted figure
############################################################

make_observed_predicted_plot <- function() {

	plot_range <- range(
		c(est$logkpl, est$mu),
		na.rm = TRUE
	)

	par(mar = c(5, 5, 3, 2))

	plot(
		est$logkpl,
		est$mu,
		pch = 16,
		xlab = "Observed logKp",
		ylab = "LOCO-CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "Selected model prediction performance"
	)

	abline(
		0,
		1,
		lty = 2
	)

	legend(
		"topleft",
		legend = paste0(
			"RMSE = ",
			round(rmse(est$logkpl, est$mu), 3),
			"; MAE = ",
			round(mae(est$logkpl, est$mu), 3),
			"; R²pred = ",
			round(r2_pred(est$logkpl, est$mu), 3)
		),
		bty = "n"
	)
}

pdf(
	path_fig_loco_cv_prediction_pdf,
	width = 7,
	height = 5
)

make_observed_predicted_plot()

dev.off()

png(
	path_fig_loco_cv_prediction_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_observed_predicted_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("LOCO-CV candidate model search complete.\n")
cat("Total final candidate models:", nrow(cv_final), "\n")
cat("Selected model:\n")
cat(selected_model, "\n")
cat("Selected model LOCO-CV RMSE:", rmse(est$logkpl, est$mu), "\n")
cat("Selected model LOCO-CV MAE:", mae(est$logkpl, est$mu), "\n")
cat("Selected model LOCO-CV R2_pred:", r2_pred(est$logkpl, est$mu), "\n")