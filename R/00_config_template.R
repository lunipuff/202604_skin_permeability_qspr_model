library("dplyr")

############################################################
# 00_config.R
# Project-wide settings for skin permeability QSPR analysis
############################################################

## Core column names
outcome_col <- "logkpl"
compound_col <- "compound_id"

## Main file paths
path_raw_data <- "data/raw/original_skin_permeability_dataset.csv"

path_clean_data <- "data/processed/skin_permeability_modeling_dataset.csv"
path_unique_compounds <- "data/processed/skin_permeability_unique_compounds.csv"

## Results paths
path_predictor_screening <- "results/predictor_screening/initial_predictor_screening.csv"
path_loco_predictions <- "results/cross_validation/loco_predictions.csv"
path_loco_performance <- "results/cross_validation/loco_model_performance.csv"
path_rowwise_predictions <- "results/cross_validation/rowwise_cv_predictions.csv"
path_rowwise_performance <- "results/cross_validation/rowwise_cv_performance.csv"

## Table paths
path_table_dataset_summary <- "tables/table1_dataset_summary.csv"
path_table_predictor_screening <- "tables/table2_predictor_screening.csv"
path_table_model_performance <- "tables/table4_model_performance.csv"
path_table_applicability_domain <- "tables/table7_applicability_domain.csv"

## Figure paths
path_fig_dataset_space <- "figures/figure2_dataset_space.pdf"
path_fig_predictor_screening <- "figures/figure3_predictor_screening.pdf"
path_fig_model_performance <- "figures/figure4_model_performance.pdf"
path_fig_applicability_domain <- "figures/figure6_applicability_domain.pdf"

############################################################
# Predictor sets
############################################################

## Full candidate predictor set available in the curated dataset
## Edit this after checking your actual column names.
candidate_predictors <- c(
  "MWa",
  "logKowb",
  "Mptc",
  "LogSaqd",
  "LogSoce",
  "Hdf",
  "Hag",
  "MVh",
  "Texpi",
  "Skin.thicknessj"
)

## Core predictors selected from best simple model
core_predictors <- c(
  "MWa",
  "logKowb",
  "Mptc",
  "LogSaqd",
  "Texpi"
)

############################################################
# Model formulas
############################################################

formula_dummy <- NULL

formula_potts_guy <- logkpl ~ MWa + logKowb

formula_linear_core <- logkpl ~ MWa + logKowb + Mptc + LogSaqd + Texpi

formula_quadratic_core <- logkpl ~ MWa + logKowb + I(logKowb^2) +
  Mptc + I(Mptc^2) +
  LogSaqd + I(LogSaqd^2) +
  Texpi + I(Texpi^2)

formula_extended_sensitivity <- logkpl ~ MWa + logKowb + I(logKowb^2) +
  Mptc + I(Mptc^2) +
  LogSaqd + I(LogSaqd^2) +
  Hdf + I(Hdf^2) +
  Hag +
  Texpi + I(Texpi^2) +
  Skin.thicknessj

## Named model list for model-comparison scripts
model_formulas <- list(
  potts_guy = formula_potts_guy,
  linear_core = formula_linear_core,
  quadratic_core = formula_quadratic_core,
  extended_sensitivity = formula_extended_sensitivity
)

############################################################
# Performance functions
############################################################

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

mae <- function(obs, pred) {
  mean(abs(obs - pred), na.rm = TRUE)
}

r2_pred <- function(obs, pred) {
  1 - sum((obs - pred)^2, na.rm = TRUE) /
    sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
}

bias <- function(obs, pred) {
  mean(pred - obs, na.rm = TRUE)
}