library("dplyr")
library("readxl")
library("knitr")
library("ggplot2")
library("bookdown")

############################################################
# 00_config.R
# Project-wide settings for skin permeability QSPR analysis
############################################################

############################################################
# File paths
############################################################



############################################################
# Predictor screening settings
############################################################

all_predictors <- c(
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

fixed_predictors <- c("MWa")

############################################################
# Functions
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

############################################################
# Benchmark model formulas and helper functions
############################################################

formula_null <- logkpl ~ 1

formula_potts_guy <- logkpl ~ MWa + logKowb

formula_linear_core <- logkpl ~ MWa + logKowb + Mptc + LogSaqd + LogSoce + Texpi



get_formula_predictors <- function(model_formula, outcome_col = "logkpl") {
	predictors <- all.vars(
		as.formula(model_formula)
	)

	predictors <- predictors[
		predictors != outcome_col
	]

	unique(predictors)
}