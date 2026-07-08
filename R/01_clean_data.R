############################################################
# 01_clean_data.R
# Clean human skin permeability dataset
############################################################

source("R/00_config.R")

############################################################
# Helper functions
############################################################

log_removed <- function(data, reason) {
	if (nrow(data) == 0) {
		return(data.frame())
	}

	data$removal_reason <- reason
	data
}

############################################################
# Read raw dataset
############################################################

rawxl <- read_excel(path_raw_excel)

write.csv(
	rawxl,
	path_interim_raw_csv,
	row.names = FALSE
)

raw <- read.csv(path_interim_raw_csv)

############################################################
# Initial formatting
############################################################

df <- raw
n_initial <- nrow(df)

# Standardize compound names
df$Compound <- tolower(df$Compound)
df$Compound <- trimws(df$Compound)

# Standardize CAS number format
df$CAS.No <- gsub("–", "-", df$CAS.No)
df$CAS.No <- gsub("--", "-", df$CAS.No)
df$CAS.No <- trimws(df$CAS.No)

# Arrange by compound name
df <- df[order(df$Compound), ]

############################################################
# Initialize removed-row log
############################################################

removed_rows <- data.frame()

############################################################
# Convert numeric columns
############################################################

numeric_cols <- c(
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
	"logkpl"
)

for (col in numeric_cols) {
	df[[col]] <- as.numeric(df[[col]])
}

############################################################
# Select relevant columns
############################################################

selected_cols <- c(
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
	"Reference"
)

df <- df[, selected_cols]

############################################################
# Remove rows without molecular volume
############################################################

missing_mvh <- is.na(df$MVh)

removed_rows <- rbind(
	removed_rows,
	log_removed(df[missing_mvh, ], "Missing molecular volume MVh")
)

df <- df[!missing_mvh, ]
n_after_missing_mvh <- nrow(df)

############################################################
# Remove rows without CAS number
############################################################

missing_cas <- is.na(df$CAS.No) | df$CAS.No == "" | df$CAS.No == "N/A"

removed_rows <- rbind(
	removed_rows,
	log_removed(df[missing_cas, ], "Missing CAS number")
)

df <- df[!missing_cas, ]
n_after_missing_cas <- nrow(df)

############################################################
# Remove known abnormal caffeine entry
############################################################

abnormal_caffeine <- df$CAS.No == "58-08-2" & df$logKowb == -0.63

removed_rows <- rbind(
	removed_rows,
	log_removed(df[abnormal_caffeine, ], "Removed abnormal caffeine logKowb entry")
)

df <- df[!abnormal_caffeine, ]
n_after_abnormal_caffeine <- nrow(df)

############################################################
# Collapse repeated entries from same compound, reference,
# temperature, skin thickness, and descriptor values
############################################################

condition_keys <- c(
	"CAS.No",
	"Compound",
	"Reference",
	"Texpi",
	"Skin.thicknessj"
)

descriptor_cols <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh"
)

merge_keys <- c(
	condition_keys,
	descriptor_cols
)

############################################################
# Check descriptor inconsistencies within compound-condition groups
############################################################

descriptor_inconsistency_log <- df %>%
	group_by(across(all_of(condition_keys))) %>%
	summarise(
		n_rows = n(),
		n_descriptor_profiles = n_distinct(
			paste(
				MWa,
				logKowb,
				Mptc,
				LogSaqd,
				LogSoce,
				Hdf,
				Hag,
				MVh,
				sep = "|"
			),
			na.rm = TRUE
		),
		across(
			all_of(descriptor_cols),
			~ n_distinct(.x, na.rm = TRUE),
			.names = "n_unique_{.col}"
		),
		.groups = "drop"
	)

descriptor_inconsistency_log <- descriptor_inconsistency_log[
	descriptor_inconsistency_log$n_descriptor_profiles > 1,
]

write.csv(
	descriptor_inconsistency_log,
	path_descriptor_inconsistency_log,
	row.names = FALSE
)

############################################################
# Merge only rows with identical compound-condition-descriptor profiles
############################################################

output <- df %>%
	group_by(across(all_of(merge_keys))) %>%
	summarise(
		logkpl = mean(logkpl, na.rm = TRUE),
		n_merged = n(),
		logkpl_sd = ifelse(n() > 1, sd(logkpl, na.rm = TRUE), NA_real_),
		.groups = "drop"
	)

############################################################
# Create compound ID
############################################################

output$compound_id <- output$CAS.No

############################################################
# Reorder columns
############################################################

output <- output[, c(
	"compound_id",
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
	"Reference",
	"n_merged",
	"logkpl_sd"
)]

############################################################
# Save cleaned dataset and removal log
############################################################

write.csv(
	output,
	path_cleaned_dataset,
	row.names = FALSE
)

write.csv(
	removed_rows,
	path_removed_rows,
	row.names = FALSE
)

############################################################
# Save cleaning summary
############################################################

n_collapsed <- sum(output$n_merged > 1, na.rm = TRUE)
n_records_collapsed <- sum(output$n_merged[output$n_merged > 1], na.rm = TRUE)

cleaning_summary <- data.frame(
	item = c(
		"initial_observations",
		"rows_after_missing_mvh_removal",
		"rows_after_missing_cas_removal",
		"rows_after_abnormal_caffeine_removal",
		"final_observations",
		"unique_compounds",
		"removed_rows",
		"collapsed_observations",
		"original_records_in_collapsed_observations",
		"compound_condition_groups_with_multiple_descriptor_profiles"
	),
	value = c(
		n_initial,
		n_after_missing_mvh,
		n_after_missing_cas,
		n_after_abnormal_caffeine,
		nrow(output),
		length(unique(output$compound_id)),
		nrow(removed_rows),
		n_collapsed,
		n_records_collapsed,
		nrow(descriptor_inconsistency_log)
	)
)

write.csv(
	cleaning_summary,
	path_cleaning_summary,
	row.names = FALSE
)

############################################################
# Summary messages
############################################################

cat("Cleaning complete.\n")
cat("Initial rows:", n_initial, "\n")
cat("Final rows:", nrow(output), "\n")
cat("Unique compounds:", length(unique(output$compound_id)), "\n")
cat("Removed rows:", nrow(removed_rows), "\n")
cat("Collapsed final observations:", n_collapsed, "\n")
cat("Original records in collapsed observations:", n_records_collapsed, "\n")
cat("Descriptor inconsistency groups:", nrow(descriptor_inconsistency_log), "\n")