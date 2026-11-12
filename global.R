# =============================================================================
# global.R
# Loaded ONCE at app startup.
# Contains: libraries, global constants, and utility functions shared
# across all modules.
# =============================================================================

library(shiny)
library(jsonlite)
library(readxl)
library(readr)
library(dplyr)
library(base64enc)

# -----------------------------------------------------------------------------
# LOGO: load BiomiX_logo3.png and encode as base64.
# This approach works regardless of how Shiny serves static files,
# because the image is embedded directly in the HTML.
# The file is looked for in www/ first, then in the app root directory.
# -----------------------------------------------------------------------------

.logo_path <- local({
  candidates <- c(
    file.path("www", "BiomiX_logo3.png"),
    "BiomiX_logo3.png"
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) found[1] else NULL
})

BIOMIX_LOGO_TAG <- if (!is.null(.logo_path)) {
  # Encode as base64 and embed directly in the img src
  img_data <- base64enc::base64encode(.logo_path)
  tags$img(
    src = paste0("data:image/png;base64,", img_data),
    height = "80px",
    style = "object-fit: contain;"
  )
} else {
  # Fallback: try the www/ URL path (standard Shiny static serving)
  tags$img(
    src = "BiomiX_logo3.png",
    height = "80px",
    style = "object-fit: contain;"
  )
}

BIOMIX_VERSION <- "2.5"

# Available omics data types (used in the 6-input dropdowns)
OMICS_TYPES <- c(
  "Undefined",
  "Transcriptomics",
  "Proteomics",
  "Methylomics",
  "Metabolomics"
)

# Integration methods available
INTEGRATION_METHODS <- c("MOFA", "Diablo", "SNF", "NEMO", "MintTea")

# Index labels for the 6 inputs (used in the final JSON)
INPUT_INDICES <- paste0("input", 1:6)

# Default placeholder for empty/unused fields in the JSON
EMPTY_VALUE <- "X"


# Supports .xlsx, .xls, .tsv, .csv, .txt
# Returns a dataframe or throws a user-readable error.
# -----------------------------------------------------------------------------

read_metadata <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))

  # The ID column is always read as text, so IDs such as "007" keep their
  # leading zeros and still match the sample names in the omics matrices.
  df <- tryCatch(
    {
      if (ext %in% c("xlsx", "xls")) {
        col_names <- names(read_excel(file_path, n_max = 0))
        read_excel(file_path,
                   col_types = ifelse(col_names == "ID", "text", "guess"))
      } else if (ext %in% c("tsv", "csv", "txt")) {
        read_delim(file_path, delim = NULL, show_col_types = FALSE,
                   col_types = cols(ID = col_character(), .default = col_guess()))
      } else {
        stop(paste(
          "Files with the extention ",
          ext,
          " are not supported. Please provide a csv, tsv, txt or excel file."
        ))
      }
    },
    error = function(e) {
      stop(paste("Error reading metadata file:", e$message))
    }
  )

  # ID column is mandatory: the pipeline matches it against the sample names
  if (!"ID" %in% colnames(df)) {
    stop("The metadata file does not contain the required 'ID' column.")
  }

  # CONDITION column is mandatory
  if (!"CONDITION" %in% colnames(df)) {
    stop("The metadata file does not contain the required 'CONDITION' column.")
  }

  # Warn if CONDITION is numeric
  if (is.numeric(df$CONDITION)) {
    warning(
      "The CONDITION column is numeric. Conditions should be character strings (e.g. 'Healthy', 'Disease')."
    )
  }

  df
}


# -----------------------------------------------------------------------------
# UTILITY: normalised metadata copy used by the pipeline.
# Whatever format the user selected (Excel, CSV, TSV, TXT), the app writes the
# already-validated table once as a TSV into the output folder. The JSON's
# METADATA_DIR points to this copy, so the pipeline scripts only ever read TSV:
#   vroom(DIR_METADATA, delim = "\t", col_types = cols(ID = col_character()))
# The user's original file is never modified.
# -----------------------------------------------------------------------------

METADATA_TSV_NAME <- "BiomiX_metadata.tsv"

metadata_tsv_path <- function(output_dir) {
  if (nchar(output_dir) == 0) return("")
  file.path(output_dir, METADATA_TSV_NAME)
}

# Returns the path of the written TSV, or throws a user-readable error.
# Callers (Run / Download JSON in app.R) catch the error and show it.
write_metadata_tsv <- function(df, output_dir) {
  path <- metadata_tsv_path(output_dir)
  if (nchar(path) == 0) {
    stop("No output folder selected.", call. = FALSE)
  }
  if (is.null(df) || nrow(df) == 0) {
    stop("No metadata loaded.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop(sprintf("The output folder does not exist: %s", output_dir), call. = FALSE)
  }

  tryCatch(
    readr::write_tsv(df, path, na = "NA"),
    error = function(e) {
      # Do not leave a half-written file behind that the pipeline could read
      if (file.exists(path) && !dir.exists(path)) unlink(path)
      stop(sprintf("Could not save the metadata as TSV to '%s': %s",
                   path, conditionMessage(e)), call. = FALSE)
    }
  )

  invisible(path)
}


# -----------------------------------------------------------------------------
# UTILITY: extract unique conditions from metadata dataframe
# -----------------------------------------------------------------------------

get_conditions <- function(df) {
  as.character(unique(df$CONDITION))
}


# -----------------------------------------------------------------------------
# UTILITY: build the COMMANDS list (6 inputs) for the final JSON.
# Arguments:
#   inputs_data  = list collected from mod_inputs
# -----------------------------------------------------------------------------

build_commands <- function(inputs_data) {
  lapply(seq_len(6), function(i) {
    d <- inputs_data[[i]]
    list(
      INDEX = INPUT_INDICES[i],
      ANALYSIS = d$analysis,
      DATA_TYPE = d$data_type,
      INTEGRATION = d$integration,
      LABEL = if (nchar(trimws(d$label)) == 0) NA_character_ else d$label,
      SELECTION = d$selection,
      DIRECTORIES = if (nchar(trimws(d$path)) == 0) EMPTY_VALUE else d$path,
      PREVIEW = d$preview
    )
  })
}


# -----------------------------------------------------------------------------
# UTILITY: build the COMMANDS_MOFA list for the final JSON.
# Arguments:
#   integration_data = list collected from mod_integration
# -----------------------------------------------------------------------------

build_commands_mofa <- function(integration_data) {
  # Convert method name to the format expected by BiomiX_BETA.r
  # e.g. "MOFA" -> "MOFA_INTEGRATION", "MintTea" -> "MINTEA_INTEGRATION"
  integration_type <- if (
    integration_data$integration_type == "NO_INTEGRATION"
  ) {
    "NO_INTEGRATION"
  } else {
    paste0(toupper(integration_data$integration_type), "_INTEGRATION")
  }

  list(
    list(nms = "Integration_type", V2 = integration_type),
    list(nms = "Data_Integration", V2 = integration_data$data_integration),
    list(
      nms = "Number_of_factors_to_calculate",
      V2 = as.character(integration_data$n_factors)
    ),
    list(
      nms = "Factor_to_explore",
      V2 = as.character(integration_data$factor_explore)
    ),
    list(
      nms = "Omics_overlay",
      V2 = as.character(integration_data$omics_overlay)
    )
  )
}


# -----------------------------------------------------------------------------
# UTILITY: assemble and write the final combined JSON to disk.
# Arguments:
#   commands          = output of build_commands()
#   commands_mofa     = output of build_commands_mofa()
#   commands_advanced = list collected from all mod_adv_* modules
#   group1, group2    = condition / control strings
#   base_dir          = BiomiX base directory
#   metadata_dir      = path to the metadata file
#   output_dir        = results output directory
#   out_path          = destination path for the JSON file
# -----------------------------------------------------------------------------

write_combined_json <- function(
  commands,
  commands_mofa,
  commands_advanced,
  group1,
  group2,
  base_dir,
  metadata_dir,
  output_dir,
  out_path = "COMBINED_COMMANDS.json"
) {
  combined <- list(
    COMMANDS = commands,
    COMMANDS_MOFA = commands_mofa,
    COMMANDS_ADVANCED = commands_advanced,
    COMMAND_LINE_ARGS = list(
      GROUP_1 = group1,
      GROUP_2 = group2,
      BASE_DIR = base_dir
    ),
    DIRECTORY_INFO = list(
      METADATA_DIR = metadata_dir,
      OUTPUT_DIR = output_dir
    )
  )

  write(
    jsonlite::toJSON(combined, pretty = TRUE, auto_unbox = TRUE, na = "null"),
    out_path
  )

  invisible(out_path)
}


# -----------------------------------------------------------------------------
# UTILITY: build the COMMANDS_ADVANCED list (3 rows) for the final JSON.
# The pipeline scripts read each option by row index, e.g.
#   COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[3]
# so the row in which each value is placed matters.
# Arguments: the lists returned by the advanced modules
#   ag = general, met = metabolomics, mf = MOFA/Diablo, sn = SNF/NEMO,
#   mt = MintTea, fil = metadata filter
# -----------------------------------------------------------------------------

build_commands_advanced <- function(ag, met, mf, sn, mt, fil) {
  list(
    # Row 1: user-configured values
    list(
      ADVANCED_OPTION_TRASCRIPTOMICS                      = as.character(ag$tx_fc),
      ADVANCED_OPTION_SUBGROUPING                         = as.character(ag$sub_zscore2),
      ADVANCED_OPTION_METABOLOMICS                        = ag$met_fc,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = met$ms1_ion_mode,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = met$ms1_pos,
      ADVANCED_OPTION_METHYLOMICS                         = as.character(ag$meth_fc),
      ADVANCED_OPTION_MOFA                                = mf$mofa_max_iter,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = as.character(met$ms2_ppm),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = met$annot_mode,
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = mf$biblio_type,
      ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = as.character(mf$pathway_padj),
      ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = mf$clinical_numerical,
      ADVANCED_OPTION_METADATA_FILTERING_1                = fil[[1]]$column,
      ADVANCED_OPTION_METADATA_FILTERING_2                = fil[[2]]$column,
      ADVANCED_OPTION_METADATA_FILTERING_3                = fil[[3]]$column,
      ADVANCED_OPTION_METADATA_FILTERING_4                = fil[[4]]$column,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = met$ms1_input_idx,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = met$ms1_annot_file,
      ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = met$ms2_dir,
      ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = mf$clinic_file_numerical,
      ADVANCED_OPTION_CLUSTERING_OPTIONS                  = ag$clust_distance,
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = as.character(mf$biblio_top_contrib),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = met$ms2_pos,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = met$ms2_ion_mode,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = met$ms2_annot_idx,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = met$ms2_annot_file,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = met$ms2_ms1_db,
      ADVANCED_OPTION_SNF_OPTIONS                         = sn$snf_neighbors,
      ADVANCED_OPTION_NEMO_OPTIONS                        = as.character(sn$nemo_neighbors),
      ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = sn$snf_enrich_cols,
      ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = sn$snf_cluster_range,
      ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = mf$diablo_design_file,
      ADVANCED_OPTION_DIABLO_OPTIONS                      = mf$diablo_feature_sel,
      ADVANCED_OPTION_MINTTEA_OPTIONS                     = as.character(mt$n_repeats),
      ADVANCED_OPTION_MINTTEA_FEATURES                    = as.character(mt$keep_x)
    ),
    # Row 2: minimum/alternative thresholds (fixed constants from global.R)
    #COMMANDS_ADVANCED_ROW2,
    list(
      ADVANCED_OPTION_TRASCRIPTOMICS                      = as.character(ag$tx_padj),
      ADVANCED_OPTION_SUBGROUPING                         = as.character(ag$sub_zscore1),
      ADVANCED_OPTION_METABOLOMICS                        = ag$met_padj,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = met$ms1_input_idx,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = met$ms1_neg,
      ADVANCED_OPTION_METHYLOMICS                         = as.character(ag$meth_padj),
      ADVANCED_OPTION_MOFA                                = mf$mofa_convergence,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = as.character(met$ms2_rt),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = met$annot_compound_col,
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = mf$biblio_n_articles,
      ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = as.character(mf$pathway_n),
      ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = mf$clinical_binary,
      ADVANCED_OPTION_METADATA_FILTERING_1                = fil[[1]]$data_type,
      ADVANCED_OPTION_METADATA_FILTERING_2                = fil[[2]]$data_type,
      ADVANCED_OPTION_METADATA_FILTERING_3                = fil[[3]]$data_type,
      ADVANCED_OPTION_METADATA_FILTERING_4                = fil[[4]]$data_type,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = as.character("X"),
      ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = mf$clinic_file_binary,
      ADVANCED_OPTION_CLUSTERING_OPTIONS                  = ag$clust_method,
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = as.character(mf$biblio_n_keywords),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = met$ms2_neg,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = met$ms2_ms1_ppm, #Tolerance MS1-MS2 peak matching
      ADVANCED_OPTION_SNF_OPTIONS                         = sn$snf_variance,
      ADVANCED_OPTION_NEMO_OPTIONS                        = as.character(sn$apply_scaling_SNF),
      ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = sn$snf_surv_cols,
      ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = sn$snf_var_heatmap,
      ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = as.character("X"),
      ADVANCED_OPTION_DIABLO_OPTIONS                      = mf$diablo_n_iter,
      ADVANCED_OPTION_MINTTEA_OPTIONS                     = as.character(mt$n_folds),
      ADVANCED_OPTION_MINTTEA_FEATURES                    = as.character(mt$design_val)
    ),
    # Row 3: maximum values / extra options (fixed constants from global.R)
    #COMMANDS_ADVANCED_ROW3
    list(
      ADVANCED_OPTION_TRASCRIPTOMICS                      = as.character(ag$tx_panel),
      ADVANCED_OPTION_SUBGROUPING                         = ag$sub_enable,  # YES/NO: remove positive controls (read as [3] by Transcriptomics)
      ADVANCED_OPTION_METABOLOMICS                        = ag$cpu_threads,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = met$ms1_ppm,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = met$ms1_db,
      ADVANCED_OPTION_METHYLOMICS                         = as.character(ag$meth_array),
      ADVANCED_OPTION_MOFA                                = mf$mofa_threshold,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = as.character(met$ms2_ms1_db_priority),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = as.character("X"),
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = ag$mofa_features,
      ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = as.character("X"),
      ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = mf$clinic_file,
      ADVANCED_OPTION_METADATA_FILTERING_1                = fil[[1]]$threshold,
      ADVANCED_OPTION_METADATA_FILTERING_2                = fil[[2]]$threshold,
      ADVANCED_OPTION_METADATA_FILTERING_3                = fil[[3]]$threshold,
      ADVANCED_OPTION_METADATA_FILTERING_4                = fil[[4]]$threshold,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = as.character("X"),
      ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = as.character("X"),
      ADVANCED_OPTION_CLUSTERING_OPTIONS                  = ag$sub_heatmap,
      ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = as.character(mf$biblio_top_contrib),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = met$ms2_col_type,
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = as.character("X"),
      ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = as.character("X"),
      ADVANCED_OPTION_SNF_OPTIONS                         = sn$snf_iterations,
      ADVANCED_OPTION_NEMO_OPTIONS                        = sn$only_common_samples,
      ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = sn$snf_condition_col,
      ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = sn$snf_max_features,
      ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = as.character("X"),
      ADVANCED_OPTION_DIABLO_OPTIONS                      = as.character("X"),
      ADVANCED_OPTION_MINTTEA_OPTIONS                     = as.character(mt$edge_threshold),
      ADVANCED_OPTION_MINTTEA_FEATURES                    = as.character("X")
    )
  )
}
