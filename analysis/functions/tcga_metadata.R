tcga_patient_barcode <- function(barcodes) {
  substr(barcodes, 1, 12)
}

normalize_stage <- function(stage) {
  stage <- toupper(trimws(as.character(stage)))
  stage <- gsub("[ABC]$", "", stage)
  dplyr::case_when(
    grepl("STAGE IV", stage) ~ "Stage IV",
    grepl("STAGE III", stage) ~ "Stage III",
    grepl("STAGE II", stage) ~ "Stage II",
    grepl("STAGE I", stage) ~ "Stage I",
    TRUE ~ NA_character_
  )
}

normalize_grade <- function(grade) {
  grade <- toupper(trimws(as.character(grade)))
  dplyr::case_when(
    grade %in% c("G1", "GRADE 1") ~ "G1",
    grade %in% c("G2", "GRADE 2") ~ "G2",
    grade %in% c("G3", "GRADE 3") ~ "G3",
    grade %in% c("G4", "GRADE 4") ~ "G4",
    TRUE ~ NA_character_
  )
}

get_optional_column <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) return(rep(NA, nrow(data)))
  data[[found[1]]]
}

derive_os <- function(clinical) {
  death <- get_optional_column(clinical, "days_to_death")
  follow_up <- get_optional_column(clinical, c("days_to_last_follow_up", "days_to_last_known_alive"))
  vital <- get_optional_column(clinical, "vital_status")
  vital_clean <- tolower(trimws(as.character(vital)))
  if (any(!is.na(vital_clean) & !vital_clean %in% c("alive", "dead"))) {
    stop("Unexpected TCGA vital_status value; survival events cannot be coded safely.")
  }
  death_clean <- suppressWarnings(as.numeric(death))
  follow_up_clean <- suppressWarnings(as.numeric(follow_up))
  if (any(vital_clean == "dead" & is.na(death_clean) & !is.na(follow_up_clean), na.rm = TRUE)) {
    stop("A deceased TCGA patient lacks a death time; refusing follow-up as death time.")
  }

  clinical |>
    dplyr::mutate(
      days_to_death_clean = death_clean,
      days_to_last_follow_up_clean = follow_up_clean,
      os_time = dplyr::if_else(vital_clean == "dead", death_clean, follow_up_clean),
      os_event = dplyr::case_when(vital_clean == "dead" ~ 1L, vital_clean == "alive" ~ 0L, TRUE ~ NA_integer_)
    )
}

add_clean_stage_grade <- function(clinical) {
  stage <- get_optional_column(
    clinical,
    c("ajcc_pathologic_stage", "tumor_stage", "pathologic_stage", "stage")
  )
  grade <- get_optional_column(
    clinical,
    c("grade", "tumor_grade", "neoplasm_histologic_grade")
  )

  clinical |>
    dplyr::mutate(
      stage_clean = normalize_stage(stage),
      grade_clean = normalize_grade(grade)
    )
}
