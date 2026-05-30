tcga_patient_barcode <- function(barcodes) {
  substr(barcodes, 1, 12)
}

normalize_stage <- function(stage) {
  stage <- toupper(trimws(as.character(stage)))
  dplyr::case_when(
    grepl("STAGE I[^V]|STAGE I$", stage) ~ "Stage I",
    grepl("STAGE II[^I]|STAGE II$", stage) ~ "Stage II",
    grepl("STAGE III", stage) ~ "Stage III",
    grepl("STAGE IV", stage) ~ "Stage IV",
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

  clinical |>
    dplyr::mutate(
      days_to_death_clean = suppressWarnings(as.numeric(death)),
      days_to_last_follow_up_clean = suppressWarnings(as.numeric(follow_up)),
      os_time = dplyr::coalesce(.data$days_to_death_clean, .data$days_to_last_follow_up_clean),
      os_event = dplyr::if_else(tolower(as.character(vital)) == "dead", 1L, 0L)
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
