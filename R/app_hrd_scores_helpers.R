library(data.table)

empty_hrd_scores_data <- function() {
  data.table(
    sample_id = character(),
    cancer_type = character(),
    HRDsum = numeric(),
    HRD_LOH = numeric(),
    LST = numeric(),
    TAI = numeric(),
    purity = numeric(),
    ploidy = numeric(),
    epi_HRD = numeric()
  )
}

prepare_hrd_scores_data <- function(ddr_data,
                                    score_seed = 260910L,
                                    subset_seed = 260911L,
                                    max_display = 500L) {
  required_cols <- c(
    "sample_id", "cancer_type", "HRDsum", "HRD_LOH",
    "LST", "TAI", "purity", "ploidy"
  )

  if (is.null(ddr_data) || !nrow(ddr_data)) {
    return(empty_hrd_scores_data())
  }

  ddr_data <- as.data.table(ddr_data)
  missing_cols <- setdiff(required_cols, names(ddr_data))
  if (length(missing_cols)) {
    stop(
      "DDR input is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  sample_data <- unique(ddr_data[, ..required_cols])
  if (!nrow(sample_data)) {
    return(empty_hrd_scores_data())
  }

  setorder(sample_data, sample_id)

  set.seed(as.integer(score_seed))
  noise_scale <- stats::sd(sample_data$HRDsum, na.rm = TRUE)
  if (!is.finite(noise_scale) || is.na(noise_scale) || noise_scale == 0) {
    noise_scale <- 1
  }
  sample_data[, epi_HRD := round(pmax(0, HRDsum + stats::rnorm(.N, sd = noise_scale / 4)), 2)]

  keep_n <- min(nrow(sample_data), as.integer(max_display))
  if (keep_n < nrow(sample_data)) {
    set.seed(as.integer(subset_seed))
    sample_data <- sample_data[sample.int(nrow(sample_data), keep_n)]
    setorder(sample_data, sample_id)
  }

  sample_data[]
}

filter_hrd_scores_data <- function(display_data,
                                   sample_ids = NULL,
                                   sample_mode = c("all", "selected")) {
  sample_mode <- match.arg(sample_mode)
  display_data <- as.data.table(copy(display_data))

  if (!nrow(display_data) || identical(sample_mode, "all")) {
    return(list(data = display_data, unavailable = character()))
  }

  if (is.null(sample_ids)) {
    sample_ids <- character()
  }

  sample_ids <- as.character(sample_ids)
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  sample_ids <- sample_ids[!duplicated(sample_ids)]

  available <- sample_ids[sample_ids %in% display_data$sample_id]
  unavailable <- sample_ids[!sample_ids %in% display_data$sample_id]

  if (!length(available)) {
    return(list(data = display_data[0], unavailable = unavailable))
  }

  list(
    data = display_data[match(available, display_data$sample_id)],
    unavailable = unavailable
  )
}

hrd_component_data <- function(display_data) {
  display_data <- as.data.table(copy(display_data))
  if (!nrow(display_data)) {
    return(data.table(sample_id = character(), component = character(), score = numeric()))
  }

  required_cols <- c("sample_id", "HRD_LOH", "LST", "TAI")
  missing_cols <- setdiff(required_cols, names(display_data))
  if (length(missing_cols)) {
    stop(
      "Display data is missing required HRD component columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  component_levels <- c("HRD_LOH", "LST", "TAI")
  melted <- melt(
    display_data[, ..required_cols],
    id.vars = "sample_id",
    variable.name = "component",
    value.name = "score"
  )
  melted[, component := factor(component, levels = component_levels)]
  melted[, sample_id := factor(sample_id, levels = unique(display_data$sample_id))]
  melted[]
}
