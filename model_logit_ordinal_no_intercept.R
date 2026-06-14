label_levels <- c("DS", "D0", "D1", "D2", "D3", "D4")
selected_covariates <- c("P_mm", "slope_deg", "U10_ms", "SP_Pa", "NDVI", "LST_day_C")

get_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
    character(1)
  )
  frame_files <- stats::na.omit(frame_files)
  if (length(frame_files) > 0L) {
    return(dirname(normalizePath(frame_files[[length(frame_files)]])))
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_match <- grep(file_arg, cmd_args)
  if (length(file_match) > 0L) {
    return(dirname(normalizePath(sub(file_arg, "", cmd_args[file_match[1L]]))))
  }

  normalizePath(getwd())
}

find_project_final_root <- function(start = get_script_dir()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    marker <- file.path(path, "data", "df_sequia_2007_2025.parquet")
    if (file.exists(marker)) {
      return(path)
    }

    nested_marker <- file.path(path, "Proyecto_final", "data", "df_sequia_2007_2025.parquet")
    if (file.exists(nested_marker)) {
      return(file.path(path, "Proyecto_final"))
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("No se encontro la raiz de Proyecto_final.")
    }

    path <- parent
  }
}

script_dir <- get_script_dir()
project_final_root <- find_project_final_root(script_dir)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Faltan paquetes de R: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

default_fit_options <- function() {
  list(
    data_path = file.path(project_final_root, "data", "df_sequia_2007_2025.parquet"),
    label_col = "clase_sequia",
    label_levels = label_levels,
    covariates = selected_covariates,
    stan_file = file.path(script_dir, "model_logit_ordinal_no_intercept.stan"),
    output_dir = file.path(project_final_root, "modelo_ordinal_bayesiano", "outputs"),
    train_year_max = 2021L,
    test_year_min = 2022L,
    test_year_max = 2023L,
    scale_covariates = TRUE,
    missing_value_threshold = -999,
    seed = 20260503L,
    chains = 10L,
    parallel_chains = 10L,
    iter_warmup = 1000L,
    iter_sampling = 2000L,
    refresh = 10L,
    adapt_delta = 0.90,
    max_treedepth = 12L,
    step_size = 0.02,
    use_inits = TRUE,
    init_seed = 20260504L,
    save_warmup = FALSE
  )
}

read_model_data <- function(data_path) {
  extension <- tolower(tools::file_ext(data_path))

  if (identical(extension, "parquet")) {
    require_packages("arrow")
    return(as.data.frame(arrow::read_parquet(data_path)))
  }

  if (extension %in% c("csv", "txt")) {
    return(utils::read.csv(data_path, stringsAsFactors = FALSE))
  }

  stop("Formato de datos no soportado: ", extension)
}

make_time_key <- function(df) {
  if ("date_ym" %in% names(df)) {
    return(as.character(df$date_ym))
  }

  if (all(c("year", "month") %in% names(df))) {
    return(sprintf("%04d-%02d", as.integer(df$year), as.integer(df$month)))
  }

  stop("El dataset debe contener date_ym o year/month.")
}

make_month_serial <- function(year, month) {
  as.integer(year) * 12L + as.integer(month)
}

make_time_delta <- function(time_serial) {
  deltas <- c(1L, diff(as.integer(time_serial)))
  if (any(deltas < 1L)) {
    stop("La secuencia temporal no esta ordenada correctamente.")
  }
  as.integer(deltas)
}

prepare_ordinal_dynamic_data <- function(
  data_path,
  label_col = "clase_sequia",
  label_levels = label_levels,
  covariates = selected_covariates,
  train_year_max = 2021L,
  test_year_min = 2022L,
  test_year_max = 2023L,
  scale_covariates = TRUE,
  missing_value_threshold = -999
) {
  df <- read_model_data(data_path)

  required_vars <- c("cell_id", "year", "month", "date_ym", label_col, covariates)
  missing_vars <- setdiff(required_vars, names(df))
  if (length(missing_vars) > 0L) {
    stop("Variables no encontradas en el dataset: ", paste(missing_vars, collapse = ", "))
  }

  sentinel_counts <- vapply(
    covariates,
    function(var_name) {
      sum(df[[var_name]] <= missing_value_threshold, na.rm = TRUE)
    },
    integer(1)
  )

  for (var_name in covariates) {
    sentinel_idx <- !is.na(df[[var_name]]) & df[[var_name]] <= missing_value_threshold
    df[[var_name]][sentinel_idx] <- NA_real_
  }

  n_before_complete_cases <- nrow(df)
  keep_idx <- stats::complete.cases(df[, c("cell_id", "year", "month", label_col, covariates)])
  df <- df[keep_idx, ]
  n_removed_complete_cases <- n_before_complete_cases - nrow(df)

  df$time_key <- make_time_key(df)
  df$time_serial <- make_month_serial(df$year, df$month)
  df <- df[order(df$cell_id, df$time_serial), ]

  y <- match(as.character(df[[label_col]]), label_levels)
  if (anyNA(y)) {
    bad_labels <- sort(unique(as.character(df[[label_col]])[is.na(y)]))
    stop("Etiquetas fuera de label_levels: ", paste(bad_labels, collapse = ", "))
  }
  df$y <- as.integer(y)
  df[[label_col]] <- factor(df[[label_col]], levels = label_levels, ordered = TRUE)

  train <- df[df$year <= train_year_max, ]
  test <- df[df$year >= test_year_min & df$year <= test_year_max, ]

  if (nrow(train) == 0L || nrow(test) == 0L) {
    stop("La particion train/test no produjo observaciones suficientes.")
  }

  train_cells <- sort(unique(train$cell_id))
  unseen_test_cells <- setdiff(unique(test$cell_id), train_cells)
  if (length(unseen_test_cells) > 0L) {
    stop("Hay celdas en prueba no vistas en entrenamiento: ", paste(head(unseen_test_cells), collapse = ", "))
  }

  cells <- train_cells
  train_times <- sort(unique(train$time_serial))
  test_times <- sort(unique(test$time_serial))
  train_time_keys <- train$time_key[match(train_times, train$time_serial)]
  test_time_keys <- test$time_key[match(test_times, test$time_serial)]

  if (min(test_times) <= max(train_times)) {
    stop("Los tiempos de prueba deben ser posteriores al entrenamiento.")
  }

  x_train_raw <- as.matrix(train[, covariates])
  x_test_raw <- as.matrix(test[, covariates])
  x_center <- rep(0, ncol(x_train_raw))
  x_scale <- rep(1, ncol(x_train_raw))

  if (isTRUE(scale_covariates)) {
    x_center <- colMeans(x_train_raw)
    x_scale <- apply(x_train_raw, 2, stats::sd)
    if (any(!is.finite(x_scale)) || any(x_scale == 0)) {
      stop("No fue posible escalar las covariables.")
    }
    x_train <- scale(x_train_raw, center = x_center, scale = x_scale)
    x_test <- scale(x_test_raw, center = x_center, scale = x_scale)
  } else {
    x_train <- x_train_raw
    x_test <- x_test_raw
  }

  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  colnames(x_train) <- covariates
  colnames(x_test) <- covariates

  cell_train <- match(train$cell_id, cells)
  cell_test <- match(test$cell_id, cells)
  time_train <- match(train$time_serial, train_times)
  time_test <- match(test$time_serial, test_times)
  month_train <- as.integer(train$month)
  month_test <- as.integer(test$month)

  train_delta <- make_time_delta(train_times)
  test_delta <- as.integer(c(test_times[1] - max(train_times), diff(test_times)))
  if (any(test_delta < 1L)) {
    stop("La secuencia de prueba no esta ordenada correctamente.")
  }

  stan_data <- list(
    N_train = nrow(train),
    N_test = nrow(test),
    I = length(cells),
    T_train = length(train_times),
    T_test = length(test_times),
    K = length(covariates),
    J = length(label_levels),
    X_train = x_train,
    y_train = as.integer(train$y),
    cell_train = as.integer(cell_train),
    time_train = as.integer(time_train),
    month_train = as.integer(month_train),
    train_delta = as.integer(train_delta),
    X_test = x_test,
    y_test = as.integer(test$y),
    cell_test = as.integer(cell_test),
    time_test = as.integer(time_test),
    month_test = as.integer(month_test),
    test_delta = as.integer(test_delta)
  )

  list(
    stan_data = stan_data,
    train = train,
    test = test,
    label_col = label_col,
    label_levels = label_levels,
    covariates = covariates,
    cells = cells,
    train_times = train_times,
    test_times = test_times,
    train_time_keys = train_time_keys,
    test_time_keys = test_time_keys,
    x_center = x_center,
    x_scale = x_scale,
    scale_covariates = scale_covariates,
    missing_value_threshold = missing_value_threshold,
    sentinel_counts = sentinel_counts,
    n_removed_complete_cases = n_removed_complete_cases,
    train_year_max = train_year_max,
    test_year_min = test_year_min,
    test_year_max = test_year_max,
    train_delta = train_delta,
    test_delta = test_delta
  )
}

make_inits <- function(stan_data, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  category_counts <- tabulate(stan_data$y_train, nbins = stan_data$J)
  category_probs <- category_counts / sum(category_counts)
  cumulative_probs <- cumsum(category_probs)[seq_len(stan_data$J - 1L)]
  cumulative_probs <- pmin(pmax(cumulative_probs, 0.01), 0.99)
  cutpoints <- stats::qlogis(cumulative_probs)

  for (j in seq_along(cutpoints)[-1L]) {
    if (cutpoints[j] <= cutpoints[j - 1L]) {
      cutpoints[j] <- cutpoints[j - 1L] + 0.05
    }
  }

  function(chain_id = 1L) {
    jitter <- 0.01 * (chain_id - 1L)

    list(
      beta = rep(0, stan_data$K),
      cutpoints = cutpoints + jitter,
      rho = 0.5,
      sigma_state = 0.1,
      state_raw = matrix(0, nrow = stan_data$I, ncol = stan_data$T_train),
      sigma_cell = 0.1,
      cell_raw = rep(0, stan_data$I),
      sigma_month = 0.1,
      month_raw = rep(0, 12L)
    )
  }
}

extract_prediction_probabilities <- function(draws_matrix, levels = label_levels) {
  draws_matrix <- as.matrix(draws_matrix)
  n_obs <- ncol(draws_matrix)
  prob <- matrix(0, nrow = n_obs, ncol = length(levels))

  for (j in seq_along(levels)) {
    prob[, j] <- colMeans(draws_matrix == j)
  }

  prob <- as.data.frame(prob)
  names(prob) <- levels
  prob
}

log_mean_exp <- function(x) {
  x_max <- max(x)
  x_max + log(mean(exp(x - x_max)))
}

extract_pointwise_log_score <- function(log_lik_matrix) {
  apply(as.matrix(log_lik_matrix), 2, log_mean_exp)
}

classification_metrics <- function(actual, prob, levels = label_levels, pointwise_log_score = NULL) {
  actual <- factor(actual, levels = levels, ordered = TRUE)
  prob <- as.data.frame(prob[, levels, drop = FALSE])
  pred <- factor(levels[max.col(as.matrix(prob), ties.method = "first")], levels = levels, ordered = TRUE)

  actual_num <- as.integer(actual)
  pred_num <- as.integer(pred)
  eps <- 1e-12
  idx <- cbind(seq_along(actual_num), actual_num)
  p_actual <- pmax(as.matrix(prob)[idx], eps)

  confusion <- table(
    actual = factor(actual, levels = levels),
    predicted = factor(pred, levels = levels)
  )

  per_class <- do.call(
    rbind,
    lapply(levels, function(level) {
      tp <- confusion[level, level]
      fp <- sum(confusion[, level]) - tp
      fn <- sum(confusion[level, ]) - tp
      precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
      recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
      f1 <- if (!is.finite(precision + recall) || (precision + recall) == 0) {
        NA_real_
      } else {
        2 * precision * recall / (precision + recall)
      }
      data.frame(class = level, support = sum(confusion[level, ]), precision = precision, recall = recall, f1 = f1)
    })
  )

  y_onehot <- stats::model.matrix(~ actual - 1)
  colnames(y_onehot) <- sub("^actual", "", colnames(y_onehot))
  y_onehot <- y_onehot[, levels, drop = FALSE]
  cumulative_actual <- t(apply(y_onehot, 1, cumsum))
  cumulative_prob <- t(apply(as.matrix(prob), 1, cumsum))
  severe_actual <- actual_num >= match("D2", levels)
  severe_pred_prob <- rowSums(prob[, c("D2", "D3", "D4"), drop = FALSE])
  severe_pred <- severe_pred_prob >= 0.5

  log_loss <- if (is.null(pointwise_log_score)) {
    -mean(log(p_actual))
  } else {
    -mean(pointwise_log_score)
  }

  metrics <- data.frame(
    metric = c(
      "n",
      "accuracy",
      "balanced_accuracy",
      "macro_f1",
      "mean_abs_ordinal_error",
      "rmse_ordinal_error",
      "adjacent_accuracy",
      "log_loss",
      "brier_score",
      "ranked_probability_score",
      "severe_brier_D2_D4",
      "severe_accuracy_D2_D4"
    ),
    value = c(
      length(actual),
      mean(pred == actual),
      mean(per_class$recall, na.rm = TRUE),
      mean(per_class$f1, na.rm = TRUE),
      mean(abs(pred_num - actual_num)),
      sqrt(mean((pred_num - actual_num)^2)),
      mean(abs(pred_num - actual_num) <= 1),
      log_loss,
      mean(rowSums((as.matrix(prob) - y_onehot)^2)),
      mean(rowSums((cumulative_prob[, -length(levels), drop = FALSE] -
        cumulative_actual[, -length(levels), drop = FALSE])^2)),
      mean((severe_pred_prob - as.numeric(severe_actual))^2),
      mean(severe_pred == severe_actual)
    )
  )

  list(
    metrics = metrics,
    per_class = per_class,
    confusion = as.data.frame(confusion),
    predicted = pred
  )
}

make_prediction_frame <- function(data, prob, pred, levels = label_levels) {
  prob <- prob[, levels, drop = FALSE]
  names(prob) <- paste0("prob_", levels)

  keep <- intersect(c("cell_id", "x", "y", "year", "month", "date_ym", "clase_sequia", "y"), names(data))
  cbind(
    data[, keep, drop = FALSE],
    predicted_class = as.character(pred),
    prob
  )
}

plot_confusion <- function(confusion_df, title) {
  require_packages("ggplot2")
  ggplot2::ggplot(confusion_df, ggplot2::aes(x = predicted, y = actual, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), size = 3) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
    ggplot2::labs(title = title, x = "Predicha", y = "Observada", fill = "n") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_time_severity <- function(prediction_df, title) {
  require_packages("ggplot2")
  prediction_df$predicted_num <- match(prediction_df$predicted_class, label_levels)
  by_time <- stats::aggregate(cbind(y, predicted_num) ~ date_ym, data = prediction_df, FUN = mean)
  by_time$date <- as.Date(paste0(by_time$date_ym, "-01"))

  ggplot2::ggplot(by_time, ggplot2::aes(x = date)) +
    ggplot2::geom_line(ggplot2::aes(y = y, color = "observada"), linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = predicted_num, color = "predicha"), linewidth = 0.7) +
    ggplot2::scale_y_continuous(breaks = seq_along(label_levels), labels = label_levels) +
    ggplot2::labs(title = title, x = "Fecha", y = "Severidad ordinal media", color = "") +
    ggplot2::theme_minimal(base_size = 11)
}

write_text <- function(x, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(x, con = con)
}

save_fit_outputs <- function(result, fit, run_dir) {
  require_packages(c("posterior", "readr", "ggplot2"))

  draws_y_train <- fit$draws("y_pred_train", format = "matrix")
  draws_y_test <- fit$draws("y_pred_test", format = "matrix")
  log_lik_train <- fit$draws("log_lik_train", format = "matrix")
  log_lik_test <- fit$draws("log_lik_test", format = "matrix")

  train_prob <- extract_prediction_probabilities(draws_y_train, result$label_levels)
  test_prob <- extract_prediction_probabilities(draws_y_test, result$label_levels)
  train_log_score <- extract_pointwise_log_score(log_lik_train)
  test_log_score <- extract_pointwise_log_score(log_lik_test)

  train_eval <- classification_metrics(
    actual = result$train[[result$label_col]],
    prob = train_prob,
    levels = result$label_levels,
    pointwise_log_score = train_log_score
  )
  test_eval <- classification_metrics(
    actual = result$test[[result$label_col]],
    prob = test_prob,
    levels = result$label_levels,
    pointwise_log_score = test_log_score
  )

  train_predictions <- make_prediction_frame(result$train, train_prob, train_eval$predicted, result$label_levels)
  test_predictions <- make_prediction_frame(result$test, test_prob, test_eval$predicted, result$label_levels)

  readr::write_csv(train_eval$metrics, file.path(run_dir, "train_metrics.csv"))
  readr::write_csv(test_eval$metrics, file.path(run_dir, "test_metrics.csv"))
  readr::write_csv(train_eval$per_class, file.path(run_dir, "train_per_class.csv"))
  readr::write_csv(test_eval$per_class, file.path(run_dir, "test_per_class.csv"))
  readr::write_csv(train_eval$confusion, file.path(run_dir, "train_confusion.csv"))
  readr::write_csv(test_eval$confusion, file.path(run_dir, "test_confusion.csv"))
  readr::write_csv(train_predictions, file.path(run_dir, "train_predictions.csv"))
  readr::write_csv(test_predictions, file.path(run_dir, "test_predictions.csv"))

  diagnostics <- fit$diagnostic_summary()
  capture.output(print(diagnostics), file = file.path(run_dir, "diagnostic_summary.txt"))

  summary_pars <- c("beta", "cutpoints", "rho", "sigma_state", "sigma_cell", "sigma_month")
  fit_summary <- fit$summary(variables = summary_pars)
  readr::write_csv(fit_summary, file.path(run_dir, "posterior_summary.csv"))

  ggplot2::ggsave(
    file.path(run_dir, "test_confusion.png"),
    plot = plot_confusion(test_eval$confusion, "Modelo bayesiano ordinal sin intercepto - matriz de confusion en prueba"),
    width = 7,
    height = 5,
    dpi = 160
  )
  ggplot2::ggsave(
    file.path(run_dir, "test_time_severity.png"),
    plot = plot_time_severity(test_predictions, "Modelo bayesiano ordinal sin intercepto - severidad media en prueba"),
    width = 8,
    height = 4,
    dpi = 160
  )

  list(
    train_metrics = train_eval$metrics,
    test_metrics = test_eval$metrics,
    posterior_summary = fit_summary,
    diagnostics = diagnostics
  )
}

fit_ordinal_dynamic_no_intercept <- function(
  data_path = file.path(project_final_root, "data", "df_sequia_2007_2025.parquet"),
  label_col = "clase_sequia",
  label_levels = label_levels,
  covariates = selected_covariates,
  stan_file = file.path(script_dir, "model_logit_ordinal_no_intercept.stan"),
  output_dir = file.path(project_final_root, "modelo_ordinal_bayesiano", "outputs"),
  train_year_max = 2021L,
  test_year_min = 2022L,
  test_year_max = 2023L,
  scale_covariates = TRUE,
  missing_value_threshold = -999,
  seed = 20260503L,
  chains = 4L,
  parallel_chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  refresh = 100L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  step_size = 0.02,
  use_inits = TRUE,
  init_seed = 20260504L,
  save_warmup = FALSE
) {
  require_packages(c("cmdstanr", "posterior", "readr", "ggplot2"))
  ensure_dir(output_dir)

  prepared <- prepare_ordinal_dynamic_data(
    data_path = data_path,
    label_col = label_col,
    label_levels = label_levels,
    covariates = covariates,
    train_year_max = train_year_max,
    test_year_min = test_year_min,
    test_year_max = test_year_max,
    scale_covariates = scale_covariates,
    missing_value_threshold = missing_value_threshold
  )

  model <- cmdstanr::cmdstan_model(stan_file)
  init_fun <- if (isTRUE(use_inits)) {
    make_inits(prepared$stan_data, seed = init_seed)
  } else {
    NULL
  }

  sample_args <- list(
    data = prepared$stan_data,
    seed = seed,
    init = init_fun,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    save_warmup = save_warmup
  )
  if (!is.null(step_size)) {
    sample_args$step_size <- step_size
  }

  fit <- do.call(model$sample, sample_args)

  run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(output_dir, paste0("ordinal_no_intercept_", run_stamp))
  ensure_dir(run_dir)

  fit$save_output_files(dir = run_dir)
  saveRDS(prepared, file.path(run_dir, "prepared_data.rds"))
  saveRDS(fit, file.path(run_dir, "fit.rds"))

  metadata <- list(
    model_name = "model_logit_ordinal_no_intercept",
    no_intercept = TRUE,
    data_path = normalizePath(data_path),
    stan_file = normalizePath(stan_file),
    label_col = label_col,
    label_levels = label_levels,
    covariates = covariates,
    train_period = paste(min(prepared$train$time_key), max(prepared$train$time_key), sep = " / "),
    test_period = paste(min(prepared$test$time_key), max(prepared$test$time_key), sep = " / "),
    train_n = nrow(prepared$train),
    test_n = nrow(prepared$test),
    cells = length(prepared$cells),
    train_times = length(prepared$train_times),
    test_times = length(prepared$test_times),
    test_time_gaps = prepared$test_delta[prepared$test_delta > 1L],
    scale_covariates = scale_covariates,
    x_center = prepared$x_center,
    x_scale = prepared$x_scale,
    sentinel_counts = prepared$sentinel_counts,
    n_removed_complete_cases = prepared$n_removed_complete_cases,
    seed = seed,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    step_size = step_size,
    beta_prior = "normal(0, 0.7)",
    cutpoint_prior = "normal(0, 3)",
    rho_prior = "normal(0, 0.5) truncated to [-0.99, 0.99]",
    sigma_priors = "exponential(1)"
  )
  saveRDS(metadata, file.path(run_dir, "run_metadata.rds"))

  output_summary <- save_fit_outputs(prepared, fit, run_dir)

  write_text(
    c(
      "Modelo bayesiano ordinal dinamico sin intercepto",
      "",
      "Cambio principal: se elimina alpha para evitar competencia de localizacion con los umbrales ordinales.",
      "La localizacion base queda absorbida por cutpoints.",
      "",
      paste("Entrenamiento:", metadata$train_period, "-", metadata$train_n, "observaciones"),
      paste("Prueba:", metadata$test_period, "-", metadata$test_n, "observaciones"),
      paste("Celdas:", metadata$cells),
      paste("Tiempos train/test:", metadata$train_times, "/", metadata$test_times),
      paste("Huecos temporales en prueba:", if (length(metadata$test_time_gaps) == 0L) "ninguno" else paste(metadata$test_time_gaps, collapse = ", ")),
      "",
      "Metricas de prueba:",
      paste(output_summary$test_metrics$metric, signif(output_summary$test_metrics$value, 5), sep = ": ")
    ),
    file.path(run_dir, "run_summary.txt")
  )

  list(
    fit = fit,
    data = prepared,
    run_dir = run_dir,
    output_summary = output_summary
  )
}

# Ejemplo de ejecucion interactiva:
# fit_options <- default_fit_options()
# fit_options$iter_warmup <- 500L
# fit_options$iter_sampling <- 500L
# modelo_bayesiano <- do.call(fit_ordinal_dynamic_no_intercept, fit_options)
#
# Para ejecutar desde terminal:
# RUN_ORDINAL_BAYES=1 Rscript model_logit_ordinal_no_intercept.R

if (identical(Sys.getenv("RUN_ORDINAL_BAYES"), "1")) {
  fit_options <- default_fit_options()
  modelo_bayesiano <- do.call(fit_ordinal_dynamic_no_intercept, fit_options)
  message("Ajuste completado en: ", modelo_bayesiano$run_dir)
}
