# linear model evaluation
library(ggplot2)
library(dplyr)
library(sf)
library(spdep)
library(lmtest)
library(car)
library(tidyr)
library(readr)
library(broom)
library(MuMIn)  


compute_rmse_mae <- function(observed, predicted, na.rm = TRUE) {
  if (length(observed) != length(predicted)) stop("⁠ observed ⁠ and ⁠ predicted ⁠ must be same length.")
  if (na.rm) {
    ok <- complete.cases(observed, predicted)
    observed <- observed[ok]
    predicted <- predicted[ok]
  }
  if (length(observed) == 0) {
    warning("No valid pairs to compute metrics; returning NA.")
    return(c(RMSE = NA_real_, MAE = NA_real_))
  }
  res <- observed - predicted
  c(RMSE = sqrt(mean(res^2)), MAE = mean(abs(res)))
}

diagnose_lm_gg <- function(model, data_sf,
                           k = 8,
                           weight_style = "W",
                           compute_local = TRUE,
                           export_summary = NULL,       
                           save_plots = FALSE,
                           plot_dir = "plots",
                           plot_prefix = "diagnostics") {
  
  if (!inherits(model, "lm")) stop("⁠ model ⁠ must be an lm object.")
  if (!inherits(data_sf, "sf")) stop("⁠ data_sf ⁠ must be an sf object.")
  if (save_plots && !dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  response <- model.response(model.frame(model))
  fitted_vals <- fitted(model)
  raw_resid <- response - fitted_vals
  std_res <- rstandard(model)
  stud_res <- rstudent(model)
  lev <- hatvalues(model)
  cookd <- cooks.distance(model)
  n <- length(fitted_vals)
  p <- length(coef(model))
  cook_cutoff <- 4 / (n - p)
  leverage_cutoff <- 2 * p / n
  
  
  shapiro_res <- if (n >= 3 && n <= 5000) {
    tryCatch(shapiro.test(std_res), error = function(e) NULL)
  } else {
    warning("Shapiro-Wilk only valid for n in [3,5000]; skipping.")
    NULL
  }
  bp_res <- tryCatch(bptest(model, studentize = TRUE), error = function(e) NULL)
  
  vif_res <- tryCatch({
    v <- car::vif(model)
    if (is.matrix(v)) v <- diag(v)
    v
  }, error = function(e) NULL)
  
 
  gl <- tryCatch(glance(model), error = function(e) NULL)
  R2 <- if (!is.null(gl) && !is.null(gl$r.squared)) gl$r.squared else summary(model)$r.squared
  adjR2 <- if (!is.null(gl) && !is.null(gl$adj.r.squared)) gl$adj.r.squared else summary(model)$adj.r.squared
  aic_val <- AIC(model)
  aicc_val <- tryCatch(AICc(model), error = function(e) {
    warning("AICc computation failed: ", e$message)
    NA_real_
  })
  
 
  coords_orig <- sf::st_coordinates(sf::st_centroid(data_sf$geometry))
  dupe_idx <- duplicated(coords_orig) | duplicated(coords_orig, fromLast = TRUE)
  coords_for_knn <- coords_orig
  if (any(dupe_idx)) {
    warning("Identical spatial coordinates detected; applying tiny jitter for neighbor search.")
    rng <- apply(coords_orig, 2, function(x) diff(range(x)))
    jitter_sd <- 1e-8 * rng
    set.seed(42)
    noise <- matrix(rnorm(nrow(coords_orig)*2, mean = 0, sd = jitter_sd), ncol = 2)
    coords_for_knn <- coords_orig + noise
  }
  knear <- spdep::knearneigh(coords_for_knn, k = k)
  nb <- spdep::knn2nb(knear)
  listw <- spdep::nb2listw(nb, style = weight_style, zero.policy = TRUE)
  
  moran_global <- tryCatch(spdep::moran.test(raw_resid, listw, zero.policy = TRUE), 
                           error = function(e) NULL)
  local_moran <- NULL
  if (compute_local && !is.null(moran_global)) {
    local_moran <- tryCatch(spdep::localmoran(raw_resid, listw, zero.policy = TRUE),
                            error = function(e) NULL)
  }
  
  
  obs_names <- rownames(model.frame(model))
  if (is.null(obs_names)) obs_names <- as.character(seq_len(n))
  diag_tbl <- tibble(
    obs = obs_names,
    fitted = as.numeric(fitted_vals),
    residual = as.numeric(raw_resid),
    std_resid = as.numeric(std_res),
    stud_resid = as.numeric(stud_res),
    leverage = as.numeric(lev),
    cooks_d = as.numeric(cookd),
    high_influence = cooks_d > cook_cutoff,
    high_leverage = leverage > leverage_cutoff
  )
  if (!is.null(local_moran)) {
    local_df <- as_tibble(local_moran)
    colnames(local_df) <- paste0("localMI_", colnames(local_df))
    diag_tbl <- bind_cols(diag_tbl, local_df)
  }
  

  error_metrics <- compute_rmse_mae(response, fitted_vals)
  

  summary_tbl <- tibble(
    R2 = as.numeric(R2),
    adjusted_R2 = as.numeric(adjR2),
    AIC = aic_val,
    AICc = aicc_val,
    RMSE = as.numeric(error_metrics["RMSE"]),
    MAE = as.numeric(error_metrics["MAE"]),
    Shapiro_W = if (!is.null(shapiro_res)) unname(shapiro_res$statistic) else NA_real_,
    Shapiro_p = if (!is.null(shapiro_res)) shapiro_res$p.value else NA_real_,
    BP_statistic = if (!is.null(bp_res)) as.numeric(bp_res$statistic) else NA_real_,
    BP_df = if (!is.null(bp_res) && !is.null(bp_res$parameter)) as.numeric(bp_res$parameter) else NA_real_,
    BP_p = if (!is.null(bp_res)) bp_res$p.value else NA_real_,
    Moran_I = if (!is.null(moran_global)) as.numeric(moran_global$estimate["Moran I statistic"]) else NA_real_,
    Moran_expectation = if (!is.null(moran_global)) as.numeric(moran_global$estimate["Expectation"]) else NA_real_,
    Moran_variance = if (!is.null(moran_global)) as.numeric(moran_global$estimate["Variance"]) else NA_real_,
    Moran_p = if (!is.null(moran_global)) as.numeric(moran_global$p.value) else NA_real_
  )
  if (!is.null(vif_res)) {
    vif_vec <- as.numeric(vif_res)
    names(vif_vec) <- paste0("VIF_", names(vif_res))
    vif_tbl <- as_tibble(t(vif_vec))
    summary_tbl <- bind_cols(summary_tbl, vif_tbl)
  }
  
  
  if (!is.null(export_summary)) {
    tryCatch({
      write_csv(summary_tbl, export_summary)
    }, error = function(e) warning("Failed to write summary table: ", e$message))
  }
  
  
  p_res_fit <- ggplot(diag_tbl, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
    theme_minimal() +
    aes(color = high_influence) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "orange"), guide = "none") +
    geom_text(data = filter(diag_tbl, high_influence),
              aes(label = obs), hjust = 1.1, vjust = 0, size = 2, check_overlap = TRUE)
  
  qqdat <- qqnorm(stud_res, plot.it = FALSE)
  qq_df <- tibble(theoretical = qqdat$x, sample = qqdat$y)
  p_qq <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Normal Q-Q (Studentized Residuals)",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  p_lev_std <- ggplot(diag_tbl, aes(x = leverage, y = std_resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed") +
    geom_vline(xintercept = leverage_cutoff, linetype = "dashed", color = "purple") +
    labs(title = "Leverage vs Standardized Residuals",
         x = "Leverage", y = "Standardized Residuals") +
    theme_minimal() +
    geom_text(data = filter(diag_tbl, high_leverage),
              aes(label = obs), hjust = 1.1, vjust = 0, size = 2, check_overlap = TRUE)
  
  diag_tbl <- diag_tbl %>% mutate(index = row_number())
  p_cook <- ggplot(diag_tbl, aes(x = index, y = cooks_d)) +
    geom_segment(aes(x = index, xend = index, y = 0, yend = cooks_d)) +
    geom_point(data = filter(diag_tbl, cooks_d > cook_cutoff), size = 2, color = "red") +
    geom_hline(yintercept = cook_cutoff, linetype = "dashed", color = "red") +
    labs(title = "Cook's Distance", x = "Observation (index)", y = "Cook's D") +
    theme_minimal() +
    geom_text(data = filter(diag_tbl, cooks_d > cook_cutoff),
              aes(label = obs), hjust = -0.1, vjust = 0, size = 2, check_overlap = TRUE)
  
  if (save_plots) {
    ggsave(file.path(plot_dir, paste0(plot_prefix, "_residuals_vs_fitted.png")),
           plot = p_res_fit, width = 6, height = 4, dpi = 300)
    ggsave(file.path(plot_dir, paste0(plot_prefix, "_qq.png")),
           plot = p_qq, width = 6, height = 4, dpi = 300)
    ggsave(file.path(plot_dir, paste0(plot_prefix, "_leverage_vs_std.png")),
           plot = p_lev_std, width = 6, height = 4, dpi = 300)
    ggsave(file.path(plot_dir, paste0(plot_prefix, "_cooks_distance.png")),
           plot = p_cook, width = 6, height = 4, dpi = 300)
  }
  
  
  cat("=== Shapiro-Wilk test (std residuals) ===\n"); print(shapiro_res)
  cat("\n=== Breusch-Pagan test (studentized) ===\n"); print(bp_res)
  cat("\n=== VIF ===\n"); print(vif_res)
  cat("\n=== R² / adjusted R² ===\n"); cat("R²:", R2, " adjusted R²:", adjR2, "\n")
  cat("\n=== AIC / AICc ===\n"); cat("AIC:", aic_val, " AICc:", aicc_val, "\n")
  cat("\n=== RMSE / MAE ===\n"); print(error_metrics)
  cat("\n=== Moran's I on residuals ===\n"); print(moran_global)
  if (!is.null(local_moran)) {
    cat("\n=== First few rows of Local Moran's I ===\n"); print(head(local_moran))
  }
  cat("\n=== High influence obs (Cook's D > ", round(cook_cutoff, 4), ") ===\n", sep = "")
  print(filter(diag_tbl, high_influence))
  cat("\n=== High leverage obs (leverage > ", round(leverage_cutoff, 4), ") ===\n", sep = "")
  print(filter(diag_tbl, high_leverage))
  
  list(
    shapiro_wilk = shapiro_res,
    breusch_pagan = bp_res,
    vif = vif_res,
    R2 = R2,
    adjusted_R2 = adjR2,
    AIC = aic_val,
    AICc = aicc_val,
    error_metrics = error_metrics,
    moran_global = moran_global,
    local_moran = local_moran,
    summary_table = summary_tbl,
    plots = list(
      residuals_vs_fitted = p_res_fit,
      qq = p_qq,
      leverage_vs_std = p_lev_std,
      cooks_distance = p_cook
    ),
    listw = listw,
    thresholds = list(cook_cutoff = cook_cutoff, leverage_cutoff = leverage_cutoff)
  )
}

#gwr evaluation
library(sf)
library(spdep)
library(GWmodel)

compute_rmse_mae <- function(observed, predicted, na.rm = TRUE) {
  if (length(observed) != length(predicted)) {
    stop("⁠ observed ⁠ and ⁠ predicted ⁠ must have the same length.")
  }
  if (na.rm) {
    ok <- complete.cases(observed, predicted)
    observed <- observed[ok]
    predicted <- predicted[ok]
  }
  if (length(observed) == 0) {
    warning("No valid pairs to compute metrics; returning NA.")
    return(c(RMSE = NA_real_, MAE = NA_real_))
  }
  res <- observed - predicted
  c(RMSE = sqrt(mean(res^2)), MAE = mean(abs(res)))
}

get_gwr_fitted <- function(gwr_res, formula, data) {
  sdf <- as.data.frame(gwr_res$SDF)
  if ("pred" %in% colnames(sdf)) {
    return(as.numeric(sdf$pred))
  }
  X <- model.matrix(formula, data = data)
  mm_cols <- colnames(X)
  sdf_names <- colnames(sdf)
  normalize <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
  norm_sdf <- vapply(sdf_names, normalize, "")
  if ("(Intercept)" %in% sdf_names) {
    intercept_sdf <- "(Intercept)"
  } else if ("X.Intercept." %in% sdf_names) {
    intercept_sdf <- "X.Intercept."
  } else {
    intercept_candidates <- grep("intercept", sdf_names, ignore.case = TRUE, value = TRUE)
    intercept_candidates <- intercept_candidates[!grepl("se($|)|se_EDF", intercept_candidates, ignore.case = TRUE)]
    if (length(intercept_candidates) >= 1) {
      intercept_sdf <- intercept_candidates[1]
    } else {
      stop("No intercept-like column found in GWR output to reconstruct fitted values.")
    }
  }
  matched <- character(length(mm_cols))
  for (i in seq_along(mm_cols)) {
    col <- mm_cols[i]
    if (col == "(Intercept)") {
      matched[i] <- intercept_sdf
      next
    }
    if (col %in% sdf_names) {
      matched[i] <- col
      next
    }
    norm_col <- normalize(col)
    exact_candidates <- sdf_names[norm_sdf == norm_col]
    if (length(exact_candidates) == 1) {
      matched[i] <- exact_candidates
      next
    }
    partials <- sdf_names[grepl(norm_col, norm_sdf)]
    if (length(partials) == 1) {
      matched[i] <- partials
      next
    }
    stop(sprintf("Could not match model matrix column '%s' to any local coefficient in GWR output. Available: %s",
                 col, paste(sdf_names, collapse = ", ")))
  }
  local_coefs <- as.matrix(sdf[, matched, drop = FALSE])
  colnames(local_coefs) <- mm_cols
  rowSums(local_coefs * X)
}

evaluate_gwr_model_csv <- function(gwr_res,
                                   global_lm,
                                   data_sf,
                                   formula,
                                   k = 8,
                                   weight_style = "W",
                                   zero_policy = TRUE,
                                   diagnostics_csv = NULL,
                                   summary_csv = NULL,
                                   montecarlo_nsim = 999) {
  if (!inherits(data_sf, "sf")) stop("⁠ data_sf ⁠ must be an sf object.")
  if (!inherits(global_lm, "lm")) stop("⁠ global_lm ⁠ must be an lm object.")
  
  mf <- model.frame(formula, data = data_sf)
  omitted <- attr(mf, "na.action")
  if (is.null(omitted)) {
    used_idx <- seq_len(nrow(data_sf))
  } else {
    used_idx <- setdiff(seq_len(nrow(data_sf)), omitted)
  }
  data_used_sf <- data_sf[used_idx, , drop = FALSE]
  response <- mf[[1]]
  n <- length(response)
  
  coords <- sf::st_coordinates(st_centroid(data_used_sf$geometry))
  dupe <- duplicated(coords) | duplicated(coords, fromLast = TRUE)
  coords_for_nb <- coords
  if (any(dupe)) {
    warning("Identical coordinates detected; jittering duplicates for neighbor search.")
    set.seed(123)
    span <- apply(coords, 2, function(x) diff(range(x)))
    jitter_sd <- 1e-8 * span
    coords_for_nb[dupe, ] <- coords[dupe, ] + matrix(rnorm(sum(dupe) * 2, mean = 0, sd = jitter_sd), ncol = 2)
  }
  knear <- spdep::knearneigh(coords_for_nb, k = k)
  nb <- spdep::knn2nb(knear)
  listw <- spdep::nb2listw(nb, style = weight_style, zero.policy = zero_policy)
  
  gwr_fitted <- get_gwr_fitted(gwr_res, formula, data_used_sf)
  sdf <- as.data.frame(gwr_res$SDF)
  reported_resid <- if ("gwr.e" %in% colnames(sdf)) as.numeric(sdf$gwr.e) else NULL
  
  if (!is.null(reported_resid)) {
    residuals_used <- reported_resid
    reconstructed_resid <- response - gwr_fitted
    diff_max <- max(abs(reconstructed_resid - reported_resid), na.rm = TRUE)
    if (diff_max > .Machine$double.eps^0.5 * 100) {
      warning("Reconstructed residuals differ from reported 'gwr.e' (max abs diff = ", signif(diff_max, 4), ").")
    }
  } else {
    residuals_used <- response - gwr_fitted
    reconstructed_resid <- residuals_used
    diff_max <- NA_real_
  }
  
  metrics <- compute_rmse_mae(response, gwr_fitted)
  
  
  moran_i <- tryCatch({
    spdep::moran.test(residuals_used, listw, zero.policy = zero_policy)
  }, error = function(e) {
    warning("Moran's I test failed: ", e$message)
    NULL
  })
  
  TSS <- sum((response - mean(response))^2, na.rm = TRUE)
  RSS <- sum((residuals_used)^2, na.rm = TRUE)
  R2_gwr <- 1 - RSS / TSS
  enp <- NULL
  if (!is.null(gwr_res$GW.diagnostics) && !is.null(gwr_res$GW.diagnostics$ENP)) {
    enp <- as.numeric(gwr_res$GW.diagnostics$ENP)
  }
  if (!is.null(enp) && !is.na(enp) && (n - enp) > 0) {
    adjR2_gwr <- 1 - (1 - R2_gwr) * (n - 1) / (n - enp)
  } else {
    p_global <- length(coef(global_lm))
    adjR2_gwr <- 1 - (1 - R2_gwr) * (n - 1) / (n - p_global)
  }
  AICc_gwr <- if (!is.null(gwr_res$GW.diagnostics) && !is.null(gwr_res$GW.diagnostics$AICc)) {
    as.numeric(gwr_res$GW.diagnostics$AICc)
  } else {
    NA_real_
  }
  
  bw_used <- bw_cv
  adaptive_flag <- if (!is.null(gwr_res$adaptive)) as.logical(gwr_res$adaptive) else NA
  
  gwr_mc <- tryCatch({
    gwr.montecarlo(formula, data = data_used_sf, bw = bw_used,
                   adaptive = ifelse(is.na(adaptive_flag), FALSE, adaptive_flag),
                   nsim = montecarlo_nsim)
  }, error = function(e) {
    warning("Monte Carlo non-stationarity test failed: ", e$message)
    NULL
  })
 
mc_pvals <- NULL
if (!is.null(gwr_mc)) {
  if (is.matrix(gwr_mc)) {
    if ("p-value" %in% colnames(gwr_mc)) {
      mc_pvals <- as.numeric(gwr_mc[, "p-value"])
      names(mc_pvals) <- rownames(gwr_mc)
    } else if (ncol(gwr_mc) == 1) {
      mc_pvals <- as.numeric(gwr_mc[, 1])
      names(mc_pvals) <- rownames(gwr_mc)
    } else {
      warning("Unexpected Monte Carlo matrix format; filling mc_pvals with NA.")
      mc_pvals <- setNames(rep(NA_real_, nrow(gwr_mc)), rownames(gwr_mc))
    }
  } else if (is.list(gwr_mc) && !is.null(gwr_mc$results)) {
    mc_df <- as.data.frame(gwr_mc$results)
    if ("Prob" %in% colnames(mc_df)) {
      mc_pvals <- mc_df$Prob
      names(mc_pvals) <- rownames(mc_df)
    } else {
      pcol <- grep("p\\.value|pvalue", tolower(colnames(mc_df)), value = TRUE)
      if (length(pcol) >= 1) {
        mc_pvals <- mc_df[[pcol[1]]]
        names(mc_pvals) <- rownames(mc_df)
      } else if (all(c("Observed", "Mean", "SD") %in% colnames(mc_df))) {
        z <- (mc_df$Observed - mc_df$Mean) / mc_df$SD
        mc_pvals <- 2 * (1 - pnorm(abs(z)))
        names(mc_pvals) <- rownames(mc_df)
      } else {
        warning("Unrecognized structure in Monte Carlo results; filling with NA.")
        mc_pvals <- setNames(rep(NA_real_, nrow(mc_df)), rownames(mc_df))
      }
    }
  } else {
    warning("Could not interpret monte carlo output format; mc_pvals set to NA.")
    
  }
}
  
 
  global_coefs <- coef(global_lm)
  predictor_names <- names(global_coefs)
  predictor_names <- predictor_names[!predictor_names %in% c("(Intercept)", "Intercept")]
  local_summary_list <- lapply(predictor_names, function(pred) {
    if (!pred %in% colnames(sdf)) {
      warning("Predictor ", pred, " not found in GWR output; skipping.")
      return(NULL)
    }
    local_coef <- as.numeric(sdf[[pred]])
    se_col <- intersect(c(paste0(pred, "_SE"), paste0(pred, "_se")), colnames(sdf))
    prop_significant <- NA_real_
    if (length(se_col) == 1) {
      local_se <- as.numeric(sdf[[se_col]])
      t_vals <- local_coef / local_se
      prop_significant <- mean(abs(t_vals) > 1.96, na.rm = TRUE)
    }
    prop_opposite_sign <- mean(sign(local_coef) != sign(global_coefs[pred]), na.rm = TRUE)
    mc_p <- if (!is.null(mc_pvals) && pred %in% names(mc_pvals)) mc_pvals[pred] else NA_real_
    data.frame(
      predictor = pred,
      global_coef = as.numeric(global_coefs[pred]),
      mean_local = mean(local_coef, na.rm = TRUE),
      sd_local = sd(local_coef, na.rm = TRUE),
      min_local = min(local_coef, na.rm = TRUE),
      max_local = max(local_coef, na.rm = TRUE),
      prop_significant = prop_significant,
      prop_opposite_sign = prop_opposite_sign,
      mc_nonstationarity_p = mc_p,
      stringsAsFactors = FALSE
    )
  })
  local_coeff_summary <- do.call(rbind, Filter(Negate(is.null), local_summary_list))
  
  
  diag_tbl <- as.data.frame(sdf)
  diag_tbl$observed <- response
  diag_tbl$gwr_fitted <- gwr_fitted
  diag_tbl$reconstructed_resid <- response - gwr_fitted
  if (!is.null(reported_resid)) diag_tbl$reported_resid <- reported_resid
  diag_tbl$residual_used <- residuals_used
  diag_tbl$abs_error <- abs(response - gwr_fitted)
  diag_tbl$sq_error <- (response - gwr_fitted)^2
  if (!is.null(moran_i)) {
    diag_tbl$global_moran_I <- as.numeric(moran_i$estimate["Moran I statistic"])
  }
  cent <- coords
  colnames(cent) <- c("centroid_x", "centroid_y")
  diag_tbl <- cbind(diag_tbl, cent)
  

  summary_tbl <- tibble::tibble(
    n_obs = length(response),
    bandwidth = bw_cv,
    adaptive = adaptive_flag,
    GWR_R2 = R2_gwr,
    GWR_adjR2 = adjR2_gwr,
    GWR_AICc = AICc_gwr,
    RMSE = as.numeric(metrics["RMSE"]),
    MAE = as.numeric(metrics["MAE"]),
    Moran_I = if (!is.null(moran_i)) as.numeric(moran_i$estimate["Moran I statistic"]) else NA_real_,
    Moran_expectation = if (!is.null(moran_i)) as.numeric(moran_i$estimate["Expectation"]) else NA_real_,
    Moran_variance = if (!is.null(moran_i)) as.numeric(moran_i$estimate["Variance"]) else NA_real_,
    Moran_p = if (!is.null(moran_i)) as.numeric(moran_i$p.value) else NA_real_,
    resid_diff_max = diff_max,
    mean_localR2 = if ("localR2" %in% colnames(sdf)) mean(sdf$localR2, na.rm = TRUE) else NA_real_
  )
    
if (!is.null(mc_pvals)) {
  for (nm in names(mc_pvals)) {
    clean_name <- gsub("[^A-Za-z0-9_]", "", nm)  # safe column name
    colname <- paste0("MC_p_", clean_name)
    summary_tbl[[colname]] <- mc_pvals[nm]
  }
  }
  
  if (!is.null(diagnostics_csv)) {
    tryCatch({
      utils::write.csv(diag_tbl, diagnostics_csv, row.names = FALSE)
    }, error = function(e) warning("Failed to write diagnostics CSV: ", e$message))
  }
  if (!is.null(summary_csv)) {
    tryCatch({
      utils::write.csv(summary_tbl, summary_csv, row.names = FALSE)
    }, error = function(e) warning("Failed to write summary CSV: ", e$message))
  }
  
  
  list(
    gwr_result = gwr_res,
    fitted = gwr_fitted,
    residuals = residuals_used,
    metrics = metrics,
    moran = moran_i,
    diagnostics_table = diag_tbl,
    summary_table = summary_tbl,
    local_coefficients_summary = local_coeff_summary,
    monte_carlo = gwr_mc,
    listw = listw
  )
}