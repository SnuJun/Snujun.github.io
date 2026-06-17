################################################################################
## main_adahard_baseline_comparison.R
## Baseline: Adaptive Hard Thresholding
## Compare: AdaHard, LPD(AdaHard), Xue_soft, Rothman, COMET_soft
################################################################################

library(dplyr)
library(MASS)
library(PDSCE)

################################################################################
## 0. Basic settings
################################################################################

wd <- "~/LPD_section_4_5/LPD_revision_v2"

result_dir <- file.path(wd, "results_revision_v2")
raw_dir <- file.path(result_dir, "raw_by_setting")
log_dir <- file.path(wd, "logs_revision_v2")

setwd(wd)

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

source("data generate.R")
source("Estimators.R")
source("getLPD.R")
source("Xue method.R")
source("Rothman method.R")
source("ICF Algorithm.R")

n <- 100
p_set <- c(30, 60, 120)
model_set <- c("m1", "m2", "m3", "m4")
iter <- 100
eps <- 1e-4
tau <- 1e-5

cov_list <- list(
  m1 = cov.mtx.m1,
  m2 = cov.mtx.m2,
  m3 = cov.mtx.m3,
  m4 = cov.mtx.m4
)

model_names <- c(
  m1 = "Banded",
  m2 = "Block",
  m3 = "Toeplitz",
  m4 = "Sparse"
)

model_levels <- c("Banded", "Block", "Toeplitz", "Sparse")

norm_levels <- c(
  "Spectral",
  "Frobenius",
  "L-infinity",
  "Element-wise max"
)

method_levels <- c(
  "AdaHard",
  "LPD_spec",
  "LPD_frob",
  "LPD_linf",
  "LPD_emax",
  "Xue_soft",
  "Rothman",
  "COMET"
)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  stop("Please provide array_id, e.g. Rscript main_adahard_baseline_comparison.R 1")
}

array_id <- as.integer(args[1])

setting_grid <- expand.grid(
  model_key = model_set,
  p_dim = p_set,
  stringsAsFactors = FALSE
)

if(array_id < 1 || array_id > nrow(setting_grid)){
  stop("array_id must be between 1 and ", nrow(setting_grid))
}

model_key <- setting_grid$model_key[array_id]
p_dim <- setting_grid$p_dim[array_id]
model_name <- model_names[model_key]

cat("============================================================\n")
cat("Array ID:", array_id, "\n")
cat("Model:", model_key, "-", model_name, "\n")
cat("n:", n, "p:", p_dim, "iter:", iter, "\n")
cat("============================================================\n")

################################################################################
## 1. Utility functions
################################################################################

norm_error <- function(est, true, norm_type){
  
  if(norm_type == "Spectral"){
    return(norm(est - true, type = "2"))
  }
  
  if(norm_type == "Frobenius"){
    return(norm(est - true, type = "F") / sqrt(nrow(true)))
  }
  
  if(norm_type == "L-infinity"){
    return(norm(est - true, type = "I"))
  }
  
  if(norm_type == "Element-wise max"){
    return(max(abs(est - true)))
  }
  
  stop("Undefined norm_type.")
}

get_relative_errors <- function(est, baseline, true){
  
  res <- sapply(norm_levels, function(nn){
    
    denom <- norm_error(baseline, true, nn)
    
    if(is.na(denom) || denom == 0){
      return(NA_real_)
    }
    
    norm_error(est, true, nn) / denom
  })
  
  data.frame(
    norm = names(res),
    rel_error = as.numeric(res)
  )
}

get_fpr_tpr <- function(est, true, tau = 1e-5){
  
  if(is.null(est) || any(is.na(est)) || any(is.infinite(est))){
    return(data.frame(
      FPR = NA_real_,
      TPR = NA_real_
    ))
  }
  
  off_diag <- row(true) != col(true)
  
  true_zero <- abs(true) < tau & off_diag
  true_nonzero <- abs(true) >= tau & off_diag
  
  est_nonzero <- abs(est) >= tau & off_diag
  
  FPR <- sum(est_nonzero & true_zero) / sum(true_zero)
  TPR <- sum(est_nonzero & true_nonzero) / sum(true_nonzero)
  
  data.frame(
    FPR = FPR,
    TPR = TPR
  )
}

get_pd_info <- function(est, eps = 1e-4){
  
  if(is.null(est) || any(is.na(est)) || any(is.infinite(est))){
    return(data.frame(
      min_eig = NA_real_,
      is_pd = FALSE,
      is_eps_pd = FALSE
    ))
  }
  
  eig_val <- eigen(est, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_val)
  
  data.frame(
    min_eig = min_eig,
    is_pd = min_eig > 0,
    is_eps_pd = min_eig > eps
  )
}

################################################################################
## 2. Estimator function
################################################################################

get_estimators <- function(X, eps = 1e-4){
  
  ## Baseline: Adaptive hard thresholding estimator
  AdaHard <- tryCatch({
    ftn.cov.adap(
      datamat = X,
      thres.ftn = thres.hard,
      thresseq = seq(0, 1, length.out = 100),
      diag.thres = FALSE
    )$estimate
  }, error = function(e){
    message("AdaHard failed: ", e$message)
    return(NULL)
  })
  
  ## Soft thresholding estimator used only for Xue/Rothman/COMET support
  Soft <- tryCatch({
    soft_thresholding_CV(X)$mat
  }, error = function(e){
    message("Soft failed: ", e$message)
    return(NULL)
  })
  
  ## LPD estimators based on AdaHard
  LPD_spec <- tryCatch({
    getLPD(AdaHard, lwr_bd = eps, type_norm = "spectral")$LPD
  }, error = function(e){
    message("LPD_spec failed: ", e$message)
    return(NULL)
  })
  
  LPD_frob <- tryCatch({
    getLPD(AdaHard, lwr_bd = eps, type_norm = "Frobenius")$LPD
  }, error = function(e){
    message("LPD_frob failed: ", e$message)
    return(NULL)
  })
  
  LPD_linf <- tryCatch({
    getLPD(AdaHard, lwr_bd = eps, type_norm = "Linf")$LPD
  }, error = function(e){
    message("LPD_linf failed: ", e$message)
    return(NULL)
  })
  
  LPD_emax <- tryCatch({
    getLPD(AdaHard, lwr_bd = eps, type_norm = "elemax")$LPD
  }, error = function(e){
    message("LPD_emax failed: ", e$message)
    return(NULL)
  })
  
  ## Xue method with soft thresholding
  Xue_soft <- tryCatch({
    ftn.cov.Xue(
      datamat   = X,
      thres.ftn = thres.soft,
      thresseq  = seq(0, 1, length.out = 100),
      mineigval = eps,
      tol       = 1e-4,
      MAXITER   = 10000,
      message   = FALSE
    )$est
  }, error = function(e){
    message("Xue_soft failed: ", e$message)
    return(NULL)
  })
  
  ## Rothman method: PDSCE::pdsoft.cv()
  Rothman <- tryCatch({
    ftn.cov.Roth(
      datamat = X,
      thresseq = seq(0, 1, length.out = 100),
      tol = 1e-7,
      MAXITER = 10000
    )$est
  }, error = function(e){
    message("Rothman failed: ", e$message)
    return(NULL)
  })
  
  ## COMET with soft-thresholding support
  COMET <- tryCatch({
    
    if(is.null(AdaHard)){
      stop("AdaHard estimator is NULL, so COMET support cannot be constructed.")
    }
    
    amat <- (abs(AdaHard) > 1e-8) * 1
    diag(amat) <- 1
    
    comet_res <- ICF(
      amat = amat,
      dat = X,
      cc = eps,
      tol = 1e-2,
      MAXiter = 10000
    )
    
    comet_res$mat
    
  }, error = function(e){
    message("COMET failed: ", e$message)
    return(NULL)
  })
  
  list(
    AdaHard = AdaHard,
    LPD_spec = LPD_spec,
    LPD_frob = LPD_frob,
    LPD_linf = LPD_linf,
    LPD_emax = LPD_emax,
    Xue_soft = Xue_soft,
    Rothman = Rothman,
    COMET = COMET
  )
}

################################################################################
## 3. Simulation: one setting per array task
################################################################################

set.seed(1234 + array_id)

output_file <- file.path(
  raw_dir,
  paste0("result_", model_key, "_p", p_dim, "_n", n, "_iter", iter, ".rds")
)

workspace_file <- file.path(
  raw_dir,
  paste0("workspace_", model_key, "_p", p_dim, "_n", n, "_iter", iter, ".RData")
)

if(file.exists(output_file)){
  cat("Output already exists. Skip:\n")
  cat(output_file, "\n")
  quit(save = "no")
}

rel_list <- list()
fpr_tpr_list <- list()
pd_list <- list()

Sigma <- cov_list[[model_key]](p_dim)

eig <- eigen(Sigma, symmetric = TRUE)
sqrt_Sigma <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

for(b in 1:iter){
  
  cat("iter:", b, "/", iter, "\n")
  
  X <- matrix(rnorm(n * p_dim), nrow = n, ncol = p_dim) %*% sqrt_Sigma
  X <- scale(X, center = TRUE, scale = FALSE)
  
  est_list <- get_estimators(X, eps = eps)
  baseline <- est_list$AdaHard
  
  for(est_name in names(est_list)){
    
    est_now <- est_list[[est_name]]
    
    ## Relative error against AdaHard
    if(is.null(est_now) || is.null(baseline)){
      rel_df <- data.frame(
        norm = norm_levels,
        rel_error = NA_real_
      )
    } else{
      rel_df <- get_relative_errors(
        est = est_now,
        baseline = baseline,
        true = Sigma
      )
    }
    
    rel_df$model_key <- model_key
    rel_df$model <- model_name
    rel_df$n <- n
    rel_df$p_dim <- p_dim
    rel_df$iter <- b
    rel_df$method <- est_name
    
    rel_list[[length(rel_list) + 1]] <- rel_df
    
    ## FPR / TPR
    fpr_tpr_now <- get_fpr_tpr(est_now, Sigma, tau = tau)
    fpr_tpr_now$model_key <- model_key
    fpr_tpr_now$model <- model_name
    fpr_tpr_now$n <- n
    fpr_tpr_now$p_dim <- p_dim
    fpr_tpr_now$iter <- b
    fpr_tpr_now$method <- est_name
    
    fpr_tpr_list[[length(fpr_tpr_list) + 1]] <- fpr_tpr_now
    
    ## PD information
    pd_now <- get_pd_info(est_now, eps = eps)
    pd_now$model_key <- model_key
    pd_now$model <- model_name
    pd_now$n <- n
    pd_now$p_dim <- p_dim
    pd_now$iter <- b
    pd_now$method <- est_name
    
    pd_list[[length(pd_list) + 1]] <- pd_now
  }
  
  ## Temporary save every 10 iterations
  if(b %% 10 == 0){
    temp_file <- file.path(
      raw_dir,
      paste0("TEMP_", model_key, "_p", p_dim, "_n", n, "_iter", b, ".rds")
    )
    
    saveRDS(
      list(
        rel = dplyr::bind_rows(rel_list),
        fpr_tpr = dplyr::bind_rows(fpr_tpr_list),
        pd = dplyr::bind_rows(pd_list),
        setting = list(
          array_id = array_id,
          model_key = model_key,
          model = model_name,
          n = n,
          p_dim = p_dim,
          iter_done = b,
          eps = eps,
          tau = tau
        )
      ),
      temp_file
    )
  }
}

all_rel <- bind_rows(rel_list) %>%
  mutate(
    norm = factor(norm, levels = norm_levels),
    method = factor(method, levels = method_levels)
  )

all_fpr_tpr <- bind_rows(fpr_tpr_list) %>%
  mutate(method = factor(method, levels = method_levels))

all_pd <- bind_rows(pd_list) %>%
  mutate(method = factor(method, levels = method_levels))

################################################################################
## 4. Summaries
################################################################################

all_rel <- bind_rows(rel_list) %>%
  mutate(
    norm = factor(norm, levels = norm_levels),
    method = factor(method, levels = method_levels)
  )

all_fpr_tpr <- bind_rows(fpr_tpr_list) %>%
  mutate(method = factor(method, levels = method_levels))

all_pd <- bind_rows(pd_list) %>%
  mutate(method = factor(method, levels = method_levels))

## Relative error: mean (sd)
rel_summary <- all_rel %>%
  group_by(model_key, model, n, p_dim, method, norm) %>%
  summarise(
    mean_re = mean(rel_error, na.rm = TRUE),
    sd_re   = sd(rel_error, na.rm = TRUE),
    rel_result = sprintf("%.3f (%.3f)", mean_re, sd_re),
    .groups = "drop"
  )

# rel_table <- rel_summary %>%
#   dplyr::select(model_key, model, n, p_dim, method, norm, rel_result) %>%
#   tidyr::pivot_wider(
#     names_from = norm,
#     values_from = rel_result
#   )

## FPR / TPR: mean (sd)
fpr_tpr_summary <- all_fpr_tpr %>%
  group_by(model_key, model, n, p_dim, method) %>%
  summarise(
    mean_FPR = mean(FPR, na.rm = TRUE),
    sd_FPR   = sd(FPR, na.rm = TRUE),
    mean_TPR = mean(TPR, na.rm = TRUE),
    sd_TPR   = sd(TPR, na.rm = TRUE),
    FPR = sprintf("%.3f (%.3f)", mean_FPR, sd_FPR),
    TPR = sprintf("%.3f (%.3f)", mean_TPR, sd_TPR),
    .groups = "drop"
  )

## PD: min eigenvalue > eps
pd_summary <- all_pd %>%
  group_by(model_key, model, n, p_dim, method) %>%
  summarise(
    PD_count = sum(is_eps_pd, na.rm = TRUE),
    PD_rate  = mean(is_eps_pd, na.rm = TRUE),
    mean_min_eig = mean(min_eig, na.rm = TRUE),
    sd_min_eig   = sd(min_eig, na.rm = TRUE),
    PD = sprintf("%d/%d", PD_count, iter),
    PD_rate_result = sprintf("%.3f", PD_rate),
    min_eig_result = sprintf("%.3f (%.3f)", mean_min_eig, sd_min_eig),
    .groups = "drop"
  )

## Final summary table
# table_final <- rel_table %>%
#   left_join(
#     fpr_tpr_summary %>%
#       dplyr::select(model_key, model, n, p_dim, method, FPR, TPR),
#     by = c("model_key", "model", "n", "p_dim", "method")
#   ) %>%
#   left_join(
#     pd_summary %>%
#       dplyr::select(
#         model_key, model, n, p_dim, method,
#         PD, PD_rate = PD_rate_result, min_eig = min_eig_result
#       ),
#     by = c("model_key", "model", "n", "p_dim", "method")
#   ) %>%
#   mutate(
#     Setting = paste0(model_key, "/", n, "/", p_dim)
#   ) %>%
#   dplyr::select(
#     Setting,
#     model_key,
#     model,
#     n,
#     p_dim,
#     method,
#     Frobenius,
#     Spectral,
#     `Element-wise max`,
#     `L-infinity`,
#     FPR,
#     TPR,
#     PD,
#     PD_rate,
#     min_eig
#   ) %>%
#   arrange(model_key, p_dim, method)
# 
# print(table_final)

################################################################################
## 5. Save
################################################################################

saveRDS(
  list(
    rel = all_rel,
    fpr_tpr = all_fpr_tpr,
    pd = all_pd,
    setting = list(
      array_id = array_id,
      model_key = model_key,
      model = model_name,
      n = n,
      p_dim = p_dim,
      iter = iter,
      eps = eps,
      tau = tau
    )
  ),
  output_file
)

save(
  all_rel,
  all_fpr_tpr,
  all_pd,
  array_id,
  model_key,
  model_name,
  n,
  p_dim,
  iter,
  eps,
  tau,
  method_levels,
  norm_levels,
  model_levels,
  file = workspace_file
)

cat("\nFinished and saved:\n")
cat(output_file, "\n")
cat(workspace_file, "\n")
