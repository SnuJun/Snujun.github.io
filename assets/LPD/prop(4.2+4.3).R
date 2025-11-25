#(1) Eigenvalue difference -----------------------------------------------------
set.seed(1)

p_set <- c(50, 100, 150, 200, 250)
n <- 100

true_cov_list <- list(cov.mtx.m1, cov.mtx.m2, cov.mtx.m3,cov.mtx.m4)
true_cov_names <- c("banded", "block", "toeplitz","sparse") 

est_list <- list(
  "hard" = function(X) hard_thresholding_CV(X)$mat,
  "blockwise" = function(X) ftn.cov.blck(X)$estimate,
  "banding" = function(X) ftn.cov.band(X)$estimate,
  "poet" = function(X) POET(t(X), K = 5, thres = "hard")$SigmaY  # note t(X)!
)

results <- list()

for (p in p_set) {
  cat(">> p =", p, "\n")
  for (i in seq_along(true_cov_list)) {
    covgen <- true_cov_list[[i]]
    Sigma <- covgen(p)
    
    # producing samples
    data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    # esimation and evaluation
    for (name in names(est_list)) {
      estimator <- est_list[[name]]
      Sigma_hat <- estimator(data)
      
      # evaluation
      ratio <- eval_ftn_1(Sigma, Sigma_hat)$ratio
      
      # storing data
      results[[length(results) + 1]] <- list(
        p = p,
        model = true_cov_names[i],
        method = name,
        eval_ratio = ratio
      )
    }
  }
}


df_results <- do.call(rbind, lapply(results, as.data.frame))
print(df_results)


library(dplyr)
df_summary <- df_results %>%
  group_by(p, model, method) %>%
  summarise(mean_ratio = mean(eval_ratio), .groups = "drop")

print(df_summary)

library(ggplot2)
library(dplyr)


# 팩터 순서 정렬 (선택)
df_results$model <- factor(df_results$model, levels = c("banded", "block", "toeplitz", "sparse"))
df_results$method <- factor(df_results$method, levels = c("hard", "blockwise", "banding", "poet"))

# 시각화
ggplot(df_results, aes(x = p, y = eval_ratio, group = method)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue") +
  facet_grid(method ~ model) +
  theme_bw(base_size = 13) +
  labs(
    title = "Eigenvalue Difference Ratio",
    x = "Dimension p",
    y = "Evaluation Ratio"
  ) +
  ylim(0, 3) + 
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "gray70")
  )



#(2) Optimal rate over elementmaxwise norm -------------------------------------
library(MASS)

set.seed(1)

p_set <- c(50, 100, 150, 200, 250)
n <- 100

true_cov_list <- list(cov.mtx.m1, cov.mtx.m2, cov.mtx.m3cov.mtx.m4) 
true_cov_names <- c("banded", "block", "toeplitz", "sparse")

est_list <- list(
  "hard" = function(X) hard_thresholding_CV(X)$mat,
  "blockwise" = function(X) ftn.cov.blck(X)$estimate,
  "banding" = function(X) ftn.cov.band(X)$estimate,
  "poet" = function(X) POET(t(X), K = 5, thres = "hard")$SigmaY
)

results_lpd <- list()

for (p in p_set) {
  for (i in seq_along(true_cov_list)) {
    covgen <- true_cov_list[[i]]
    Sigma <- covgen(p)
    
    # Generate sample
    data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    
    for (name in names(est_list)) {
      estimator <- est_list[[name]]
      Sigma_hat <- estimator(data)
      
      # LPD projection
      LPD_result <- getLPD(Sigma_hat, lwr_bd = 1e-3, type_norm = "elemax")
      Sigma_lpd <- LPD_result$LPD
      
      # Ratio evaluation
      ratio <- eval_ftn_2(Sigma, Sigma_hat, Sigma_lpd)
      
      # Store
      results_lpd[[length(results_lpd) + 1]] <- list(
        p = p,
        model = true_cov_names[i],
        method = name,
        ratio = ratio
      )
    }
  }
}

df_lpd <- do.call(rbind, lapply(results_lpd, as.data.frame))


library(ggplot2)

df_lpd$model <- factor(df_lpd$model, levels = c("banded", "block", "toeplitz"))#, "sparse"
df_lpd$method <- factor(df_lpd$method, levels = c("hard", "blockwise", "banding", "poet"))

ggplot(df_lpd, aes(x = p, y = ratio)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  facet_grid(method ~ model) +
  theme_bw(base_size = 13) +
  labs(
    title = "Ratio of Errors(LPD)",
    x = "Dimension p",
    y = "‖LPD - Σ‖ / ‖Estimator - Σ‖"
  ) +
  coord_cartesian(ylim = c(0, 1.5)) 
