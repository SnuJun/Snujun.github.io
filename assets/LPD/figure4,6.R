################################################################################
## 0. 기본 세팅
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)

wd <- "C:/Users/user/Desktop/LPD 0413"
setwd(wd)

source("data generate.R")
source("Estimators.R")
source("getLPD.R")

n <- 100
p <- 120
eps <- 1e-4
BIG_MU <- 12

################################################################################
## 1. 모델 리스트
################################################################################

cov_list <- list(
  m1 = cov.mtx.m1,
  m2 = cov.mtx.m2,
  m3 = cov.mtx.m3,
  m4 = cov.mtx.m4
)

model_names <- c("Banded", "Block", "Toeplitz", "Sparse")
model_levels <- c("Banded", "Block", "Toeplitz", "Sparse")

################################################################################
## 2. 공통 initial estimator 생성
################################################################################

init_list <- list()

set.seed(1234)

for(i in 1:4){
  
  Sigma <- cov_list[[i]](p)
  
  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma
  )
  
  A <- hard_thresholding_CV(X)$mat
  
  init_list[[i]] <- A
}

################################################################################
## 3. mu 계산 함수
################################################################################

get_mu_info <- function(A, eps){
  
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lam1 <- min(eig)
  lamp <- max(eig)
  bar_lam <- mean(eig)
  
  mu_spec <- max(eps, (lam1 + lamp) / 2)
  
  mu_frob <- bar_lam + sum((eig - bar_lam)^2) / (length(eig) * (bar_lam - lam1))
  
  M1 <- max(rowSums(abs(A)))
  
  A_tmp <- A - 2 * diag(diag(A))
  A_tmp[row(A_tmp) != col(A_tmp)] <- abs(A_tmp[row(A_tmp) != col(A_tmp)])
  M2 <- max(rowSums(A_tmp))
  
  if(lam1 + M2 > 0){
    mu_linf_type <- "inf"
    mu_linf <- NA_real_
  } else if(abs(lam1 + M2) < 1e-8){
    mu_linf_type <- "interval"
    mu_linf <- (M1 - M2) / 2
  } else{
    mu_linf_type <- "point"
    mu_linf <- (M1 - M2) / 2
  }
  
  a_diag <- diag(A)
  a_max <- max(a_diag)
  a_min <- min(a_diag)
  a_off <- max(abs(A[row(A) != col(A)]))
  
  if((a_max - a_min) / 2 > a_off){
    mu_ele <- (a_max + a_min) / 2
  } else{
    mu_ele <- a_min + a_off
  }
  
  return(list(
    spec = as.numeric(mu_spec),
    frob = as.numeric(mu_frob),
    linf = as.numeric(mu_linf),
    linf_type = mu_linf_type,
    ele = as.numeric(mu_ele)
  ))
}

################################################################################
## 4. mu 정보 정리
################################################################################

plot_data <- data.frame()

for(i in 1:4){
  
  A <- init_list[[i]]
  mu_info <- get_mu_info(A, eps)
  
  tmp <- data.frame(
    model = model_names[i],
    
    spec_start = mu_info$spec,
    spec_end   = BIG_MU,
    
    linf       = mu_info$linf,
    linf_type  = mu_info$linf_type,
    
    frob = mu_info$frob,
    ele  = mu_info$ele
  )
  
  plot_data <- rbind(plot_data, tmp)
}

plot_data <- plot_data %>%
  mutate(
    model = factor(model, levels = model_levels),
    linf_plot = ifelse(linf_type == "inf", BIG_MU, linf),
    linf_start = ifelse(linf_type == "interval", linf, NA_real_),
    linf_end = ifelse(linf_type == "interval", BIG_MU, linf_plot)
  )

print(plot_data)

################################################################################
## 5. Figure 1: error-optimal mu 시각화
################################################################################

spec_seg <- plot_data %>%
  transmute(
    model,
    norm = "Spectral",
    start = spec_start,
    end = spec_end
  )

linf_seg <- plot_data %>%
  filter(linf_type == "interval") %>%
  transmute(
    model,
    norm = "L-infinity",
    start = linf_start,
    end = linf_end
  )

seg_df <- bind_rows(spec_seg, linf_seg)

point_df <- plot_data %>%
  transmute(
    model,
    Frobenius = frob,
    `Element-wise max` = ele
  ) %>%
  pivot_longer(
    cols = c(Frobenius, `Element-wise max`),
    names_to = "norm",
    values_to = "mu"
  )

linf_point_df <- plot_data %>%
  filter(linf_type %in% c("point", "inf")) %>%
  transmute(
    model,
    norm = "L-infinity",
    mu = linf_plot
  )

point_df <- bind_rows(point_df, linf_point_df)

p1 <- ggplot() +
  
  geom_segment(
    data = seg_df,
    aes(x = start, xend = end, y = norm, yend = norm, color = norm),
    size = 1.2
  ) +
  
  geom_point(
    data = spec_seg,
    aes(x = start, y = norm, color = norm, shape = norm),
    size = 3
  ) +
  
  geom_point(
    data = point_df,
    aes(x = mu, y = norm, color = norm, shape = norm),
    size = 3
  ) +
  

  
  facet_wrap(~ model, nrow = 2) +
  
  scale_y_discrete(
    limits = rev(c("Spectral", "Frobenius", "L-infinity", "Element-wise max"))
  ) +
  
  scale_color_manual(
    values = c(
      "Spectral" = "red",
      "Frobenius" = "darkgreen",
      "L-infinity" = "blue",
      "Element-wise max" = "purple"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Spectral" = 16,
      "Frobenius" = 15,
      "L-infinity" = 17,
      "Element-wise max" = 18
    )
  ) +
  
  labs(
    title = expression("Optimal " * mu * " minimizing " * 
                         "||A-" * Phi(A) * "||"),
    x = expression(mu),
    y = NULL,
    color = "Norm",
    shape = "Norm"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p1)

################################################################################
## 6. fixed mu LPD 함수
################################################################################

get_LPD_fixed_mu <- function(A, mu, eps = 1e-4){
  
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  lam1 <- min(eig)
  
  if(mu <= eps || abs(mu - lam1) < 1e-10){
    return(NULL)
  }
  
  alpha <- (mu - eps) / (mu - lam1)
  
  if(alpha <= 0 || alpha >= 1){
    return(NULL)
  }
  
  LPD <- alpha * A + (1 - alpha) * mu * diag(nrow(A))
  
  return(LPD)
}

################################################################################
## 7. Spectral norm error curve 계산
################################################################################

mu_grid <- seq(0, BIG_MU, length.out = 700)

error_curve_df <- data.frame()
vline_df <- data.frame()

for(i in 1:4){
  
  A <- init_list[[i]]
  mu_info <- get_mu_info(A, eps)
  
  mu_S <- mu_info$spec
  
  tmp_curve <- lapply(mu_grid, function(mu){
    
    LPD <- get_LPD_fixed_mu(A, mu, eps)
    
    if(is.null(LPD)){
      err <- NA_real_
    } else{
      err <- norm(A - LPD, type = "2")
    }
    
    data.frame(
      model = model_names[i],
      mu = mu,
      error = err
    )
  }) %>% bind_rows()
  
  best_row <- tmp_curve %>%
    filter(!is.na(error)) %>%
    slice_min(error, n = 1, with_ties = FALSE)
  
  mu_inf_plot <- ifelse(mu_info$linf_type == "inf", BIG_MU, mu_info$linf)
  
  tmp_vline <- data.frame(
    model = model_names[i],
    type = c("Empirical", "mu_S", "mu_inf"),
    mu = c(best_row$mu, mu_S, mu_inf_plot)
  )
  
  error_curve_df <- bind_rows(error_curve_df, tmp_curve)
  vline_df <- bind_rows(vline_df, tmp_vline)
}

error_curve_df$model <- factor(error_curve_df$model, levels = model_levels)
vline_df$model <- factor(vline_df$model, levels = model_levels)

vline_df$type <- factor(
  vline_df$type,
  levels = c("Empirical", "mu_S", "mu_inf")
)

################################################################################
## 7. Norm-specific error curves 계산
################################################################################

mu_grid <- seq(0, BIG_MU, length.out = 700)

get_error_by_norm <- function(A, LPD, norm_type){
  if(norm_type == "Spectral"){
    return(norm(A - LPD, type = "2"))
  } else if(norm_type == "Frobenius"){
    return(norm(A - LPD, type = "F") / sqrt(nrow(A)))  # scaled Frobenius
  } else if(norm_type == "L-infinity"){
    return(norm(A - LPD, type = "I"))
  } else if(norm_type == "Element-wise max"){
    return(max(abs(A - LPD)))
  }
}

error_curve_df <- data.frame()
vline_df <- data.frame()

for(i in 1:4){
  
  A <- init_list[[i]]
  mu_info <- get_mu_info(A, eps)
  
  mu_linf_plot <- ifelse(mu_info$linf_type == "inf", BIG_MU, mu_info$linf)
  
  mu_theory_df <- data.frame(
    model = model_names[i],
    norm_type = c("Spectral", "Frobenius", "L-infinity", "Element-wise max"),
    mu_opt = c(mu_info$spec, mu_info$frob, mu_linf_plot, mu_info$ele)
  )
  
  tmp_curve <- lapply(mu_grid, function(mu){
    
    LPD <- get_LPD_fixed_mu(A, mu, eps)
    
    if(is.null(LPD)){
      return(data.frame(
        model = model_names[i],
        norm_type = c("Spectral", "Frobenius", "L-infinity", "Element-wise max"),
        mu = mu,
        error = NA_real_
      ))
    }
    
    data.frame(
      model = model_names[i],
      norm_type = c("Spectral", "Frobenius", "L-infinity", "Element-wise max"),
      mu = mu,
      error = c(
        get_error_by_norm(A, LPD, "Spectral"),
        get_error_by_norm(A, LPD, "Frobenius"),
        get_error_by_norm(A, LPD, "L-infinity"),
        get_error_by_norm(A, LPD, "Element-wise max")
      )
    )
  }) %>% bind_rows()
  
  error_curve_df <- bind_rows(error_curve_df, tmp_curve)
  vline_df <- bind_rows(vline_df, mu_theory_df)
}

error_curve_df <- error_curve_df %>%
  mutate(
    model = factor(model, levels = model_levels),
    norm_type = factor(
      norm_type,
      levels = c("Spectral", "Frobenius", "L-infinity", "Element-wise max")
    )
  )

vline_df <- vline_df %>%
  mutate(
    model = factor(model, levels = model_levels),
    norm_type = factor(
      norm_type,
      levels = c("Spectral", "Frobenius", "L-infinity", "Element-wise max")
    )
  )

################################################################################
## Figure 2: Norm-specific error curves over mu
## 각각의 norm별로 따로 출력
################################################################################

make_error_plot <- function(target_norm){
  
  df_curve <- error_curve_df %>%
    filter(norm_type == target_norm)
  
  df_vline <- vline_df %>%
    filter(norm_type == target_norm)
  
  p <- ggplot(df_curve, aes(x = mu, y = error)) +
    geom_line(size = 0.9, color = "black") +
    geom_vline(
      data = df_vline,
      aes(xintercept = mu_opt),
      color = "red",
      linetype = "dashed",
      size = 0.9
    ) +
    facet_wrap(~ model, nrow = 2, scales = "free_y") +
    labs(
      title = paste0(target_norm, " norm error curve over mu"),
      x = expression(mu),
      y = expression("||A - " * Phi[mu](A) * "||")
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ## arrow only for Spectral norm
  if(target_norm == "Spectral"){
    
    arrow_df <- df_curve %>%
      filter(!is.na(error)) %>%
      group_by(model) %>%
      summarise(
        y_arrow =
          min(error, na.rm = TRUE) +
          0.20 * (max(error, na.rm = TRUE) - min(error, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      left_join(
        df_vline %>% dplyr::select(model, mu_opt),
        by = "model"
      ) %>%
      mutate(
        x_start = mu_opt,
        x_end   = mu_opt + 1.2
      )
    
    p <- p +
      geom_segment(
        data = arrow_df,
        aes(
          x = x_start,
          xend = x_end,
          y = y_arrow * 1.5,
          yend = y_arrow * 1.5
        ),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.25, "cm")),
        color = "red",
        size = 1.3
      )
  }
  
  return(p)
}

p2_spec <- make_error_plot("Spectral")
p2_frob <- make_error_plot("Frobenius")
p2_linf <- make_error_plot("L-infinity")
p2_ele  <- make_error_plot("Element-wise max")

print(p2_spec)
print(p2_frob)
print(p2_linf)
print(p2_ele)




# install.packages("patchwork")  # 처음 한 번만
library(patchwork)

p_all <- (p2_spec | p2_frob) /
  (p2_linf | p2_ele)

print(p_all)


################################################################################
## Figure 3: convergence-rate-preserving mu 시각화 (∞ 표시 제거)
################################################################################

spec_rate_seg <- plot_data %>%
  transmute(
    model,
    norm = "Spectral",
    start = spec_start,
    end = BIG_MU
  )

frob_rate_seg <- plot_data %>%
  transmute(
    model,
    norm = "Frobenius",
    start = frob,
    end = BIG_MU
  )

linf_rate_seg <- plot_data %>%
  filter(linf_type == "interval") %>%
  transmute(
    model,
    norm = "L-infinity",
    start = linf_start,
    end = BIG_MU
  )

rate_seg_df <- bind_rows(
  spec_rate_seg,
  frob_rate_seg,
  linf_rate_seg
)

rate_start_points <- bind_rows(
  spec_rate_seg %>% transmute(model, norm, mu = start),
  frob_rate_seg %>% transmute(model, norm, mu = start)
)

linf_rate_point <- plot_data %>%
  filter(linf_type %in% c("point", "inf")) %>%
  transmute(
    model,
    norm = "L-infinity",
    mu = linf_plot
  )

ele_rate_point <- plot_data %>%
  transmute(
    model,
    norm = "Element-wise max",
    mu = spec_start
  )

rate_point_df <- bind_rows(
  rate_start_points,
  linf_rate_point,
  ele_rate_point
)

rate_seg_df$model <- factor(rate_seg_df$model, levels = model_levels)
rate_point_df$model <- factor(rate_point_df$model, levels = model_levels)

rate_seg_df$norm <- factor(
  rate_seg_df$norm,
  levels = c("Spectral", "Frobenius", "L-infinity", "Element-wise max")
)

rate_point_df$norm <- factor(
  rate_point_df$norm,
  levels = c("Spectral", "Frobenius", "L-infinity", "Element-wise max")
)

p_rate <- ggplot() +
  
  geom_segment(
    data = rate_seg_df,
    aes(x = start, xend = end, y = norm, yend = norm, color = norm),
    size = 1.2
  ) +
  
  geom_point(
    data = rate_point_df,
    aes(x = mu, y = norm, color = norm, shape = norm),
    size = 3
  ) +
  
  facet_wrap(~ model, nrow = 2) +
  
  scale_y_discrete(
    limits = rev(c("Spectral", "Frobenius", "L-infinity", "Element-wise max"))
  ) +
  
  scale_color_manual(
    values = c(
      "Spectral" = "red",
      "Frobenius" = "darkgreen",
      "L-infinity" = "blue",
      "Element-wise max" = "purple"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Spectral" = 16,
      "Frobenius" = 15,
      "L-infinity" = 17,
      "Element-wise max" = 18
    )
  ) +
  
  labs(
    title = expression("Rate-preserving " * mu *
                         " for " *
                         "||" * Phi[mu](hat(Sigma)) - Sigma * "||"),
    x = expression(mu),
    y = NULL,
    color = "Norm",
    shape = "Norm"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_rate)
