library(Matrix)

# Model 1: sigma_ij = (1 - |i - j|/10)+
cov.mtx.m1 <- function(p) {
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- max(0, 1 - abs(i - j) / 10)
    }
  }
  return(Sigma)
}
# Model 2: Block matrix
cov.mtx.m2 <- function(p) {
  K <- p / 10  
  if (K != floor(K)) {
    stop("p must be splitted with 10")
  }
  subset_size <- 10
  I_k_max <- seq(10, p, by = 10)  
  
  Sigma <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      indicator1 <- sum(sapply(I_k_max, function(max_idx) {
        (i <= max_idx && j <= max_idx) && ((i > max_idx - subset_size) && (j > max_idx - subset_size))
      }))
      indicator2 <- sum(sapply(I_k_max[-length(I_k_max)], function(max_idx) {
        (i == max_idx + 1 && j == max_idx + 1)
      }))
      
      Sigma[i, j] <- 0.6 * (i == j) + 0.4 * indicator1 + 0.4 * indicator2
    }
  }
  
  return(Sigma)
}
# Model 3: AR(1)
cov.mtx.m3 <- function(p, rho = 0.7) {
  mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      mat[i, j] <- rho^abs(i - j)
    }
  }
  return(mat)
}
# Model 4: Sparse model
cov.mtx.m4 <- function(p, q = 0.2) {
  B <- matrix(0, p, p)
  for (i in 1:(p - 1)) {
    rho <- runif(1, 0.3, 0.8)
    z   <- rbinom(p - i, 1, q)          # Bern(q)
    B[i, (i + 1):p] <- rho * z          # upper
  }
  B <- B + t(B)                            # symmetrize
  mat <- diag(1, p) + B                  # I + B
  
  ## epsilon = max(-lambda_min(I+B), 0) + 0.01
  lam_min <- min(eigen(mat, symmetric = TRUE, only.values = TRUE)$values)
  eps <- max(-lam_min, 0) + 0.01
  mat <- mat + eps * diag(p)
  
  
  
  return(mat)
}

#visualization------------------------------------------------------------------
p <- 100
Sigma1 <- cov.mtx.m1(p)
Sigma2 <- cov.mtx.m2(p)
Sigma3 <- cov.mtx.m3(p)
Sigma4 <- cov.mtx.m4(p)

image(Sigma1, main = "Model 1 Covariance Matrix")
image(Sigma2, main = "Model 2 Covariance Matrix")
image(Sigma3, main = "Model 3 Covariance Matrix")
image(Sigma4, main = "Model 4 Covariance Matrix")
