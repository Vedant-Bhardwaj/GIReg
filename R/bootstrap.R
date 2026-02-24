
# 2) Fully parallel bootstrap (pairs bootstrap)
bootstrap_theta <- function(y, X1,X2, c, w, tau_seq,
                            B = 200,ncores = max(1L, parallel::detectCores() - 1L), basis="exper",basis_type = "spline",method = c("cqr", "qr"),alpha=c(0.05,0.1,0.01)) {
  stopifnot(nrow(X2) == length(y), length(y) == length(c), length(y) == length(w))
  stopifnot(nrow(X1) == length(y), length(y) == length(c), length(y) == length(w))
  
  fit0 <- compute_theta_hat(y, X1,X2, tau_seq, c, w, basis,basis_type,method=method )
  theta_hat <- fit0$theta_hat
  print(theta_hat)
  gir_r2    <- fit0$gir_r2
  theta_col <- matrix(theta_hat, ncol = 1)
  rownames(theta_col) <- colnames(as.matrix(X2))
  colnames(theta_col) <- "Coefficients"
  print(theta_col)
  cat("GIR_R2 (original sample):", gir_r2, "\n")
  p <- length(theta_hat)
  theta_boot <- matrix(NA_real_, nrow = B, ncol = p, dimnames = list(NULL, rownames(theta_col)))
  
  # for (b in 1:B) {
  #   set.seed(202 + b)  
  #   idx <- sample.int(length(y), replace = TRUE)
  #   #print(idx)
  #  
  #   theta_boot[b, ] <- compute_theta_hat(
  #     y[idx], X[idx, , drop = FALSE], tau_seq, c[idx], w[idx],
  #     ncores_est = ncores_est  # parallel **within** each replicate over taus
  #   )
  # }
  
  cl <- parallel::makeCluster(ncores/2-2)
  #cl <- parallel::makeCluster(4)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(3)
  parallel::clusterEvalQ(cl, {
    library(quantreg); library(splines)
    library(fda); library(orthopolynom)  # uncomment if you use those bases
  })
  
  theta_boot <- foreach::foreach(
    b = 1:B, .combine = cbind,
   # .errorhandling = "pass",
    .packages = c("quantreg"),   # add more if needed
    .export   = c("compute_theta_hat","fit_qr_multi","predict_cdf_qr_multi","GIR_R2",
                  "poly_features","legendre_features","fourier_features","spline_features","est_cqr","est_qr","wq","quantreg","predict.qr_like")
  ) %dopar% {
    # multiplier bootstrap
    w_b <- w * rexp(length(w), rate = 1)
    compute_theta_hat(y, X1,X2, tau_seq, c, w_b,
                      basis = basis, basis_type = basis_type,method=method )$theta_hat
  }
  theta_boot <-t(theta_boot )
  print(theta_boot)
  colnames(theta_boot) <- rownames(theta_col)

  # your requested layout: B × (2p)
  results <- cbind(
    theta_hat = matrix(rep(theta_hat, each = B), nrow = B,
                       dimnames = list(NULL, names(theta_hat))),
    theta_boot = theta_boot
  )
  #print(B)
 # print(p)
  # “basic/root” CI exactly like you wrote
  #print(matrix(theta_hat, nrow = B, ncol = p, byrow = TRUE))
 #print(theta_boot)
 #print(results)
  root   <- theta_boot - matrix(theta_hat, nrow = B, ncol = p, byrow = TRUE)
  # Assume: 'alpha' is a numeric vector like c(0.10, 0.05)
  #          meaning you want 90% and 95% confidence intervals.
  print(names(theta_hat))
  q_root_all <- lapply(alpha, function(a) {
    probs <- c(a / 2, 1 - a / 2)
    q_root <- apply(root, 2, quantile, probs = probs, na.rm = TRUE)
    
    lower_q <- as.numeric(theta_hat - q_root[2, ])
    upper_q <- as.numeric(theta_hat - q_root[1, ])
    names(lower_q) <- names(upper_q) <- names(theta_hat)
    
    ci <- rbind(lower = lower_q, upper = upper_q)
    colnames(ci) <- colnames(as.matrix(X2))  # force it
    ci
  })
  
  # Name each CI by confidence level, e.g. ci_basic_90, ci_basic_95
  names(q_root_all) <- paste0("ci_basic_", 100 * (1 - alpha))
  
  # Standard error
  se_boot <- apply(theta_boot, 2, sd)
  
  # Final output
  list(
    theta_hat = theta_col,
    theta_boot = theta_boot,
    #results = results,
    #se_boot = se_boot,
    ci_basic = q_root_all,
    R_squared=gir_r2, 
    n = length(y)
    
  )
  
}
