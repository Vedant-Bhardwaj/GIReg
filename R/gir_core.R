
GIR_R2 <- function(y, M, cdf_hat, theta_hat){
  Lplus_hat <- y*(((2*cdf_hat - 1) - M%*%theta_hat)^2)
  C <- (y*(2*cdf_hat - 1))/mean(y)
  Lplus_tilde <- y*(((2*cdf_hat - 1) - C)^2)
  return(1 - (mean(Lplus_hat)/mean(Lplus_tilde)))
}
compute_theta_hat <- function(y,X1, X2, tau_seq, c, ew,basis,basis_type,method="qr",alpha=c(0.01,0.05,0.1)) {


  # number of observations
  obs <- nrow(X2)
  # build the expanded design
  #M <- build_design_matrix(X)
  M<-as.matrix(X2)
  #print(M)
  print("build")
  # (A) fit quantile regressions over tau_seq
  #fits <- fit_qr_multi(M, y, tau_seq)
  #M <- build_design_matrix(X_raw)
  if (method == "cqr") {

    fits_qr<- est_cqr(y, c = c, x = X1, w = ew,
                      nqreg = length(tau_seq),
                      c1 = 0.1, c2 = 0.05, nsteps = 3)

  } else {   method == "qr"

    fits_qr <- fit_qr_multi(X1, y, tau_grid = tau_seq,weights = ew, basis,basis_type)
  }


  print("FIts")
  # (B) predict the CDF at the observed points
  cdf_hat <- predict_cdf_qr_multi(fits = fits_qr$fits, X_new = fits_qr$X_poly[, -1, drop = FALSE], y)

  #print("FIts")
  # (B) predict the CDF at the observed points

  #cdf_hat <- predict_cdf_qr_multi(fits = fits_qr$fits, X_new = fits_qr$X_poly, y0 = y)

  # (C) compute theta_hat
  #    solve( t(M * y) %*% M ) %*% colSums( M * y * (2*cdf_hat - 1) )
  # use crossprod for readability:
  y<-y*ew
  theta_hat <- solve( t(M * y) %*% M) %*%
    (colSums( M * y * (2 * cdf_hat - 1) ))
  names(theta_hat)<-colnames(M)

  gir_r2 <- GIR_R2(y = y, M = M, cdf_hat = cdf_hat, theta_hat = theta_hat)

  list(
    theta_hat = as.numeric(theta_hat),  # keep it numeric vector-like if you prefer
    gir_r2 = as.numeric(gir_r2),
    cdf_hat = cdf_hat   # optional: drop if you don’t need it outside
  )
}

