fit_qr_multi <- function(
    X, y, tau_grid,
    weights = rep(1, nrow(as.data.frame(X))),
    basis = NULL,
    basis_type = c("none","poly","legendre","fourier","spline"),
    basis_args = list()
) {
  X_df <- as.data.frame(X)
  n    <- nrow(X_df); p <- ncol(X_df)
  stopifnot(length(y) == n, length(weights) == n)
  basis_type <- match.arg(basis_type)

  # ---- resolve which columns to expand ----
  if (is.null(basis) || length(basis) == 0 || basis_type == "none") {
    X_poly <- X_df
  } else {
    if (is.character(basis)) {
      miss <- setdiff(basis, names(X_df))
      if (length(miss)) stop("basis columns not found in X: ", paste(miss, collapse=", "))
      basis_idx <- match(basis, names(X_df))
    } else {
      basis_idx <- as.integer(basis)
      if (any(basis_idx < 1 | basis_idx > ncol(X_df))) stop("basis indices out of bounds.")
    }

    X_keep   <- X_df[, setdiff(seq_len(ncol(X_df)), basis_idx), drop = FALSE]
    X_expand <- X_df[, basis_idx, drop = FALSE]

    defaults <- switch(
      basis_type,
      poly     = list(degree = max(1L, as.integer(sqrt(n) / (2 * max(1L, p))))),
      legendre = list(degree = max(1L, as.integer(sqrt(n) / max(1L, p)))),
      fourier  = list(nbasis = 2L),
      spline   = list(df     = max(3L, as.integer(sqrt(n) / (2 * max(1L, p))))),
      list()
    )
    call_args <- utils::modifyList(defaults, basis_args)

    B <- switch(
      basis_type,
      poly     = do.call(poly_features,     c(list(as.matrix(X_expand)), call_args)),
      legendre = do.call(legendre_features, c(list(as.matrix(X_expand)), call_args)),
      fourier  = do.call(fourier_features,  c(list(as.matrix(X_expand)), call_args)),
      spline   = do.call(spline_features,   c(list(as.matrix(X_expand)), call_args)),
      stop("Unknown basis_type: ", basis_type)
    )
    B <- as.data.frame(B)

    if (is.null(colnames(B))) {
      colnames(B) <- paste0(basis_type, "_", seq_len(ncol(B)))
    } else {
      colnames(B) <- paste0(basis_type, "_", colnames(B))
    }

    X_poly <- cbind(X_keep, B)
  }

  df <- data.frame(y = y, X_poly)

  # ---- serial over taus ----
  fits <- lapply(tau_grid, function(tau) {
    quantreg::rq(y ~ ., data = df, tau = tau, weights = weights)
  })

  names(fits) <- paste0("tau", tau_grid)
  list(fits = fits, X_poly = X_poly)
}


predict_cdf_qr_multi <- function(fits, X_new, y0) {
  # must expand X_new the same way
  #X_new_poly <- poly_features(X_new, degree = as.integer(sqrt(obs)/(2*p)))
  #X_new_poly<- legendre_features(X, degree = as.integer(sqrt(obs)/p))
  #X_new_poly<- fourier_features(X, nbasis  = 2)


  #X_new_poly <-spline_features(X, df = as.integer(sqrt(obs)/(2*p)))
  #X_new_poly <-add_interactions_2_3_4(X)
  X_new_poly <-X_new

  n_pts    <- nrow(X_new_poly)

  # matrix of quantile predictions: rows = tau, cols = points
  Q <- sapply(fits, function(f) {
    predict(f, newdata = data.frame(X_new_poly))
  })
  Q <- t(Q)  # now length(tau_grid) × n_pts
  # print(fits)
  if (length(y0) == 1) y0 <- rep(y0, n_pts)
  #print(Q)
  #print(y0)
  cdf <- vapply(seq_len(n_pts), function(i) {
    mean(y0[i] > Q[, i])
  }, numeric(1))
  #print(cdf)
  return(cdf)
}
