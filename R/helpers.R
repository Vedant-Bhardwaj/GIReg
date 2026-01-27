
# ——————————————
# 0) Utility: Expand each column of X into 1st, 2nd and 3rd powers
poly_features <- function(X, degree = 5) {
  X <- as.data.frame(X)
  polys <- lapply(names(X), function(var) {
    x <- X[[var]]
    # create a data.frame of x^1, x^2, x^3
    df <- sapply(1:degree, function(d) x^d)
    colnames(df) <- paste0(var, "_", 1:degree)
    as.data.frame(df)
  })
  do.call(cbind, polys)
}

legendre_features <- function(X,
                              degree       = 5,
                              normalized   = TRUE,
                              include_P0   = FALSE) {

  stopifnot(degree >= 1,
            is.data.frame(X) || is.matrix(X))



  X <- as.data.frame(X)

  ## ----------------------------------------------------------------
  ## Pre‑compute the list of Legendre polynomials (as functions)
  ## ----------------------------------------------------------------
  # legendre.polynomials() returns a list where element [[k+1]] is P_k
  P_list <- orthopolynom::legendre.polynomials(n = degree,
                                               normalized = normalized)

  # Keep P0 if requested; otherwise drop it
  idx_needed <- seq.int(from = if (include_P0) 0 else 1,
                        to   = degree)
  P_funs <- lapply(P_list[idx_needed + 1L], as.function)

  ## ----------------------------------------------------------------
  ## For each variable: evaluate every P_k at that column
  ## ----------------------------------------------------------------
  polys <- lapply(names(X), function(var) {
    x <- X[[var]]
    Pmat <- vapply(
      P_funs,
      FUN      = function(f) f(x),
      FUN.VALUE = numeric(length(x))
    )
    colnames(Pmat) <- sprintf("%s_P%d", var, idx_needed)
    as.data.frame(Pmat)
  })

  do.call(cbind, polys)
}


fourier_features <- function(X,
                             nbasis        = 5,         # total # basis fns in Fourier system
                             rangeval      = c(0, 3),   # domain of each variable
                             period        = diff(rangeval),
                             dropind       = NULL,      # e.g. dropind = 1L to omit constant
                             rescale_0_1   = FALSE,      # rescale each col to [0,1]?
                             ...) {                     # extra args passed to create.fourier.basis
  ## ---------------------------------------------------------------------------
  stopifnot(nbasis >= 1,
            is.data.frame(X) || is.matrix(X))
  X <- as.data.frame(X)

  ## -- optional rescaling so that rangeval stays [0,1] ------------------------
  if (isTRUE(rescale_0_1)) {
    X[] <- lapply(X, function(col) {
      rng <- range(col, na.rm = TRUE)
      if (diff(rng) == 0) rep(0, length(col)) else (col - rng[1]) / diff(rng)
    })
    rangeval <- c(0, 1)          # keep domain consistent after rescale
    period   <- diff(rangeval)
  }

  ## -- build a single Fourier basis object -----------------------------------
  fbasis <- fda::create.fourier.basis(rangeval = rangeval,
                                      nbasis   = nbasis,
                                      period   = period,
                                      dropind  = dropind,
                                      ...)     # axes, names, etc. :contentReference[oaicite:0]{index=0}

  ## -- evaluate that basis for EACH variable ----------------------------------
  mat_list <- lapply(names(X), function(var) {
    x <- X[[var]]
    if (any(x < rangeval[1] | x > rangeval[2])) {
      stop(sprintf("Values of '%s' lie outside 'rangeval'.", var))
    }
    B <- fda::eval.basis(x, fbasis)                      # n × nbasis'
    colnames(B) <- sprintf("%s_F%d", var, seq_len(ncol(B)) - 1L)
    as.data.frame(B)
  })

  do.call(cbind, mat_list)
}


spline_features <- function(X, df = 6, degree = 3, include_intercept = FALSE) {
  stopifnot(is.data.frame(X) || is.matrix(X))
  X <- as.data.frame(X)

  spline_list <- lapply(names(X), function(var) {
    x <- X[[var]]
    spline_mat <- splines::bs(
      x,
      df = df,
      degree = degree,
      intercept = include_intercept
    )
    colnames(spline_mat) <- sprintf("%s_S%d", var, seq_len(ncol(spline_mat)))
    as.data.frame(spline_mat)
  })

  do.call(cbind, spline_list)
}


`%||%` <- function(x, y) if (is.null(x)) y else x
