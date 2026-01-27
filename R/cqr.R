wq <- function(y,w,prob){
  o <- order(y); y <- y[o]; w <- w[o]
  a <- 0; i <- 0; s <- sum(w); res <- c()
  for(q in prob){
    while(a < q){ i <- i + 1; a <- a + w[i]/s }
    res <- c(res, y[i])
  }
  res
}

est_cqr <- function(y, c, x, w, nqreg, c1, c2, nsteps){
  x    <- cbind(1, x)
  p    <- ncol(x)
  coef <- matrix(, p, nqreg)
  ncensored <- (y > c)
  pred <- fitted(glm(ncensored ~ x - 1, weights = w, family = binomial(logit)))
  coefqr <- est_qr(y, x, w, nqreg)

  for(i in 1:nqreg){
    quant <- i/nqreg - 0.5/nqreg
    fit   <- x %*% coefqr[, i]
    if(sum(fit < c) > 0){
      delta1 <- wq(pred[pred > (1 - quant)], w[pred > (1 - quant)], c1)
      temps  <- lm(y[pred >= delta1] ~ x[pred >= delta1,] - 1, weights = w[pred >= delta1])$coef
      reg1   <- x[, (1:p)[is.na(temps) == FALSE]]
      temp   <- temps
      temp[(1:p)[is.na(temps) == FALSE]] <- quantreg(y[pred >= delta1], reg1[pred >= delta1,], w[pred >= delta1], quant)
      temp[(1:p)[is.na(temps) == TRUE]]  <- (-100)
      step <- 3
      while(step <= nsteps){
        pred1 <- x %*% temp
        delta2 <- wq(pred1[pred1 > c], w[pred1 > c], c2)
        temps  <- lm(y[pred1 >= delta2] ~ x[pred1 >= delta2,] - 1, weights = w[pred1 >= delta2])$coef
        reg1   <- x[, (1:p)[is.na(temps) == FALSE]]
        temp   <- temps
        temp[(1:p)[is.na(temps) == FALSE]] <- quantreg(y[pred1 >= delta2], reg1[pred1 >= delta2,], w[pred1 >= delta2], quant)
        temp[(1:p)[is.na(temps) == TRUE]]  <- (-100)
        step <- step + 1
      }
      coef[, i] <- temp
    } else {
      coef[, i] <- coefqr[, i]
    }
  }

  # --- only difference: wrap columns as "qr_like" objects with predict() support
  feature_names <- colnames(x) %||% paste0("V", seq_len(ncol(x)))  # ensure names
  fits <- lapply(seq_len(nqreg), function(i){
    structure(list(
      coef = coef[, i],
      feature_names = feature_names
    ), class = "qr_like")
  })
  names(fits) <- paste0("tau", (1:nqreg)/nqreg - 0.5/nqreg)
  #names(fits) <- paste0("tau", tau_grid)
  # print(fits)
  list(fits = fits, X_poly = x)
}
