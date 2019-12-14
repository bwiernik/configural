#' Estimate the sampling error variance for criterion profile analysis parameters
#'
#' @param Rxx An intercorrelation matrix among the predictor variables
#' @param rxy A vector of predictor–criterion correlations
#' @param n The sample size
#' @param se_var_mat The method used to calculate the sampling covariance matrix for the unique elements of Rxx and rxy. Can be "normal" to estimate the covariance matrix using "n" or a matrix of sampling covariance values.
#' @param adjust Method to adjust R-squared for overfitting. See \code{\link{adjust_Rsq}} for details.
#'
#' @return A list containing sampling covariance matrices or sampling erorr variance estimates for CPA parameters
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' var_error_cpa(mindfulness$rho[1:5, 1:5], mindfulness$rho[1:5, 6], n = 17060)
var_error_cpa <- function(Rxx, rxy, n, se_var_mat = "normal", adjust = c("fisher", "pop", "cv")) {
  adjust <- match.arg(adjust)
  Rxx <- as.matrix(Rxx)
  p <- ncol(Rxx)
  Rinv <- solve(Rxx)
  Id <- diag(p)
  one <- rep(1, p)
  Q <- Id - one %*% t(one) / p

  if (p == 1) {
    v.beta <- v.R <- ((1 - rxy^2)^2) / (n - 3)
    v.R2 <- 4 * rxy^2 * v.R
    v.bstar <- v.lev2 <- v.lev <- v.pat2 <- v.pat <-
      v.levpat2 <- v.levpat <- v.deltalev <- v.deltapat <- NULL
  } else {
    sR <- rbind(cbind(Rxx, rxy), c(rxy, 1))

    if (inherits(se_var_mat, "matrix") ) {
      Sigma <- se_var_mat
    } else {
      Sigma <- cor_covariance(sR, n)
    }

    Kpc <- as.matrix(transition(p))
    if (ncol(Kpc) == 1) Kpc <- t(Kpc)

    # Reorder matrix
    rxx.nms <- matrix(0, p, p)
    rxy.nms <- c(rep(0, p + 1))
    for (i in 1:p) for (j in 1:p) rxx.nms[i, j] <- paste("rx", i, "rx", j, sep = "")
    for (i in 1:p + 1) rxy.nms[i] <- paste("rx", i, "y", sep = "")
    nm.mat <- rbind(cbind(rxx.nms, rxy.nms[-(p + 1)]), rxy.nms)
    old.ord <- nm.mat[lower.tri(nm.mat)]
    new.ord <- c(rxx.nms[lower.tri(rxx.nms)], rxy.nms)

    # Beta
    beta <- Rinv %*% rxy
    dbeta.drxx <- -2 * (t(beta) %x% Rinv) %*% t(Kpc)
    dbeta.drxy <- Rinv
    j.beta <- cbind(dbeta.drxx, dbeta.drxy)
    j.beta <- j.beta[, match(old.ord, new.ord)]
    v.beta <- j.beta %*% Sigma %*% t(j.beta)

    # BStar
    bstar <- Q %*% beta
    j.bstar <- Q %*% j.beta
    v.bstar <- j.bstar %*% Sigma %*% t(j.bstar)

    # R Squared
    R2 <- t(rxy) %*% beta
    dR2.drxx <- t(rxy) %*% dbeta.drxx
    dR2.drxy <- 2 * t(rxy) %*% Rinv
    j.R2 <- cbind(dR2.drxx, dR2.drxy)
    j.R2 <- j.R2[, match(old.ord, new.ord)]
    v.R2 <- j.R2 %*% Sigma %*% j.R2
    v.R <- v.R2 / (4 * R2)

    # Level Squared
    lev <- (t(rxy) %*% one) / sqrt(t(one) %*% Rxx %*% one)
    lev2 <- lev^2
    dlev2.drxx <- -1 * ((t(rxy) %*% one) / (t(one) %*% Rxx %*% one))^2 %*% (t(one) %x% t(one)) %*% t(Kpc)
    dlev2.drxy <- 2 * ((t(rxy) %*% one) / (t(one) %*% Rxx %*% one)) %*% t(one)
    j.lev2 <- cbind(dlev2.drxx, dlev2.drxy)
    j.lev2 <- j.lev2[, match(old.ord, new.ord)]
    v.lev2 <- j.lev2 %*% Sigma %*% j.lev2
    v.lev <- v.lev2 / (4 * lev2)

    # Pattern Squared
    pat <- (t(rxy) %*% bstar) / sqrt(t(bstar) %*% Rxx %*% bstar)
    pat2 <- pat^2
    dpat2.drxx <-
      -1 * ((t(rxy) %*% bstar) / (t(bstar) %*% Rxx %*% bstar)) %*%
      ((t(beta) %x% (t(rxy) %*% Q %*% Rinv)) + ((t(rxy) %*% Q %*% Rinv) %x% t(beta)) -
         ((t(rxy) %*% bstar) / (t(bstar) %*% Rxx %*% bstar)) %*%
         ((t(bstar) %x% t(bstar)) - (t(beta) %x% (t(bstar) %*% Rxx %*% Q %*% Rinv)))) %*% t(Kpc)
    dpat2.drxy <-
      2 * ((t(rxy) %*% bstar) / (t(bstar) %*% Rxx %*% bstar)) %*%
      (t(rxy) %*% Q %*% Rinv + t(bstar) -
         (solve((t(bstar) %*% Rxx %*% bstar)) %*% (t(rxy) %*% bstar %*% t(bstar) %*% Rxx %*% Q %*% Rinv) ))
    j.pat2 <- cbind(dpat2.drxx, dpat2.drxy)
    j.pat2 <- j.pat2[, match(old.ord, new.ord)]
    v.pat2 <- j.pat2 %*% Sigma %*% j.pat2
    v.pat <- v.pat2 / (4 * pat2)

    # Level-Pattern Squared
    levpat <- (t(one) %*% Rxx %*% bstar) / sqrt(t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one)
    levpat2 <- levpat^2
    dlevpat2.drxx <-
      (1 / (t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one)) %*%
      (((t(bstar) %x% t(one)) - (t(beta) %x% (t(one) %*% Rxx %*% Q  %*% Rinv))) -
         (((t(one) %*% Rxx %*% one %*% t(one) %*% Rxx %*% bstar) /
             (t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one)) %*%
            ((t(bstar) %x% t(bstar)) - (t(beta) %x% (t(bstar) %*% Rxx %*% Q  %*% Rinv)) -
               ((t(bstar) %*% Rxx %*% Q %*% Rinv) %x% t(bstar))) +
            (((t(one) %*% Rxx %*% bstar) / (t(one) %*% Rxx %*% one) ) %*% (t(one) %x% t(one))) )) %*% t(Kpc)
    dlevpat2.drxy <-
      (solve(t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one) %*% (t(one) %*% Rxx %*% Q %*% Rinv)) -
        2 * (solve((t(bstar) %*% Rxx %*% bstar %*% t(one) %*% Rxx %*% one)^2) %*%
               (t(one) %*% Rxx %*% one %*% t(bstar) %*% Rxx %*% one)) %*%
            (t(bstar) %*% Rxx %*% Q %*% Rinv)
    j.levpat2 <- cbind(dlevpat2.drxx, dlevpat2.drxy)
    j.levpat2 <- j.levpat2[, match(old.ord, new.ord)]
    v.levpat2 <- j.levpat2 %*% Sigma %*% j.levpat2
    v.levpat <- v.levpat2 / (4 * levpat2)

    # Delta_Level
    v.deltalev <-
      4 * R2 * v.R + 4 * lev2 * v.lev -
      8 * (sqrt(lev2) * (R2^3 + R2^2 * (lev2 - 3) +
                           R2 * (2 - lev2) - lev2) /
             (2 * R2^1.5 * (1 - R2) * (1 - lev2))) *
      sqrt(v.R * v.lev * R2 * lev2)

    # Delta_Pattern
    v.deltapat <-
      4 * R2 * v.R + 4 * pat2 * v.pat -
      8 * (sqrt(pat2) * (R2^3 + R2^2 * (pat2 - 3) +
                           R2 * (2 - pat2) - pat2) /
             (2 * R2^1.5 * (1 - R2) * (1 - pat2))) *
      sqrt(v.R * v.pat * R2 * pat2)

    # Adjusted R Squared
    R2_adj <- adjust_Rsq(R2, n, p, adjust)
    if (adjust == "fisher") {
      d2adj.R2 <- ((n - 1) / (n - p - 1))^2
    } else if (adjust == "pop") {
      d2adj.R2 <- ((n^2 + n * (-1 * p - 4 * R2 - 1.3) + 3 * p + 12 * R2 - 5.1)/((n - p - 2.3) * (n - p - 1)))^2
    } else {
      d2adj.R2 <- (
        ((2 - n + 2 * p) *
           ((n - 1) * (n - p - 2.3)^2 * (n - p - 1)^3 -
              (3 - n + p) *
              ((1 - n) * (1 - n + p) *
                 ((1 - n + p) * (2.3 - n + p) +
                    (n - 3) * (n - .3 - p - 2 * R2) * (R2 - 1))^2 -
                 2 * (n - 3)^2 * p * (n - .3 - p - 2 * R2)^2 * (R2 - 1)^2) +
              (n - 3) * (n - 1) *
              (n - p - 2.3) * (n - p - 1)^2 *
              (n - .3 - p - 2 * R2) * (R2 - 1)) *
           (n^2 - 5.1 + 3 * p + n * (-1.3 - p - 4 * R2) + 12 * R2) +
           (p * (1 - n + p) * (2.3 - n + p) -
              (2 - n + 2 * p) * ((1 - n + p) * (2.3 - n + p) +
                                   (n - 3) * (n - .3 - p - 2 * R2) * (R2 - 1))) *
           ((n - 3) * (n - 1) * (n - 2.3 - p) * (n - 1 - p)^2 *
              (n - .3 - p - 2 * R2) - 2 * (n - 3) * (n - 1) * (n - 2.3 - p) *
              (n - 1 - p)^2 * (R2 - 1) -
              2 * (3 - n + p) *
              (-2 * (n - 3)^2 * p * (n - .3 - p - 2 * R2)^2 *
                 (R2 - 1) + 4 * (n - 3)^2 * p * (n - .3 - p - 2 * R2) *
                 (R2 - 1)^2 + (n - 1) * (n - 1 - p) *
                 (n^2 - 5.1 + 3 * p + n * (-1.3 - p - 4 * R2) + 12 * R2) *
                 (1.4 + p^2 + (n^2 - 1.3 * n - 5.1) * R2 +
                    (6 - 2 * n) * R2^2 + p * (.3 + n * (-1 - R2) + 3 * R2))))) /
          ((1 - n) * (1 - n + p)^4 * (2.3 - n + p)^3 *
             (p + (n - 2 * (1 + p)) *
                (1 + ((n - 3) * (n - .3 - p - 2 * R2) * (R2 - 1)) /
                   ((n - 2.3 - p) * (n - p - 1))))^2)
      )^2
    }
    v.R2_adj <- d2adj.R2 * v.R2
    v.R_adj <- v.R2_adj / (4 * R2_adj)

    # Adjusted Pattern Squared
    pat_adj <- sqrt(lev2) * sqrt(levpat2) + sqrt(max((1 - levpat2) * (R2_adj - lev2), 0))
    pat2_adj <- pat_adj^2
    # TODO: Replace this with a proper delta method se_var (as R2_adj above)
    v.pat_adj <- v.pat / (pat2 / pat2_adj)
    v.pat2_adj <- v.pat2 / (pat2 / pat2_adj)^2

    # Adjusted Delta_Level
    v.deltalev_adj <-
      4 * R2_adj * v.R_adj + 4 * lev2 * v.lev -
      8 * (sqrt(lev2) * (R2_adj^3 + R2_adj^2 * (lev2 - 3) +
                           R2_adj * (2 - lev2) - lev2) /
             (2 * R2^1.5 * (1 - R2_adj) * (1 - lev2))) *
      sqrt(v.R_adj * v.lev * R2_adj * lev2)

    # Adjusted Delta_Pattern
    v.deltapat_adj <-
      4 * R2_adj * v.R_adj + 4 * pat2_adj * v.pat_adj -
      8 * (sqrt(pat2_adj) * (R2_adj^3 + R2_adj^2 * (pat2_adj - 3) +
                               R2_adj * (2 - pat2_adj) - pat2_adj) /
             (2 * R2_adj^1.5 * (1 - R2_adj) * (1 - pat2_adj))) *
      sqrt(v.R_adj * v.pat_adj * R2_adj * pat2_adj)

    if (is.null(colnames(Rxx))) var.names <- 1:p else var.names <- colnames(Rxx)
    rownames(v.beta) <- colnames(v.beta) <- paste0("beta_", var.names)
    rownames(v.bstar) <- colnames(v.bstar) <- paste0("bstar_", var.names)
  }

  out <- list(
    beta = v.beta,
    bstar = v.bstar,
    r.total = v.R,
    r.total.squared = v.R2,
    r.level = v.lev,
    r.level.squared = v.lev2,
    r.pattern = v.pat,
    r.pattern.squared = v.pat2,
    r.level.pattern = v.levpat,
    r.level.pattern.squared = v.levpat2,
    delta.r.squared.level = v.deltalev,
    delta.r.squared.pattern = v.deltapat,
    adjusted.r.total = v.R_adj,
    adjusted.r.total.squared = v.R2_adj,
    adjusted.r.pattern = v.pat_adj,
    adjusted.r.pattern.squared = v.pat2_adj,
    adjusted.delta.r.squared.level = v.deltalev_adj,
    adjusted.delta.r.squared.pattern = v.deltapat_adj
  )

  class(out) <- "var_cpa"

  return(out)
}
