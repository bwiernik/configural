#' Locate extrema of fungible OLS regression weights
#'
#' @param theta The value of the R-squared decrement used to generate a family of fungible coefficients.
#' @param Rxx An intercorrelation matrix among the predictor variables
#' @param rxy A vector of predictor–criterion correlations
#' @param Nstarts The maximum number
#' @param MaxMin Should the cosine between the OLS and alternative weights be maximized ("max") to find the maximally similar coefficients or minimized ("min") to find the maximally dissimilar coefficients?
#' @param silent Should current optimization values be printed to the console (`FALSE`) or suppressed (`TRUE`)?
#'
#' @return A list containing the alternative weights and other fungible weights estimation parameters
#' @author Adapted from [fungible::fungibleExtrema()] by Niels Waller
#'
#' @export
#'
#' @keywords internal
.fungible_extrema <-
  function(theta, Rxx, rxy, Nstarts = 1000, MaxMin = c("min", "max"), silent = FALSE) {
    MaxMin <- match.arg(MaxMin)

    norm <- function(x) x / as.numeric(sqrt(t(x) %*% x))
    vcos <- function(x, y) t(norm(x)) %*% norm(y)
    GenU <- function(mat, u) qr.Q(qr(cbind(u, mat)))[, -1]

    gradf <- function(sv) {
      z <- sv[1:(p - 1)]
      zpz <- t(z) %*% z
      L <- sv[p]
      k <- r * u + U %*% z * e
      a <- V %*% Linv.sqrt %*% k
      apa <- t(a) %*% a
      ab <- as.numeric(t(a) %*% b.unit)
      dfdz <-
        (e * t(U) %*% Linv.sqrt %*% t(V)) %*%
        (-as.numeric(ab * (apa)^-1.5) * a + as.numeric((apa^-0.5)) * b.unit)
      if (MaxMin == "min") {
        dfdz <- dfdz + 4 * L * as.numeric(zpz - 1) * z
        dfdL <- (zpz - 1)^2
      } else {
        dfdz <- dfdz - 4 * L * as.numeric(zpz - 1) * z
        dfdL <- -(zpz - 1)^2
      }
      c(dfdz, dfdL)
    }

    minv <- function(sv) {
      z <- sv[1:(p - 1)]
      zpz <- t(z) %*% z
      L <- sv[p]
      k <- r * u + U %*% z * sqrt(1 - r^2)
      k.star <- Linv.sqrt %*% k
      a <- V %*% k.star
      len.a <- sqrt(t(a) %*% a)
      f <- (as.numeric(len.a)^-1) * t(a) %*% b.unit + L * (zpz - 1)^2
      if (MaxMin == "max") {
        f <- (as.numeric(len.a)^-1) * t(a) %*% b.unit - L * (zpz - 1)^2
      } else {
        f <- (as.numeric(len.a)^-1) * t(a) %*% b.unit + L * (zpz - 1)^2
      }
      return(f)
    }

    p <- ncol(Rxx)
    b <- crossprod(solve(Rxx), rxy)
    b.unit <- norm(b)
    R2 <- as.numeric(crossprod(rxy, b))
    r <- sqrt(1 - theta / R2)
    e <- sqrt(1 - r^2)

    VLV <- eigen(Rxx)
    V <- VLV$vectors
    L <- diag(VLV$values)
    Linv.sqrt <- solve(sqrt(L))
    u <- (sqrt(L) %*% t(V) %*% b) / as.numeric(sqrt(t(b) %*% Rxx %*% b))
    mat <- matrix(stats::rnorm(p * (p - 1)), p, p - 1)
    U <- GenU(mat, u)

    FNSCALE <- 1
    if (MaxMin == "max") FNSCALE <- FNSCALE * -1
    minf <- 99 * FNSCALE
    breakflag <- 0
    iter <- 0
    solutions <- 0
    while (iter < Nstarts) {
      sv <- c(stats::rnorm(p - 1) / sqrt(p - 1), 1000)
      tmp <- try(suppressWarnings(
        optimx::optimr(par = sv, fn = minv, gr = gradf, method = "BFGS",
                       control = list(fnscale = FNSCALE, maxit = 500,
                                      parscale = c(rep(1, p - 1), 1)))))
      if (abs(tmp$value) > 1) tmp$convergence <- 1
      if (tmp$convergence == 0) {
        iter <- iter + 1
        if ((FNSCALE * tmp$value <= FNSCALE * minf)) {
          solutions <- solutions + 1
          minf <- tmp$value
          out <- tmp
          z <- out$par[1:(p - 1)]
          k <- r * u + U %*% z * sqrt(1 - r^2)
          a <- a.tilde <- V %*% Linv.sqrt %*% k
          scaling.weight <- (t(rxy) %*% a) / (t(a) %*% Rxx %*% a)
          a <- as.numeric(scaling.weight) * a
          if (max(abs(gradf(out$par))) <= 1e-05) breakflag <- 1
        }
        if (!silent & iter %% 100 == 0) {
          cat("Solution count:", iter,
              "  Current function val:", minf,
              "  Gradient:", max(abs(gradf(out$par))), "\n")
        }
      }
      if (breakflag) break
    }

    if (sign(a[1]) != sign(a.tilde[1])) a <- a * -1

    R.b <- sqrt(R2)
    R2.b <- R2
    R.a <- R.b * r
    R2.a <- R.a^2
    s <- vcos(a, b)
    solution.gradient <- gradf(out$par)
    names(a) <- names(b) <- colnames(Rxx)

    out <-
      list(a = a, b = b, cos.ab = s, k = k, z = z, u = u, r.yhata.yhatb = r,
           theta = theta, multiple.r.b = R.b, r.squared.b = R2.b,
           multiple.r.a = R.a, r.squared.a = R2.a,
           gradient = solution.gradient, MaxMin = MaxMin)

    return(out)
  }

#' Locate extrema of fungible weights for regression and related models
#'
#' Generates fungible regression weights (Waller, 2008) and related results using the method by Waller and Jones (2010).
#'
#' @param object A fitted model object. Currently supported classes are: "cpa"
#' @param theta A vector of values to decrement from R-squared to compute families of fungible coefficients.
#' @param Nstarts Maximum number of (max) minimizations from random starting configurations.
#' @param MaxMin Should the cosine between the observed and alternative weights be maximized ("max") to find the maximally similar coefficients or minimized ("min") to find the maximally dissimilar coefficients?
#' @param silent Should current optimization values be printed to the console (`FALSE`) or suppressed (`TRUE`)?
#' @param ... Additional arguments
#'
#' @return A list containing the alternative weights and other fungible weights estimation parameters
#'
#' @author Niels Waller, Jeff Jones, Brenton M. Wiernik. Adapted from [fungible::fungibleExtrema()].
#'
#' @export
#'
#' @references
#' Waller, N. G. (2008).
#' Fungible weights in multiple regression.
#' _Psychometrika, 73_(4), 691–703. <https://doi.org/10/c6h5qn>
#'
#' Waller, N. G., & Jones, J. A. (2009).
#' Locating the extrema of fungible regression weights.
#' _Psychometrika, 74_(4), 589–602. <https://doi.org/10/c3wbtd>
#'
#' @examples
#' \dontrun{
#'   mind <- cpa_mat(mindfulness ~ ES + A + C + Ex + O,
#'                   cov_mat = mindful_rho,
#'                   n = 17060)
#'   mind_fung <- fungible(mind)
#' }
fungible <- function(object, theta = .005, Nstarts = 1000,
                     MaxMin = c("min", "max"), silent = FALSE, ...) {
  UseMethod("fungible")
}

#' Locate extrema of fungible criterion profile patterns
#'
#' Identify maximally similar or dissimilar criterion patterns in criterion profile analysis
#'
#' @param object A fitted model object of class "cpa".
#' @param theta A vector of values to decrement from R-squared to compute families of fungible coefficients.
#' @param Nstarts Maximum number of (max) minimizations from random starting configurations.
#' @param MaxMin Should the cosine between the observed and alternative weights be maximized ("max") to find the maximally similar coefficients or minimized ("min") to find the maximally dissimilar coefficients?
#' @param silent Should current optimization values be printed to the console (`FALSE`) or suppressed (`TRUE`)?
#' @param ... Additional arguments
#'
#' @references
#' Wiernik, B. M., Wilmot, M. P., Davison, M. L., & Ones, D. S. (2019).
#' _Meta-analytic criterion profile analysis_.
#' Manuscript submitted for publication, University of South Florida.
#'
#' @return A list containing the alternative weights and other fungible weights estimation parameters
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   mind <- cpa_mat(mindfulness ~ ES + A + C + Ex + O,
#'                   cov_mat = mindful_rho,
#'                   n = 17060)
#'   mind_fung <- fungible(mind)
#' }
fungible.cpa <- function(object, theta = .005, Nstarts = 1000,
                         MaxMin = c("min", "max"), silent = FALSE, ...) {
  if (!inherits(object, "cpa")) stop("'object' must have class 'cpa'")
  MaxMin <- match.arg(MaxMin)

  R2 <- object$fit["Total", "R-squared"]
  lev2 <- object$fit["Level", "R-squared"]
  Rxx <- object$model$Rxx
  rxy <- object$model$rxy
  theta <- theta[theta > 0 & theta < R2]

  fun_ex <- vector("list", length(theta))
  names(fun_ex) <- paste("theta =", theta)

  for (i in 1:length(theta)) {
    fun_ex[[i]] <- .fungible_extrema(theta[[i]], Rxx = Rxx, rxy = rxy,
                                     Nstarts = Nstarts, MaxMin = MaxMin, silent = silent)
  }
  for (i in 1:length(theta)) {
    fun_ex[[i]]$astar <- fun_ex[[i]]$a - mean(fun_ex[[i]]$a)
    attributes(fun_ex[[i]]$astar) <- attributes(fun_ex[[i]]$a)
    fun_ex[[i]]$pat2_a <- (t(rxy) %*% fun_ex[[i]]$astar)^2 / (t(fun_ex[[i]]$astar) %*% Rxx %*% fun_ex[[i]]$astar)
    fun_ex[[i]]$pat_a <- sqrt(fun_ex[[i]]$pat2_a)
    fun_ex[[i]]$delta_lev_a <- fun_ex[[i]]$r.squared.a - fun_ex[[i]]$pat2_a
    fun_ex[[i]]$delta_pat_a <- fun_ex[[i]]$r.squared.a - lev2
  }

  b <- matrix(as.numeric(object$beta[,1]), 1)
  rownames(b) <- "b (OLS)"
  bstar <- matrix(as.numeric(object$bstar[,1]), 1)
  rownames(bstar) <- "bstar"
  fit_bstar <- object$fit[1:2, c(1, 5, 10, 10)]
  rownames(fit_bstar) <- c("Level", "Pattern (bstar)")
  colnames(fit_bstar) <- c("R", "R-squared", "Delta R-squared (Level)", "Delta R-squared (Pattern)")
  fit_bstar[2, 4] <- fit_bstar[1, 4]
  fit_bstar[1, 3:4] <- NA_integer_

  a <- t(sapply(fun_ex, function(x) x$a))
  astar <- t(sapply(fun_ex, function(x) x$astar))
  fit_astar <- t(sapply(fun_ex, function(x) x[c("pat_a", "pat2_a", "delta_lev_a", "delta_pat_a")]))
  fungible_params <- t(sapply(fun_ex, function(x) x[c("theta", "r.yhata.yhatb", "cos.ab", "k", "z", "u")]))
  gradient <- t(sapply(fun_ex, function(x) x$gradient))

  rownames(a) <- paste0("a (theta = ", theta, ")")
  rownames(astar) <- paste0("astar (theta = ", theta, ")")
  rownames(fit_astar) <- paste0("Pattern (astar, theta = ", theta, ")")
  rownames(fungible_params) <- rownames(gradient) <- paste0("theta = ", theta)

  coefficients_ols <- rbind(b, a)
  coefficients_cpa <- rbind(bstar, astar)
  fit_cpa <- rbind(fit_bstar, fit_astar)
  dim_names <- dimnames(fit_cpa)
  fit_cpa <- apply(fit_cpa, 2, as.numeric)
  dimnames(fit_cpa) <- dim_names

  out <- list(coefficients_cpa = coefficients_cpa, coefficients_ols = coefficients_ols,
              fit_cpa = fit_cpa, fungible_params = fungible_params,
              gradient = gradient, MaxMin = MaxMin, vcov = object$vcov)
  class(out) <- "fungible_extrema"
  return(out)

}

print.fungible_extrema <- function(object,
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  cat("\nFungible weights analysis\n")
  cat(  "=========================\n\n")
  if (object$MaxMin == "max") {
    cat("  Identifying most similar coefficients\n\n")
  } else {
    cat("  Identifying most dissimilar coefficients\n\n")
  }

  if (exists("coefficients_ols", object)) {
    cat("OLS Regression:\n\n")

    cat("Regression coefficients:\n")
    print(round(object$coefficients_ols, digits))
    cat("\n\n")
  }

  if (exists("coefficients_cpa", object)) {
    cat("Criterion Profile Analysis:\n\n")

    cat("Criterion pattern:\n")
    print(round(object$coefficients_cpa, digits))
    cat("\n\n")

    cat("Model comparison:\n")
    print(round(object$fit_cpa, digits))
    cat("\n\n")
  }

  invisible(object)
}

#' Locate extrema of fungible OLS regression weights
#'
#' Identify maximally similar or dissimilar sets of fungible standardized regression coefficients from an OLS regression model
#'
#' @param object A fitted model object of class "lm" or "summary.lm".
#' @param theta A vector of values to decrement from R-squared to compute families of fungible coefficients.
#' @param Nstarts Maximum number of (max) minimizations from random starting configurations.
#' @param MaxMin Should the cosine between the observed and alternative weights be maximized ("max") to find the maximally similar coefficients or minimized ("min") to find the maximally dissimilar coefficients?
#' @param silent Should current optimization values be printed to the console (`FALSE`) or suppressed (`TRUE`)?
#' @param ... Additional arguments
#'
#' @references
#' Waller, N. G., & Jones, J. A. (2009).
#' Locating the extrema of fungible regression weights.
#' _Psychometrika, 74_(4), 589–602. <https://doi.org/10/c3wbtd>
#'
#' @return A list containing the alternative weights and other fungible weights estimation parameters
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   lm_mtcars <- lm(mpg ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb,
#'                   data = mtcars)
#'   lm_mtcars_fung <- fungible(lm_mtcars)
#' }
fungible.lm <- function(object, theta = .005, Nstarts = 1000,
                        MaxMin = c("min", "max"), silent = FALSE, ...) {
  if (!inherits(object, "lm")) stop("'object' must have class 'lm'")
  if (!is.null(stats::model.offset(object))) stop("models with offsets not yet supported")
  MaxMin <- match.arg(MaxMin)

  X <- stats::model.matrix(object)
  if ("(Intercept)" %in% colnames(X)) X <- X[, -1]
  y <- object$model[, 1]
  w <- stats::model.weights(object)
  if (is.null(w)) w <- rep(1, nrow(X))

  if (is.null(stats::na.action(object))) {
    corr <- wt_cor(cbind(X, y), wt = w, use = "listwise")
    Rxx <- corr[1:ncol(X), 1:ncol(X)]
    rxy <- corr[1:ncol(X), ncol(X) + 1]
  }
  else if (stats::na.action(object) %in% c("na.omit", "na.exclude")) {
    corr <- wt_cor(cbind(X, y), wt = w, use = "listwise")
    Rxx <- corr[1:ncol(X), 1:ncol(X)]
    rxy <- corr[1:ncol(X), ncol(X) + 1]
  } else {
    corr <- wt_cor(cbind(X, y), wt = w)
    Rxx <- corr[1:ncol(X), 1:ncol(X)]
    rxy <- corr[1:ncol(X), ncol(X) + 1]
  }

  R2 <- t(rxy) %*% solve(Rxx) %*% rxy

  theta <- theta[theta > 0 & theta < R2]

  fun_ex <- vector("list", length(theta))
  names(fun_ex) <- paste("theta =", theta)

  for (i in 1:length(theta)) {
    fun_ex[[i]] <- .fungible_extrema(theta[[i]], Rxx = Rxx, rxy = rxy,
                                     Nstarts = Nstarts, MaxMin = MaxMin, silent = silent)
  }

  b <- t(rxy) %*% solve(Rxx)
  rownames(b) <- "b (OLS)"

  a <- t(sapply(fun_ex, function(x) x$a))
  fungible_params <- t(sapply(fun_ex, function(x) x[c("theta", "r.yhata.yhatb", "cos.ab", "k", "z", "u")]))
  gradient <- t(sapply(fun_ex, function(x) x$gradient))

  rownames(a) <- paste0("a (theta = ", theta, ")")
  rownames(fungible_params) <- rownames(gradient) <- paste0("theta = ", theta)

  coefficients_ols <- rbind(b, a)

  out <- list(coefficients_ols = coefficients_ols,
              fungible_params = fungible_params,
              gradient = gradient, MaxMin = MaxMin)
  class(out) <- "fungible_extrema"
  return(out)

}
