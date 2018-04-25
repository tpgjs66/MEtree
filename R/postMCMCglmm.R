paramNamesMCMCglmm <- function(object, ...) {
  fNames <- as.character(attr(terms(object$Fixed$formula), "variables"))[-c(1, 2)]

  rNames <- NULL

  if (!is.null(object$Random$formula)) {
    rNames <- as.character(attr(terms(object$Random$formula), "variables"))[-c(1)]
  }

  if ("(Intercept)" %in% colnames(object[["Sol"]])) {
    fNames <- c("(Intercept)", fNames)
  }

  list(fixed = fNames, random = rNames)
}

ranefLevels <- function(object, data, ...) {
  n <- paramNamesMCMCglmm(object)$random
  res <- lapply(n, function(n) {
    levels(data[, n])
  })
  names(res) <- n
  return(res)
}

.extractEffects <- function(object, use = c("all", "mean"),
                            which = c("fixed", "random"), ...) {

  use <- match.arg(use)
  which <- match.arg(which)

  b <- as.matrix(object[["Sol"]])

  eff <- switch(which,
                fixed = {
                  # this does not work because the intercept term contains
                  # special characters (), that screw with the regular expressions
                  # cannot use fixed matching because factor variables can
                  # expand with arbitrarily named levels
                  #unlist(lapply(paramNamesMCMCglmm(object)$fixed, function(n) {
                  #  grep(paste0("^", n, ".*$"), colnames(b), value = TRUE)
                  #}))
                  object$X@Dimnames[[2]]
                },
                random = {
                  # regex <- paste(paramNamesMCMCglmm(object)$random, collapse = "|")
                  # regex <- paste0("^(", regex, ")\\..*$")

                  regex <- paramNamesMCMCglmm(object)$random[2]
                  grep(regex, colnames(b), value = TRUE)
                }
  )

  b <- b[, eff, drop=FALSE]

  switch(use,
         all = t(b),
         mean = as.matrix(colMeans(b)))
}


fixef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "fixed", ...)
}

ranef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "random", ...)
}

stdranef <- function(object, which, type = c("lp", "response"), ...) {
  type <- match.arg(type)

  if (is.null(object$Z)) stop("Z matrix must be saved")
  ## z <- object$Z
  z <- diag(ncol(object$Z))
  colnames(z) <- colnames(object$Z)

  re <- paramNamesMCMCglmm(object)$random

  if (missing(which)) which <- c(re, list(re))

  stopifnot(is.list(which))

  if (is.numeric(unlist(which))) {
    stopifnot(all(unlist(which) %in% seq_along(re)))
    which <- lapply(which, function(i) re[i])
  } else {
    stopifnot(all(unlist(which) %in% re))
  }

  index <- lapply(which, function(n) {
    n <- paste(n, collapse = "|")
    regex <- paste0("^(", n, ")\\..*$")
    index <- grep(regex, colnames(z))
    return(index)
  })

  # tricky, note that coefficients and predicted probabilities
  # do not come out with the same dimensions, the are transposed
  # samples on columns for coefs, samples on rows for predicted probs
  res <- switch(type,
                lp = {
                  yhat <- lapply(index, function(i) ranef(object, use = "all")[i, ])
                  lapply(yhat, function(m) matrix(apply(m, 2, sd)))
                }, response = {
                  yhat <- lapply(index, function(i) {
                    tmp <- z
                    if (length(i) < ncol(tmp)) {
                      tmp[, -i] <- 0L # zero out all random effects we are not interested in
                    }
                    predict2(object, X = NULL, Z = tmp, use = "all", type = type)
                  })
                  lapply(yhat, function(m) sapply(m, function(n) apply(n, 1, sd)))
                })

  names(res) <- which

  M <- do.call(rbind, lapply(res, colMeans))

  finalres <- list(M = M, Data = res)
  class(finalres) <- "postMCMCglmmRE"

  return(finalres)
}
