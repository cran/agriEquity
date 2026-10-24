.chk <- function(x, name, lower = -Inf, upper = Inf, scalar = TRUE) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be finite numeric with no missing values.", name),
         call. = FALSE)
  }
  if (scalar && length(x) != 1L) {
    stop(sprintf("`%s` must be a single number.", name), call. = FALSE)
  }
  if (any(x < lower) || any(x > upper)) {
    stop(sprintf("`%s` must lie between %s and %s.", name,
                 format(lower), format(upper)), call. = FALSE)
  }
  invisible(x)
}

.expand <- function(x, n, name) {
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) {
    stop(sprintf("`%s` must have length 1 or %d.", name, n), call. = FALSE)
  }
  x
}

.weights <- function(w, n) {
  if (is.null(w)) return(rep(1, n))
  .chk(w, "weights", lower = 0, scalar = FALSE)
  if (length(w) != n) stop("`weights` must have one value per unit.",
                           call. = FALSE)
  if (sum(w) <= 0) stop("`weights` must not all be zero.", call. = FALSE)
  w
}

.wmean <- function(x, w) sum(w * x) / sum(w)

# Weighted fractional rank (ties share their weighted average rank)
.frank <- function(x, w) {
  o <- order(x)
  cw <- cumsum(w[o])
  r <- numeric(length(x))
  r[o] <- (cw - w[o] / 2) / sum(w)
  key <- match(x, sort(unique(x)))
  sw <- tapply(w, key, sum)
  ra <- ifelse(sw > 0, tapply(r * w, key, sum) / sw, tapply(r, key, mean))
  as.numeric(ra[key])
}

# Concentration index of y ranked by `rank_by` (weighted covariance formula)
.conc <- function(y, rank_by, w) {
  m <- .wmean(y, w)
  if (m == 0) return(NA_real_)
  r <- .frank(rank_by, w)
  2 * sum(w * (y - m) * (r - .wmean(r, w))) / sum(w) / m
}

.gini <- function(y, w) .conc(y, y, w)

.fgt <- function(y, z, alpha, w) {
  g <- pmax((z - y) / z, 0)
  if (alpha == 0) return(.wmean(as.numeric(y < z), w))
  .wmean(g^alpha, w)
}

# Weighted quantile groups 1..k of x
.wq_groups <- function(x, w, k) {
  r <- .frank(x, w)
  pmin(k, floor(r * k) + 1L)
}

.fmt <- function(x, digits = 3) {
  out <- formatC(x, format = "f", digits = digits, big.mark = ",")
  out[is.na(x)] <- "NA"
  out
}
.pct <- function(x, digits = 1) {
  out <- paste0(.fmt(100 * x, digits), "%")
  out[is.na(x)] <- "NA"
  out
}
