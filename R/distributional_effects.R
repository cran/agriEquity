#' Distributional effects of subsidies on income, inequality and poverty
#'
#' Compares farm household income without support with income including
#' the benefits under one or more subsidy scenarios. For each scenario it
#' reports mean income, the Gini coefficient, the redistributive effect,
#' the Foster-Greer-Thorbecke poverty measures \eqn{FGT_0} (headcount),
#' \eqn{FGT_1} (poverty gap) and \eqn{FGT_2} (severity) at the given poverty
#' line, the share of households lifted out of poverty, and the cost per
#' household lifted. A quantile table shows mean benefit and the percentage
#' income gain in each quantile of pre-subsidy income.
#'
#' @param income Numeric vector of pre-subsidy household income.
#' @param benefits A numeric vector of benefits, or a named list (or data
#'   frame) of benefit vectors, one per scenario.
#' @param weights Optional survey weights.
#' @param poverty_line Poverty line in the same units as income.
#' @param n_quantiles Number of income quantiles in the quantile table.
#'
#' @return A list of class `distributional_effects` with `summary` (one row
#'   per scenario plus the no-support baseline) and `by_quantile` (mean
#'   benefit and percentage income gain by pre-subsidy income quantile and
#'   scenario).
#'
#' @references Foster, J., Greer, J. and Thorbecke, E. (1984). A class of
#'   decomposable poverty measures. *Econometrica*, 52(3), 761-766.
#'   \doi{10.2307/1913475}
#'
#' @examples
#' de <- with(farm_subsidy, distributional_effects(
#'   farm_income,
#'   list(fertiliser = fertiliser_subsidy, transfer = income_transfer,
#'        power = power_subsidy),
#'   weights = weight, poverty_line = 60000))
#' de
#' @export
distributional_effects <- function(income, benefits, weights = NULL,
                                   poverty_line, n_quantiles = 5) {
  .chk(income, "income", scalar = FALSE)
  n <- length(income)
  if (is.numeric(benefits)) benefits <- list(scenario = benefits)
  benefits <- as.list(benefits)
  if (is.null(names(benefits)) || any(names(benefits) == "")) {
    names(benefits) <- paste0("scenario_", seq_along(benefits))
  }
  for (nm in names(benefits)) {
    .chk(benefits[[nm]], nm, scalar = FALSE)
    if (length(benefits[[nm]]) != n) {
      stop("Each benefit vector must have one value per household.",
           call. = FALSE)
    }
  }
  w <- .weights(weights, n)
  .chk(poverty_line, "poverty_line", lower = 1e-12)
  .chk(n_quantiles, "n_quantiles", lower = 1)
  gx <- .gini(income, w)
  poor0 <- income < poverty_line
  stat <- function(y, b) {
    c(mean_income = .wmean(y, w), mean_benefit = .wmean(b, w),
      gini = .gini(y, w), redistributive_effect = gx - .gini(y, w),
      fgt0 = .fgt(y, poverty_line, 0, w), fgt1 = .fgt(y, poverty_line, 1, w),
      fgt2 = .fgt(y, poverty_line, 2, w),
      lifted_out_share = .wmean(as.numeric(poor0 & y >= poverty_line), w),
      cost_per_household_lifted = {
        lifted <- sum(w * (poor0 & y >= poverty_line))
        if (lifted > 0) sum(w * b) / lifted else NA_real_
      })
  }
  rows <- rbind(no_support = stat(income, rep(0, n)),
                t(vapply(benefits, function(b) stat(income + b, b),
                         numeric(9))))
  summary <- data.frame(scenario = rownames(rows), rows, row.names = NULL)
  q <- .wq_groups(income, w, n_quantiles)
  bq <- do.call(rbind, lapply(names(benefits), function(nm) {
    b <- benefits[[nm]]
    do.call(rbind, lapply(seq_len(n_quantiles), function(k) {
      i <- q == k
      data.frame(scenario = nm, quantile = k,
                 mean_income = .wmean(income[i], w[i]),
                 mean_benefit = .wmean(b[i], w[i]),
                 benefit_share = sum(w[i] * b[i]) / sum(w * b),
                 income_gain_pct = 100 * sum(w[i] * b[i]) /
                   sum(w[i] * income[i]))
    }))
  }))
  out <- list(summary = summary, by_quantile = bq,
              poverty_line = poverty_line)
  class(out) <- "distributional_effects"
  out
}

#' @export
print.distributional_effects <- function(x, ...) {
  cat("Distributional effects of support (poverty line ",
      .fmt(x$poverty_line, 0), ")\n\n", sep = "")
  s <- x$summary
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], .fmt)
  print(s, row.names = FALSE)
  invisible(x)
}
