.group_var <- function(data, group, rank_var, n_quantiles, w) {
  if (!is.null(group)) {
    if (!group %in% names(data)) {
      stop("Column `", group, "` not found in `data`.", call. = FALSE)
    }
    return(as.character(data[[group]]))
  }
  if (is.null(rank_var)) {
    stop("Supply `group` or `rank_var`.", call. = FALSE)
  }
  q <- .wq_groups(data[[rank_var]], w, n_quantiles)
  paste0("Q", q)
}

#' Benefit incidence analysis of agricultural support
#'
#' Shows how a subsidy or support payment is distributed across groups of
#' farm households: landholding classes, gender of the farmer, regions, or
#' quantiles of income or landholding. For each group the function reports
#' its (weighted) share of farm households, its share of total benefits,
#' the parity ratio (benefit share divided by population share; values
#' below 1 mean the group receives less than its population share), the
#' coverage rate (share of households receiving any benefit), the mean
#' benefit per household and per beneficiary, and, if supplied, the mean
#' benefit per hectare and benefits as a share of group income.
#'
#' @param data A data frame with one row per farm household.
#' @param benefit Name of the column with the benefit received (currency
#'   per household per year; zero for non-recipients).
#' @param group Name of a grouping column (for example landholding class).
#'   Alternatively use `rank_var`.
#' @param rank_var Name of a column (for example income or land) used to
#'   form weighted quantile groups when `group` is `NULL`.
#' @param n_quantiles Number of quantile groups for `rank_var`.
#' @param weights Optional name of a survey-weight column.
#' @param area Optional name of an operated-area column (hectares).
#' @param income Optional name of an income column.
#'
#' @return A data frame of class `benefit_incidence`, one row per group,
#'   with columns `group`, `households`, `population_share`,
#'   `benefit_share`, `parity_ratio`, `coverage`,
#'   `mean_benefit_per_household`, `mean_benefit_per_beneficiary` and,
#'   when available, `benefit_per_hectare` and `benefit_income_ratio`.
#'
#' @examples
#' benefit_incidence(farm_subsidy, "fertiliser_subsidy",
#'                   group = "size_class", weights = "weight",
#'                   area = "land_ha", income = "farm_income")
#' benefit_incidence(farm_subsidy, "income_transfer",
#'                   rank_var = "farm_income", weights = "weight")
#' @export
benefit_incidence <- function(data, benefit, group = NULL, rank_var = NULL,
                              n_quantiles = 5, weights = NULL, area = NULL,
                              income = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  miss <- setdiff(c(benefit, rank_var, weights, area, income), names(data))
  if (length(miss)) {
    stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  b <- data[[benefit]]
  .chk(b, benefit, lower = 0, scalar = FALSE)
  w <- .weights(if (is.null(weights)) NULL else data[[weights]], nrow(data))
  g <- .group_var(data, group, rank_var, n_quantiles, w)
  lev <- if (is.null(group)) paste0("Q", seq_len(n_quantiles)) else
    unique(g)
  lev <- lev[lev %in% g]
  tb <- sum(w * b)
  tw <- sum(w)
  rows <- lapply(lev, function(l) {
    i <- g == l
    wi <- w[i]
    bi <- b[i]
    out <- data.frame(
      group = l, households = sum(wi), population_share = sum(wi) / tw,
      benefit_share = if (tb > 0) sum(wi * bi) / tb else NA_real_,
      coverage = sum(wi * (bi > 0)) / sum(wi),
      mean_benefit_per_household = sum(wi * bi) / sum(wi),
      mean_benefit_per_beneficiary = if (any(bi > 0))
        sum(wi * bi) / sum(wi * (bi > 0)) else 0,
      stringsAsFactors = FALSE)
    if (!is.null(area)) {
      a <- data[[area]][i]
      out$benefit_per_hectare <- if (sum(wi * a) > 0)
        sum(wi * bi) / sum(wi * a) else NA_real_
    }
    if (!is.null(income)) {
      y <- data[[income]][i]
      out$benefit_income_ratio <- if (sum(wi * y) > 0)
        sum(wi * bi) / sum(wi * y) else NA_real_
    }
    out
  })
  res <- do.call(rbind, rows)
  res$parity_ratio <- res$benefit_share / res$population_share
  first <- c("group", "households", "population_share", "benefit_share",
             "parity_ratio")
  res <- res[, c(first, setdiff(names(res), first))]
  rownames(res) <- NULL
  class(res) <- c("benefit_incidence", "data.frame")
  attr(res, "benefit") <- benefit
  res
}

#' Progressivity of agricultural support
#'
#' Measures whether a benefit is distributed progressively with respect to
#' a welfare ranking variable such as household income or landholding.
#' \describe{
#'   \item{Concentration index \eqn{C_b}}{of the benefit, ranking households
#'     by `rank_by`. Negative values mean the benefit is concentrated among
#'     the poorer (or smaller) households in absolute terms.}
#'   \item{Gini \eqn{G_x}}{of pre-subsidy income.}
#'   \item{Kakwani index}{\eqn{K = G_x - C_b} (Kakwani, 1977), positive when
#'     the benefit is distributed more equally than income, i.e. relatively
#'     progressive.}
#'   \item{Redistributive effect}{\eqn{RE = G_x - G_{x+b}}, decomposed into
#'     the vertical (Reynolds-Smolensky) effect \eqn{G_x - C_{x+b}} and
#'     re-ranking \eqn{G_{x+b} - C_{x+b}}.}
#' }
#' Concentration and Gini indices are computed with survey weights using
#' the covariance formula. Optional bootstrap confidence intervals resample
#' households with replacement (set a seed beforehand for reproducibility).
#'
#' @param benefit Numeric vector of benefits per household.
#' @param income Numeric vector of pre-subsidy income per household.
#' @param rank_by Ranking variable for the concentration index; defaults to
#'   `income`. Use landholding to assess progressivity by farm size.
#' @param weights Optional survey weights.
#' @param n_boot Number of bootstrap replications (0 for none).
#' @param conf Confidence level for bootstrap intervals.
#'
#' @return A list of class `progressivity_index` with `estimates` (named
#'   vector of the indices) and, if `n_boot > 0`, `ci` (bootstrap
#'   percentile intervals) and an `interpretation` string.
#'
#' @references Kakwani, N. C. (1977). Measurement of tax progressivity: an
#'   international comparison. *The Economic Journal*, 87(345), 71-80.
#'   \doi{10.2307/2230485}
#'
#'   Reynolds, M. and Smolensky, E. (1977). *Public Expenditures, Taxes,
#'   and the Distribution of Income: The United States, 1950, 1961, 1970*.
#'   Academic Press, New York.
#'
#' @examples
#' with(farm_subsidy, progressivity_index(fertiliser_subsidy, farm_income,
#'                                        weights = weight))
#' set.seed(1)
#' with(farm_subsidy, progressivity_index(income_transfer, farm_income,
#'                                        weights = weight, n_boot = 99))
#' @export
progressivity_index <- function(benefit, income, rank_by = income,
                                weights = NULL, n_boot = 0, conf = 0.95) {
  n <- length(benefit)
  .chk(benefit, "benefit", lower = 0, scalar = FALSE)
  .chk(income, "income", scalar = FALSE)
  .chk(rank_by, "rank_by", scalar = FALSE)
  if (length(income) != n || length(rank_by) != n) {
    stop("`benefit`, `income` and `rank_by` must have equal length.",
         call. = FALSE)
  }
  w <- .weights(weights, n)
  .chk(n_boot, "n_boot", lower = 0)
  .chk(conf, "conf", 0.5, 0.999)
  est_fun <- function(idx) {
    b <- benefit[idx]
    x <- income[idx]
    r <- rank_by[idx]
    ww <- w[idx]
    post <- x + b
    cb <- .conc(b, r, ww)
    gx <- .gini(x, ww)
    gp <- .gini(post, ww)
    cp <- .conc(post, x, ww)
    c(concentration = cb, gini_income = gx, kakwani = gx - cb,
      gini_post = gp, redistributive_effect = gx - gp,
      reynolds_smolensky = gx - cp, reranking = gp - cp,
      benefit_share_bottom40 = {
        rr <- .frank(r, ww)
        sum(ww * b * (rr <= 0.4)) / sum(ww * b)
      })
  }
  est <- est_fun(seq_len(n))
  out <- list(estimates = est)
  if (n_boot > 0) {
    reps <- vapply(seq_len(n_boot), function(i)
      est_fun(sample.int(n, n, replace = TRUE)), numeric(length(est)))
    a <- (1 - conf) / 2
    out$ci <- t(apply(reps, 1, stats::quantile, probs = c(a, 1 - a),
                      na.rm = TRUE))
  }
  out$interpretation <- paste0(
    if (est[["concentration"]] < 0) "Pro-poor in absolute terms"
    else "Not pro-poor in absolute terms",
    "; ",
    if (est[["kakwani"]] > 0) "progressive relative to income"
    else "regressive relative to income", ".")
  class(out) <- "progressivity_index"
  out
}

#' @export
print.progressivity_index <- function(x, ...) {
  cat("Progressivity of support\n\n")
  tab <- data.frame(index = names(x$estimates),
                    estimate = .fmt(x$estimates, 4), row.names = NULL)
  if (!is.null(x$ci)) {
    tab$lower <- .fmt(x$ci[, 1], 4)
    tab$upper <- .fmt(x$ci[, 2], 4)
  }
  print(tab, row.names = FALSE)
  cat("\n", x$interpretation, "\n", sep = "")
  invisible(x)
}
