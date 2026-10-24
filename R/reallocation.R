# ---- internal model ---------------------------------------------------------
.prod_gain <- function(s, r, k) r * k * log1p(s / k)

.ede <- function(y, w, eps) {
  if (any(y <= 0)) return(NA_real_)
  if (abs(eps - 1) < 1e-12) return(exp(.wmean(log(y), w)))
  .wmean(y^(1 - eps), w)^(1 / (1 - eps))
}

.setup <- function(income, current, weights, budget, return_rate, saturation,
                   env_cost, income_floor) {
  .chk(income, "income", scalar = FALSE)
  n <- length(income)
  .chk(current, "current", lower = 0, scalar = FALSE)
  if (length(current) != n) {
    stop("`current` must have one value per household.", call. = FALSE)
  }
  w <- .weights(weights, n)
  if (is.null(budget)) budget <- sum(w * current)
  .chk(budget, "budget", lower = 1e-12)
  .chk(return_rate, "return_rate", lower = 0, scalar = FALSE)
  .chk(env_cost, "env_cost", lower = 0, scalar = FALSE)
  if (is.null(saturation)) saturation <- budget / sum(w)
  .chk(saturation, "saturation", lower = 1e-12, scalar = FALSE)
  y <- income
  if (!is.null(income_floor)) {
    .chk(income_floor, "income_floor", lower = 1e-12)
    y <- pmax(y, income_floor)
  }
  list(y = y, cur = current, w = w, B = budget, n = n,
       r = .expand(return_rate, n, "return_rate"),
       k = .expand(saturation, n, "saturation"),
       e = .expand(env_cost, n, "env_cost"), ybar = .wmean(y, w))
}

.optimal <- function(m, eps, lambda, cap) {
  if (eps > 0 && any(m$y <= 0)) {
    stop("Incomes must be positive when `inequality_aversion` > 0; ",
         "set `income_floor`.", call. = FALSE)
  }
  ub <- pmin(cap, m$B / m$w)
  ub[!is.finite(ub)] <- m$B
  if (sum(m$w * ub) < m$B * (1 - 1e-9)) {
    stop("`max_per_farm` is too low to allocate the whole budget.",
         call. = FALSE)
  }
  h <- function(s) {
    yy <- m$y + s + .prod_gain(s, m$r, m$k)
    mu_y <- if (eps == 0) 1 else (yy / m$ybar)^(-eps)
    mu_y * (1 + m$r / (1 + s / m$k)) - lambda * m$e
  }
  s_of <- function(mu) {
    h0 <- h(rep(0, m$n))
    hu <- h(ub)
    s <- numeric(m$n)
    s[hu >= mu] <- ub[hu >= mu]
    act <- which(h0 > mu & hu < mu)
    if (length(act)) {
      lo <- rep(0, length(act))
      hi <- ub[act]
      for (it in 1:60) {
        mid <- (lo + hi) / 2
        yy <- m$y[act] + mid + .prod_gain(mid, m$r[act], m$k[act])
        mu_y <- if (eps == 0) 1 else (yy / m$ybar)^(-eps)
        hv <- mu_y * (1 + m$r[act] / (1 + mid / m$k[act])) -
          lambda * m$e[act]
        up <- hv > mu
        lo[up] <- mid[up]
        hi[!up] <- mid[!up]
      }
      s[act] <- (lo + hi) / 2
    }
    s
  }
  h0 <- h(rep(0, m$n))
  hu <- h(ub)
  mu_hi <- max(h0)
  mu_lo <- min(hu) - 1
  for (it in 1:100) {
    mu <- (mu_lo + mu_hi) / 2
    if (sum(m$w * s_of(mu)) > m$B) mu_lo <- mu else mu_hi <- mu
  }
  s <- s_of(mu_lo)
  tot <- sum(m$w * s)
  if (tot > 0) s <- s * m$B / tot
  s
}

.metrics <- function(s, m, eps, lambda, poverty_line, target) {
  q <- .prod_gain(s, m$r, m$k)
  post <- m$y + s + q
  B <- sum(m$w * s)
  pg <- sum(m$w * q)
  ec <- sum(m$w * m$e * s)
  out <- c(budget = B, productivity_gain = pg,
           productivity_per_unit = if (B > 0) pg / B else NA,
           environmental_cost = ec, net_social_value = pg - ec,
           gini_post = .gini(post, m$w),
           concentration = .conc(s, m$y, m$w),
           ede_income = .ede(post, m$w, eps),
           welfare_gain_ede = .ede(post, m$w, eps) - .ede(m$y, m$w, eps))
  if (!is.null(target)) {
    out["target_share"] <- if (B > 0) sum(m$w * s * target) / B else NA
  }
  if (!is.null(poverty_line)) {
    out["poverty_headcount"] <- .fgt(post, poverty_line, 0, m$w)
    out["poverty_gap"] <- .fgt(post, poverty_line, 1, m$w)
  }
  out
}

.prop <- function(x, m) {
  if (sum(m$w * x) <= 0) stop("Allocation key sums to zero.", call. = FALSE)
  x * m$B / sum(m$w * x)
}

# ---- exported ---------------------------------------------------------------

#' Budget-neutral reallocation of agricultural support
#'
#' Redistributes a fixed support budget among farm households under
#' alternative allocation rules and evaluates every allocation on the same
#' productivity, equity and sustainability metrics, so that reforms can be
#' compared at equal fiscal cost.
#'
#' **Rules.** `"current"` (the existing distribution, rescaled to the
#' budget), `"per_farm"` (equal amount per household), `"per_hectare"`
#' (proportional to operated area), `"per_hectare_capped"` (proportional to
#' area up to `area_cap` hectares), `"targeted"` (proportional to a
#' household priority score, for example higher for marginal farmers, women
#' and resource-poor regions) and `"optimal"`.
#'
#' **Household model.** A household receiving support \eqn{s_i} gains the
#' transfer itself plus an additional productivity response with
#' diminishing returns,
#' \deqn{q_i(s_i) = r_i k_i \ln(1 + s_i / k_i),}
#' where \eqn{r_i} is the marginal output value of the first unit of
#' support (for example from better input use) and \eqn{k_i} a saturation
#' scale. Support also causes an environmental cost \eqn{e_i s_i} (for
#' example the damage value of nutrient surplus or groundwater depletion
#' induced by input or power subsidies).
#'
#' **Optimal rule.** Shares \eqn{s_i \ge 0} with \eqn{\sum_i w_i s_i = B}
#' maximise
#' \deqn{\sum_i w_i U(y_i + s_i + q_i(s_i)) - \lambda \sum_i w_i e_i s_i,}
#' where \eqn{U} has constant inequality aversion \eqn{\varepsilon}
#' (Atkinson, 1970) normalised so that the marginal welfare weight at mean
#' income is 1. \eqn{\varepsilon = 0} values all income equally (pure
#' efficiency); larger \eqn{\varepsilon} gives more weight to poorer
#' households. \eqn{\lambda} is the weight on environmental cost. The
#' solution equalises marginal social value across households and is found
#' by nested bisection on the budget multiplier.
#'
#' **Metrics.** Budget, productivity gain (total \eqn{q_i}), productivity
#' per unit of budget, environmental cost, net social value (productivity
#' gain minus environmental cost), Gini of post-support income,
#' concentration index of support by income, Atkinson equally distributed
#' equivalent income and its gain, share of support reaching `target`
#' households, and optionally poverty headcount and gap.
#'
#' @param income Pre-support household income.
#' @param current Current support received per household.
#' @param weights Optional survey weights.
#' @param budget Budget to allocate (weighted total); defaults to the
#'   current total, which makes the comparison budget-neutral.
#' @param area Optional operated area (ha) for area-based rules.
#' @param priority Optional non-negative priority score per household for
#'   the targeted rule.
#' @param target Optional logical vector marking priority households
#'   (smallholders, women, resource-poor regions) for the target share
#'   metric.
#' @param return_rate Marginal productivity return \eqn{r_i} (scalar or one
#'   per household).
#' @param saturation Saturation scale \eqn{k_i}; defaults to the budget per
#'   household.
#' @param env_cost Environmental cost per unit of support \eqn{e_i}.
#' @param inequality_aversion Atkinson inequality aversion
#'   \eqn{\varepsilon \ge 0}.
#' @param env_weight Weight \eqn{\lambda} on environmental cost.
#' @param area_cap Area cap (ha) for the capped per-hectare rule.
#' @param max_per_farm Upper limit of support per household in the optimal
#'   rule.
#' @param poverty_line Optional poverty line.
#' @param income_floor Optional floor applied to incomes (needed when some
#'   incomes are zero or negative and `inequality_aversion` > 0).
#' @param rules Allocation rules to evaluate.
#'
#' @return A list of class `budget_reallocation` with `metrics` (one row per
#'   rule), `allocations` (data frame of support per household under each
#'   rule) and `params`.
#'
#' @references Atkinson, A. B. (1970). On the measurement of inequality.
#'   *Journal of Economic Theory*, 2(3), 244-263.
#'   \doi{10.1016/0022-0531(70)90039-6}
#'
#' @examples
#' fs <- farm_subsidy
#' prio <- with(fs, 1 + (size_class %in% c("marginal", "small")) +
#'                0.5 * (gender == "female") +
#'                0.5 * (region %in% c("Rainfed-East", "Hill-Tribal")))
#' br <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
#'                           weights = fs$weight, area = fs$land_ha,
#'                           priority = prio,
#'                           target = fs$size_class %in% c("marginal", "small"),
#'                           return_rate = fs$marginal_return,
#'                           env_cost = fs$env_cost_per_rupee,
#'                           inequality_aversion = 1, env_weight = 1,
#'                           poverty_line = 60000)
#' br
#' @export
budget_reallocation <- function(income, current, weights = NULL,
                                budget = NULL, area = NULL, priority = NULL,
                                target = NULL, return_rate = 0,
                                saturation = NULL, env_cost = 0,
                                inequality_aversion = 1, env_weight = 1,
                                area_cap = 2, max_per_farm = Inf,
                                poverty_line = NULL, income_floor = NULL,
                                rules = c("current", "per_farm", "per_hectare",
                                          "per_hectare_capped", "targeted",
                                          "optimal")) {
  m <- .setup(income, current, weights, budget, return_rate, saturation,
              env_cost, income_floor)
  rules <- match.arg(rules, several.ok = TRUE)
  .chk(inequality_aversion, "inequality_aversion", lower = 0)
  .chk(env_weight, "env_weight", lower = 0)
  .chk(area_cap, "area_cap", lower = 1e-12)
  if (!is.numeric(max_per_farm) || length(max_per_farm) != 1L ||
      is.na(max_per_farm) || max_per_farm <= 0) {
    stop("`max_per_farm` must be a positive number.", call. = FALSE)
  }
  if (!is.null(area)) {
    .chk(area, "area", lower = 0, scalar = FALSE)
    if (length(area) != m$n) stop("`area` must have one value per household.",
                                  call. = FALSE)
  }
  if (!is.null(priority)) {
    .chk(priority, "priority", lower = 0, scalar = FALSE)
    if (length(priority) != m$n) {
      stop("`priority` must have one value per household.", call. = FALSE)
    }
  }
  if (!is.null(target)) {
    if (!is.logical(target) || length(target) != m$n || anyNA(target)) {
      stop("`target` must be a logical vector, one value per household.",
           call. = FALSE)
    }
  }
  if (!is.null(poverty_line)) .chk(poverty_line, "poverty_line", lower = 0)
  if (is.null(area)) rules <- setdiff(rules, c("per_hectare",
                                               "per_hectare_capped"))
  if (is.null(priority)) rules <- setdiff(rules, "targeted")

  alloc <- list()
  for (rl in rules) {
    alloc[[rl]] <- switch(rl,
      current = .prop(m$cur, m),
      per_farm = rep(m$B / sum(m$w), m$n),
      per_hectare = .prop(area, m),
      per_hectare_capped = .prop(pmin(area, area_cap), m),
      targeted = .prop(priority, m),
      optimal = .optimal(m, inequality_aversion, env_weight, max_per_farm))
  }
  met <- t(vapply(alloc, .metrics, numeric(length(.metrics(
    alloc[[1]], m, inequality_aversion, env_weight, poverty_line, target))),
    m = m, eps = inequality_aversion, lambda = env_weight,
    poverty_line = poverty_line, target = target))
  metrics <- data.frame(rule = rownames(met), met, row.names = NULL)
  out <- list(metrics = metrics, allocations = as.data.frame(alloc),
              params = list(budget = m$B,
                            inequality_aversion = inequality_aversion,
                            env_weight = env_weight))
  class(out) <- "budget_reallocation"
  out
}

#' @export
print.budget_reallocation <- function(x, ...) {
  cat("Budget-neutral reallocation of support (budget ",
      .fmt(x$params$budget, 0), "; inequality aversion ",
      x$params$inequality_aversion, "; environmental weight ",
      x$params$env_weight, ")\n\n", sep = "")
  s <- x$metrics
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], .fmt)
  print(s, row.names = FALSE)
  invisible(x)
}

#' Productivity-equity-sustainability trade-off frontier
#'
#' Solves the welfare-optimal budget-neutral allocation of
#' [budget_reallocation()] over a grid of inequality-aversion and
#' environmental-weight values and records the resulting productivity gain,
#' environmental cost, post-support Gini, target-group share and poverty.
#' Plotting productivity against equity for each environmental weight shows
#' how much output must be given up for a more equal, or a less
#' environmentally damaging, allocation of the same budget, and where the
#' current allocation lies relative to the frontier.
#'
#' @inheritParams budget_reallocation
#' @param aversion Values of inequality aversion \eqn{\varepsilon}.
#' @param env_weights Values of the environmental weight \eqn{\lambda}.
#'
#' @return A list of class `tradeoff_frontier` with `frontier` (one row per
#'   combination of `aversion` and `env_weights`) and `current` (metrics of
#'   the current allocation). A `plot()` method draws productivity gain
#'   against the Gini of post-support income.
#'
#' @examples
#' fs <- farm_subsidy[1:400, ]
#' tf <- tradeoff_frontier(fs$farm_income, fs$fertiliser_subsidy,
#'                         weights = fs$weight,
#'                         return_rate = fs$marginal_return,
#'                         env_cost = fs$env_cost_per_rupee,
#'                         aversion = c(0, 1, 3), env_weights = c(0, 2))
#' tf
#' plot(tf)
#' @export
tradeoff_frontier <- function(income, current, weights = NULL, budget = NULL,
                              target = NULL, return_rate = 0,
                              saturation = NULL, env_cost = 0,
                              aversion = c(0, 0.5, 1, 2, 3),
                              env_weights = c(0, 1, 2),
                              max_per_farm = Inf, poverty_line = NULL,
                              income_floor = NULL) {
  m <- .setup(income, current, weights, budget, return_rate, saturation,
              env_cost, income_floor)
  .chk(aversion, "aversion", lower = 0, scalar = FALSE)
  .chk(env_weights, "env_weights", lower = 0, scalar = FALSE)
  grid <- expand.grid(inequality_aversion = aversion, env_weight = env_weights)
  res <- lapply(seq_len(nrow(grid)), function(i) {
    e <- grid$inequality_aversion[i]
    l <- grid$env_weight[i]
    s <- .optimal(m, e, l, max_per_farm)
    .metrics(s, m, e, l, poverty_line, target)
  })
  frontier <- cbind(grid, do.call(rbind, res))
  cur <- .metrics(.prop(m$cur, m), m, 1, 0, poverty_line, target)
  out <- list(frontier = frontier, current = cur)
  class(out) <- "tradeoff_frontier"
  out
}

#' @export
print.tradeoff_frontier <- function(x, ...) {
  cat("Productivity-equity-sustainability frontier\n\n")
  f <- x$frontier
  keep <- intersect(c("inequality_aversion", "env_weight", "productivity_gain",
                      "environmental_cost", "net_social_value", "gini_post",
                      "target_share", "poverty_headcount"), names(f))
  f <- f[, keep]
  f[-(1:2)] <- lapply(f[-(1:2)], .fmt)
  print(f, row.names = FALSE)
  cat("\nCurrent allocation: productivity gain ",
      .fmt(x$current[["productivity_gain"]]), ", environmental cost ",
      .fmt(x$current[["environmental_cost"]]), ", Gini ",
      .fmt(x$current[["gini_post"]]), "\n", sep = "")
  invisible(x)
}

#' Plot a productivity-equity-sustainability frontier
#'
#' @param x A [tradeoff_frontier()] object.
#' @param ... Further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' fs <- farm_subsidy[1:300, ]
#' plot(tradeoff_frontier(fs$farm_income, fs$fertiliser_subsidy,
#'                        weights = fs$weight, return_rate = 0.3,
#'                        aversion = c(0, 2), env_weights = 0))
#' @export
plot.tradeoff_frontier <- function(x, ...) {
  f <- x$frontier
  lv <- sort(unique(f$env_weight))
  cols <- seq_along(lv) + 1
  graphics::plot(range(c(f$gini_post, x$current[["gini_post"]])),
                 range(c(f$productivity_gain,
                         x$current[["productivity_gain"]])),
                 type = "n", xlab = "Gini of post-support income",
                 ylab = "Productivity gain", ...)
  for (i in seq_along(lv)) {
    g <- f[f$env_weight == lv[i], ]
    g <- g[order(g$inequality_aversion), ]
    graphics::lines(g$gini_post, g$productivity_gain, col = cols[i], lwd = 2)
    graphics::points(g$gini_post, g$productivity_gain, col = cols[i],
                     pch = 19)
  }
  graphics::points(x$current[["gini_post"]], x$current[["productivity_gain"]],
                   pch = 8, cex = 1.5)
  graphics::legend("topright",
                   legend = c(paste("env. weight", lv), "current"),
                   col = c(cols, 1), pch = c(rep(19, length(lv)), 8),
                   bty = "n", cex = 0.8)
  invisible(x)
}
