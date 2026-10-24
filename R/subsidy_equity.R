.theil_decomp <- function(b, g, w) {
  mu <- .wmean(b, w)
  if (mu <= 0) return(c(total = NA, between = NA, within = NA))
  tl <- function(x, ww) {
    m <- .wmean(x, ww)
    if (m <= 0) return(0)
    r <- x / m
    .wmean(ifelse(r > 0, r * log(r), 0), ww)
  }
  total <- tl(b, w)
  between <- 0
  within <- 0
  W <- sum(w)
  for (l in unique(g)) {
    i <- g == l
    mg <- .wmean(b[i], w[i])
    sg <- sum(w[i]) * mg / (W * mu)
    if (mg > 0) {
      between <- between + sg * log(mg / mu)
      within <- within + sg * tl(b[i], w[i])
    }
  }
  c(total = total, between = between, within = within)
}

#' Multi-dimensional equity scorecard for agricultural support
#'
#' Assesses a subsidy along several equity dimensions at once, for
#' example landholding class, gender of the farmer and region. For each
#' dimension and group it reports population and benefit shares, the parity
#' ratio (benefit share / population share), coverage, and mean benefit per
#' household relative to a reference group. Inequality of benefits across
#' households is measured by the Theil index and decomposed into
#' between-group and within-group components for each dimension, so the
#' share of benefit inequality explained by farm size, gender or region can
#' be compared.
#'
#' A summary row per dimension gives the minimum parity ratio (the most
#' under-served group), the between-group Theil share and the target-group
#' parity when `target` groups are named.
#'
#' @param data A data frame with one row per household.
#' @param benefit Name of the benefit column.
#' @param dimensions Character vector of grouping columns (for example
#'   `c("size_class", "gender", "region")`).
#' @param weights Optional name of a survey-weight column.
#' @param reference Optional named list giving, for each dimension, the
#'   reference group for relative benefits (default: the group with the
#'   highest mean benefit).
#' @param target Optional named list giving, for each dimension, the
#'   priority groups (for example smallholders, women, resource-poor
#'   regions) whose combined parity ratio is reported.
#'
#' @return A list of class `subsidy_equity` with `groups` (data frame by
#'   dimension and group) and `summary` (one row per dimension).
#'
#' @examples
#' se <- subsidy_equity(farm_subsidy, "fertiliser_subsidy",
#'                      dimensions = c("size_class", "gender", "region"),
#'                      weights = "weight",
#'                      target = list(size_class = c("marginal", "small"),
#'                                    gender = "female",
#'                                    region = c("Rainfed-East",
#'                                               "Hill-Tribal")))
#' se
#' @export
subsidy_equity <- function(data, benefit, dimensions, weights = NULL,
                           reference = NULL, target = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  miss <- setdiff(c(benefit, dimensions, weights), names(data))
  if (length(miss)) {
    stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  b <- data[[benefit]]
  .chk(b, benefit, lower = 0, scalar = FALSE)
  w <- .weights(if (is.null(weights)) NULL else data[[weights]], nrow(data))
  tb <- sum(w * b)
  tw <- sum(w)
  groups <- list()
  smry <- list()
  for (dname in dimensions) {
    g <- as.character(data[[dname]])
    lev <- unique(g)
    tab <- do.call(rbind, lapply(lev, function(l) {
      i <- g == l
      data.frame(dimension = dname, group = l,
                 population_share = sum(w[i]) / tw,
                 benefit_share = if (tb > 0) sum(w[i] * b[i]) / tb else NA,
                 coverage = sum(w[i] * (b[i] > 0)) / sum(w[i]),
                 mean_benefit = sum(w[i] * b[i]) / sum(w[i]),
                 stringsAsFactors = FALSE)
    }))
    tab$parity_ratio <- tab$benefit_share / tab$population_share
    ref <- if (!is.null(reference[[dname]])) reference[[dname]] else
      tab$group[which.max(tab$mean_benefit)]
    if (!ref %in% tab$group) {
      stop("Reference group `", ref, "` not found in `", dname, "`.",
           call. = FALSE)
    }
    tab$relative_to_reference <- tab$mean_benefit /
      tab$mean_benefit[tab$group == ref]
    groups[[dname]] <- tab
    th <- .theil_decomp(b, g, w)
    tp <- NA_real_
    if (!is.null(target[[dname]])) {
      tg <- target[[dname]]
      if (!all(tg %in% tab$group)) {
        stop("Target group(s) not found in `", dname, "`.", call. = FALSE)
      }
      k <- tab$group %in% tg
      tp <- sum(tab$benefit_share[k]) / sum(tab$population_share[k])
    }
    smry[[dname]] <- data.frame(
      dimension = dname, groups = nrow(tab), reference = ref,
      min_parity_group = tab$group[which.min(tab$parity_ratio)],
      min_parity_ratio = min(tab$parity_ratio),
      target_parity_ratio = tp,
      theil_total = th[["total"]],
      theil_between_share = if (isTRUE(th[["total"]] > 0))
        th[["between"]] / th[["total"]] else NA_real_,
      stringsAsFactors = FALSE)
  }
  out <- list(groups = do.call(rbind, groups),
              summary = do.call(rbind, smry), benefit = benefit)
  rownames(out$groups) <- NULL
  rownames(out$summary) <- NULL
  class(out) <- "subsidy_equity"
  out
}

#' @export
print.subsidy_equity <- function(x, ...) {
  cat("Equity scorecard for `", x$benefit, "`\n\n", sep = "")
  g <- x$groups
  num <- vapply(g, is.numeric, logical(1))
  g[num] <- lapply(g[num], .fmt)
  print(g, row.names = FALSE)
  cat("\nSummary by dimension\n")
  s <- x$summary
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], .fmt)
  print(s, row.names = FALSE)
  invisible(x)
}
