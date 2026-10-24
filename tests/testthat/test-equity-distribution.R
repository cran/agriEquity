test_that("Theil decomposition adds up and scorecard works", {
  se <- subsidy_equity(farm_subsidy, "fertiliser_subsidy",
                       c("size_class", "region"), weights = "weight",
                       target = list(size_class = c("marginal", "small")))
  g <- se$groups
  expect_equal(c(tapply(g$benefit_share, g$dimension, sum)),
               c(region = 1, size_class = 1))
  d <- farm_subsidy
  th <- agriEquity:::.theil_decomp(d$fertiliser_subsidy, d$region, d$weight)
  expect_equal(th[["total"]], th[["between"]] + th[["within"]])
  s <- se$summary
  expect_lt(s$target_parity_ratio[s$dimension == "size_class"], 1)
  expect_true(all(s$theil_between_share >= 0 & s$theil_between_share <= 1))
  expect_output(print(se), "Summary by dimension")
  expect_error(subsidy_equity(farm_subsidy, "fertiliser_subsidy", "region",
                              reference = list(region = "Mars")))
})

test_that("distributional effects compute FGT correctly", {
  y <- c(50, 80, 120, 200)
  de <- distributional_effects(y, list(a = c(60, 0, 0, 0), b = rep(0, 4)),
                               poverty_line = 100)
  s <- de$summary
  expect_equal(s$fgt0, c(0.5, 0.25, 0.5))
  expect_equal(s$fgt1[1], mean(pmax((100 - y) / 100, 0)))
  expect_equal(s$lifted_out_share[2:3], c(0.25, 0))
  de2 <- distributional_effects(y, c(60, 0, 0, 0), poverty_line = 100)
  expect_equal(de2$summary$lifted_out_share[2], 0.25)
  expect_equal(de2$summary$cost_per_household_lifted[2], 60)
  expect_equal(nrow(de$by_quantile), 10)
  expect_output(print(de), "no_support")
  expect_error(distributional_effects(y, 1:3, poverty_line = 1))
})
