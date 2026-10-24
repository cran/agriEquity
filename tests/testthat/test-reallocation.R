fs <- farm_subsidy[1:300, ]

test_that("all rules are budget-neutral", {
  br <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
                            weights = fs$weight, area = fs$land_ha,
                            priority = rep(1, 300),
                            return_rate = fs$marginal_return,
                            env_cost = fs$env_cost_per_rupee)
  B <- sum(fs$weight * fs$fertiliser_subsidy)
  expect_equal(br$metrics$budget, rep(B, 6), tolerance = 1e-8)
  expect_true(all(br$allocations >= 0))
  expect_equal(br$allocations$targeted, br$allocations$per_farm)
  # optimal should beat every rule on its own objective proxy (EDE gain
  # with no environmental cost when lambda = 0 and r = 0)
  br0 <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
                             weights = fs$weight, inequality_aversion = 2,
                             env_weight = 0)
  m <- br0$metrics
  expect_equal(which.max(m$welfare_gain_ede), which(m$rule == "optimal"))
  expect_output(print(br), "Budget-neutral")
})

test_that("pure efficiency allocation goes to highest net returns", {
  y <- rep(100, 4)
  br <- budget_reallocation(y, rep(10, 4), return_rate = c(0.1, 0.2, 0.3,
                                                           0.9),
                            saturation = 1e6, inequality_aversion = 0,
                            env_weight = 0, rules = "optimal")
  a <- br$allocations$optimal
  expect_equal(which.max(a), 4L)
  expect_equal(sum(a), 40, tolerance = 1e-8)
})

test_that("inequality aversion equalises post-support income", {
  y <- c(50, 100, 150, 200)
  br <- budget_reallocation(y, rep(25, 4), inequality_aversion = 5,
                            env_weight = 0, rules = "optimal")
  post <- y + br$allocations$optimal
  expect_equal(post[1:2], c(125, 125), tolerance = 1e-3)
  expect_equal(br$allocations$optimal[3:4], c(0, 0), tolerance = 1e-6)
})

test_that("environmental weight reduces environmental cost", {
  a <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
                           weights = fs$weight, env_cost = fs$env_cost_per_rupee,
                           return_rate = 0.2, env_weight = 0, rules = "optimal")
  b <- budget_reallocation(fs$farm_income, fs$fertiliser_subsidy,
                           weights = fs$weight, env_cost = fs$env_cost_per_rupee,
                           return_rate = 0.2, env_weight = 3, rules = "optimal")
  expect_lt(b$metrics$environmental_cost, a$metrics$environmental_cost)
})

test_that("caps and input validation", {
  br <- budget_reallocation(c(10, 20, 30), c(5, 5, 5), max_per_farm = 6,
                            inequality_aversion = 3, env_weight = 0,
                            rules = "optimal")
  expect_true(all(br$allocations$optimal <= 6 + 1e-8))
  expect_error(budget_reallocation(c(10, 20), c(5, 5), max_per_farm = 1,
                                   rules = "optimal"))
  expect_error(budget_reallocation(c(-5, 20), c(5, 5), rules = "optimal"))
  expect_silent(budget_reallocation(c(-5, 20), c(5, 5), income_floor = 1,
                                    rules = "optimal"))
  expect_error(budget_reallocation(1:3, 1:3, target = c(TRUE, FALSE)))
})

test_that("trade-off frontier runs and plots", {
  tf <- tradeoff_frontier(fs$farm_income, fs$fertiliser_subsidy,
                          weights = fs$weight, return_rate = fs$marginal_return,
                          env_cost = fs$env_cost_per_rupee,
                          aversion = c(0, 2), env_weights = c(0, 3))
  expect_equal(nrow(tf$frontier), 4)
  f <- tf$frontier
  expect_lte(f$gini_post[f$inequality_aversion == 2 & f$env_weight == 0],
             f$gini_post[f$inequality_aversion == 0 & f$env_weight == 0])
  expect_output(print(tf), "Current allocation")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(tf))
})
