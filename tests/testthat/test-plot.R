data("meuse.all", package = "gstat")

df_spdf <- meuse.all
sp::coordinates(df_spdf) <- ~ x + y
df <- meuse.all
N <- 1000

test_that("fitting model work", {
  fit_with_coords <<- DSSP(
    formula = log(zinc) ~ 1, data = df, N = N, function(x) -2 * log(1 + x),
    pars = c(0.001, 0.001), coords = ~ x + y
  )
  
  expect_true(class(fit_with_coords)=="dsspMod")
  
  
  fit_spdf <<- DSSP(
    formula = log(zinc) ~ 1, data = df_spdf, N = N, function(x) -2 * log(1 + x),
    pars = c(0.001, 0.001)
  )
  
  
  expect_true(class(fit_spdf)=="dsspMod")
  
})

test_that("plot runs without errors", {
  expect_output(plot(fit_with_coords), regexp = NA)
  expect_output(plot(fit_spdf), regexp = NA)
})
