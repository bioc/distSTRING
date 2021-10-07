data(hiv)

test_that("cds2aa() outputs AAStringSet", {
  expect_true(class(cds2aa(hiv)) == "AAStringSet")
})
