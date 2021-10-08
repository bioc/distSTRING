data(hiv)

test_that("aastring2dist()", {
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix())
    expect_true(h$sitesUsed[1,1] == 91)
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(), threads=2)
    expect_true(h$sitesUsed[1,1] == 91)
})
