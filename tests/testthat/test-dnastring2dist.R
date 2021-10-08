data(hiv)

test_that("compareCodons()", {
    myscore <- iupacMatrix()
    myscore[1,1] <- -1
    h <- hiv |> dnastring2dist(score=myscore)
    expect_true(h$sitesUsed[1,1] == 150)
    h <- hiv |> dnastring2dist(model="IUPAC")
    expect_true(h$sitesUsed[1,1] == 273)
    h <- hiv |> dnastring2dist(model="K80")
    expect_true(h$sitesUsed[1,1] == 273)
})
