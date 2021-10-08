data(hiv)

test_that("subString()", {
    expect_true(as.character(subString(hiv,c(1,2),c(3,4))[1]) == "GTATAG")
    expect_true(as.character(subString(Biostrings::DNAStringSet("ATG"),1,2))
    == "AT")
    expect_true(as.character(subString(Biostrings::DNAString("ATG"),1,2))
    == "AT")
})
