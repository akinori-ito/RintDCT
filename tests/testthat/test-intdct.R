test_that("DCT/IDCT pair works", {
  N <- 2048
  Ch <- DCT4(N/2)
  test_intdct <- function(N,Ch) {
    x <- floor(runif(N)*1000)
    y <- intDCT(x,Ch)
    xx <- intIDCT(y,Ch)
    if (any(xx!=x)) {
      return(FALSE)
    }
    TRUE
  }
  expect_equal(test_intdct(N,Ch),TRUE)
})
