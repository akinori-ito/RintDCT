test_that("DCT/IDCT pair works", {
  N <- 1024
  L <- 100
  test_intdct <- function(N) {
    x <- matrix(floor(runif(N*L)*1000),nrow=L,ncol=N)
    y <- intDCT(x)
    xx <- intIDCT(y)
    if (any(xx!=x)) {
      return(FALSE)
    }
    TRUE
  }
  expect_equal(test_intdct(N),TRUE)
})
