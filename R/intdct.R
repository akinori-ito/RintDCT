#' DCT4
#' Generates DCT-IV matrix.
#' @export
#' @param n matrix dimension
#' @return n by n matrix for DCT-IV
DCT4 <- function(n) {
  x <- matrix(0,nrow=n,ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      x[i,j] <- sqrt(2/n)*cos(pi*(i-0.5)*(j-0.5)/n)
    }
  }
  x
}

#' intDCT: the Integer DCT
#' @export
#' @param x the input integer vector
#' @param Ch optional DCT-IV matrix, which should be half size of x
#' @return transformed integer vector
#' @examples
#' \dontrun{
#' library(RintDCT)
#' x <- floor(runif(100)*100)
#' Ch <- DCT4(50)
#' y <- intDCT(x,Ch)
#' }
#'
intDCT <- function(x,Ch=NULL) {
  n <- length(x)
  h <- n/2
  h1 <- h+1
  if (is.null(Ch)) Ch <- DCT4(h)
  y <- rep(0,n)
  z <- floor(Ch %*% x[h1:n])+x[1:h]
  y[1:h] <- floor(Ch %*% z)-x[h1:n]
  y[h1:n] <- -floor(Ch %*% y[1:h])+z
  y
}

#' intIDCT: the Integer Inverse DCT
#' @export
#' @param y the input integer vector
#' @param Ch optional DCT-IV matrix, which should be half size of x
#' @return transformed integer vector
#' @examples
#' \dontrun{
#' library(RintDCT)
#' y <- floor(runif(100)*100)
#' Ch <- DCT4(50)
#' x <- intDCT(y,Ch)
#' }
#'
intIDCT <- function(y,Ch=NULL) {
  n <- length(y)
  h <- n/2
  h1 <- h+1
  if (is.null(Ch)) Ch <- DCT4(h)
  xx <- rep(0,n)
  z <- floor(Ch %*% y[1:h])+y[h1:n]
  xx[h1:n] <- floor(Ch %*% z)-y[1:h]
  xx[1:h] <- -floor(Ch %*% xx[h1:n])+z
  xx
}

