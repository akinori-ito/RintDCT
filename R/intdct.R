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

intDCT2 <- function(x1,x2,Ch=NULL) {
  n <- length(x1)
  if (is.null(Ch)) Ch <- DCT4(n)
  y1 <- rep(0,n)
  y2 <- rep(0,n)
  z <- floor(Ch %*% x2)+x1
  y1 <- floor(Ch %*% z)-x2
  y2 <- -floor(Ch %*% y1)+z
  matrix(c(y1,y2),nrow=2,ncol=n,byrow=TRUE)
}

#' intDCT: the Integer DCT
#' @export
#' @param x the input integer matrix
#' @return transformed integer vector
#' @examples
#' \dontrun{
#' library(RintDCT)
#' x <- matrix(floor(runif(100)*100),ncol=4)
#' Ch <- DCT4(50)
#' y <- intDCT(x,Ch)
#' }
#'
intDCT <- function(x) {
  dd <- dim(x)
  L <- dd[1]
  N <- dd[2]
  was_odd <- FALSE
  Ch <- DCT4(N)
  if (L %% 2 == 1) {
    L <- L+1
    x <- rbind(x,rep(0,N))
    was_odd <- TRUE
  }
  y <- matrix(0,nrow=L,ncol=N)
  i <- 1
  while (i < L) {
    y[i:(i+1),] <- intDCT2(x[i,],x[i+1,],Ch)
    i <- i+2
  }
  if (was_odd) {
    return(y[1:(L-1),])
  }
  return(y)
}

intIDCT2 <- function(y1,y2,Ch=NULL) {
  n <- length(y1)
  if (is.null(Ch)) Ch <- DCT4(n)
  x1 <- rep(0,n)
  x2 <- rep(0,n)
  z <- floor(Ch %*% y1)+y2
  x2 <- floor(Ch %*% z)-y1
  x1 <- -floor(Ch %*% x2)+z
  matrix(c(x1,x2),nrow=2,ncol=n,byrow=TRUE)
}
#' intIDCT: the Integer Inverse DCT
#' @export
#' @param x the input integer matrix
#' @return transformed integer vector
#' @examples
#' \dontrun{
#' library(RintDCT)
#' y <- matrix(floor(runif(100)*100),ncol=4)
#' Ch <- DCT4(50)
#' x <- intDCT(y,Ch)
#' }
#'
intIDCT <- function(x) {
  dd <- dim(x)
  L <- dd[1]
  N <- dd[2]
  was_odd <- FALSE
  Ch <- DCT4(N)
  if (L %% 2 == 1) {
    L <- L+1
    x <- rbind(x,rep(0,N))
    was_odd <- TRUE
  }
  y <- matrix(0,nrow=L,ncol=N)
  i <- 1
  while (i < L) {
    y[i:(i+1),] <- intIDCT2(x[i,],x[i+1,],Ch)
    i <- i+2
  }
  if (was_odd) {
    return(y[1:(L-1),])
  }
  return(y)
}


