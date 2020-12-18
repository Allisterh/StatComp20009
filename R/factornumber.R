#' @title Determining The Number of Factors In Approximate Factor Model
#'
#' @description This function approximates the number of factors in an approximate factor model
#'  for large N by T matrices using the methods and criteria found in Bai and Ng (2002)
#'
#' @usage getnfac(x,kmax,criteria)
#'
#' @param x A matrix containing the data.
#'
#' @param kmax An integer with the maximum number of common factors to search over. This methodology
#' is weak to underestimation of the number of common factors so setting this value higher is preferred.
#'
#' @param criteria a character vector of length one with values of either PC1,PC2,PC3, or eigen.
#'  Choosing eigen makes the number of factors equal to the number of columns whose sum of eigenvalues is less than  or equal to .5.
#'
#' @details This function approximates the number of factors in an approximate
#' factor model. Using the penalty functions PC(1),PC(2),PC(3).
#'
#' @return ic Integer of the approximate number of factors based off of the chosen
#' penalty function
#'
#' @return lambda A matrix of the estimated factor loadings associated with common factors.
#'
#' @return Fhat A matrix of the estimated common components
#'
#'
#' @references Jushan Bai and Serena Ng. 'Determining the Number of Factors in
#' Approximate Factor Models.' Econometrica 70.1 (2002): 191-221. Print.
#'
#'@export
getnfac <- function(x, kmax = NULL, criteria = NULL) {
  
  #ic : Integer of the approximate number of factors based off of the chosen
  #penalty function
  #
  #lambda : A matrix of the estimated factor loadings associated with common factors.
  #
  # Fhat : A matrix of the estimated common components
  
  if (!(kmax %% 1 == 0)){
    stop(" k1 must be an integer.")
  }
  if (is.null(criteria)){
    warning("criteria is NULL, setting criteria to PC3")
    criteria <- "PC3"
  }
  
  Tn <- dim(x)[1]
  
  N <- dim(x)[2]
  
  NT <- N * Tn
  
  NT1 <- N + Tn
  
  CT <- matrix(0, 1, kmax)
  
  ii <- seq(1:kmax)
  
  GCT <- min(N, Tn)
  
  
  if(is.null(criteria)){
    warning("Criteria is NULL, setting to BIC3")
    criteria <- "PC3"
  }
  if (criteria == "PC1") {
    CT[1, ] <- log(NT/NT1) * ii * NT1/NT
  }
  
  if (criteria == "PC2") {
    CT[1, ] <- (NT1/NT) * log(GCT) * ii
  }
  
  if (criteria == "PC3") {
    CT[1, ] <- ii * log(GCT)/GCT
  }
  
  IC1 <- matrix(0, dim(CT)[1], I(kmax + 1))
  
  Sigma <- matrix(0, 1, I(kmax + 1))
  
  XX <- tcrossprod(x)
  
  eig <- svd(t(XX))
  
  Fhat0 <- eig$u
  
  eigval <- as.matrix(eig$d)
  
  Fhat1 <- eig$v
  
  sumeigval <- apply(eigval, 2, cumsum)/sum(eigval)
  
  Sigma[kmax + 1] <- sum(x * x)/NT
  
  IC1[, kmax + 1] <- Sigma[kmax + 1]
  
  if (criteria != "eigen") {
    for (i in kmax:1) {
      
      Fhat <- Fhat0[, 1:i]
      
      lambda <- crossprod(Fhat, x)
      
      chat <- Fhat %*% lambda
      
      ehat = x - chat
      
      Sigma[i] <- sum(ehat * ehat)/NT
      
      IC1[, i] <- (Sigma[i]) + Sigma[kmax + 1]*CT[, i]
    }
    if (criteria == "eigen") {
      
      for (j in 1:I(nrow(sumeigval))) {
        
        if (sumeigval[j] >= 0.5) {
          ic1 <- j
          break
        }
      }
      
    }
    
    ic1 <- which(IC1==min(IC1))
    
    ic1 <- ifelse(ic1 <= kmax, ic1 * 1, ic1 * 0)
    
  }
  
  if (ic1 == 0) {
    
    Fhat = matrix(0, T, N)
    
    lambda = matrix(0, N, T)
    
    chat = matrix(0, T, N)
  } else {
    
    Fhat <- Fhat0[, 1:ic1]
    
    lambda <- crossprod(x, Fhat)
    
    chat <- Fhat %*% t(lambda)
  }
  
  
  output <- list(ic = ic1, lambda = lambda, Fhat = Fhat)
  
  return(output)
} 

#'@title Determining The Number of Factors Using Yao and Lam's Method
#'
#'@description This function approximates the number of factors in a factor model
#'  for large N by T matrices using the methods  found in Lam and Yao(2012).
#'
#'@usage nof(x,k0)
#'
#'@param x A matrix containing the data.
#'
#'@param k0 An integer represents the number of lags to sum together. 
#'
#' @return r An integer of the approximate number of factors.
#'
#' @references Lam C , Yao Q . Factor modeling for high-dimensional time series: 
#' Inference for the number of factors[J]. 
#' LSE Research Online Documents on Economics, 2012, 40(40):694-726.
#'@export
nof<-function(x,k0){
  L<-dim(x)[1]
  N<-dim(x)[2]
  s<-matrix(0,N,N)
  for (q in 1:k0) {
    s.q<-t(x[(q+1):L,])%*%x[1:(L-q),]/(N-q)
    s<-s+s.q%*%t(s.q)
  }
  e.value<-eigen(s)$values[1:(L*2/3)]
  er<-e.value[-1]/e.value[-length(e.value)]
  r<-which(er==min(er))
  return(r)
}


#' @title dt.d
#' @name dt.d
#' @description A dataset used to illustrate the performance of \code{getnfac} and \code{nof}.
#' @examples
#' \dontrun{
#' data(dt.d)
#' r<-nof(dt.d,5)
#' print(r)
#' }
NULL

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions \code{rw.R}) and Cpp functions (\code{rwRcpp}.
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats runif
#' @useDynLib StatComp20009
NULL