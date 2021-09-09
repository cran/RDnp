#' @title Test for Independence in High-Dimensional Data  
#'
#' @description
#' A Test Statistic for Independence in High-Dimensional Datasets 
#' @details
#' \code{Dnp_Test} function tests the complete independence in high-dimensional
#' data sets. This statistic was proposed by Najarzadeh (2021).
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats cor pnorm
#' @param X the data. It must be matrix.
#' @export
#'
#' @return a list with 2 elements:
#' \item{TestValue}{The value of test statistic}
#' \item{pval}{The p value}
#' \item{robust}{Logical. Indicates whether the results are based on robust
#'  statistic. Here, it returns \code{robust=FALSE}}
#' @references Najarzadeg, D (2021). Testing independece in high-dimensional
#'  multivariate normal data, Communication in Statistics: Theory and Methods.
#'  50 (14): 3421-3435.
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' @examples
#' 
#' # Under H0
#' library(MASS)
#' data_H0<-mvrnorm(n = 20,mu = rep(0,30),Sigma = diag(30))
#' Dnp_Test(data_H0)
#' 
#' # Under H1
#' library(MASS)
#' data_H1<-mvrnorm(n = 20,mu = rep(0,30),Sigma = (diag(30)+1))
#' Dnp_Test(data_H1)


Dnp_Test<-function(X){    
  n=nrow(X)
  p=ncol(X)
  
  R=cor(X)  
  
  Dnp=0  
  for (i in 2:p){
    k=i-1
    for(j in 1:k){
      Dnp=Dnp+abs(R[i,j])   
    }
  }
  
  C=(p*(p-1)/2)  
  GAMMA=(gamma(1)*gamma((n-1)/2))/(sqrt(pi)*gamma(n/2))
  
  MUnp=C*GAMMA   
  TAUKarenp=C*(1/(n-1)-GAMMA^2)
  
  Dnp.star=(Dnp-MUnp)/sqrt(TAUKarenp)      
  p_value=pnorm(Dnp.star,lower.tail=FALSE)
  
  return(list(TestValue=Dnp.star, pval=p_value, robust=FALSE))
}