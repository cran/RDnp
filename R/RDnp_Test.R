#' @title Robust Test for Independence in High-Dimensional Data  
#'
#' @description
#' A Robust Test Statistic for Independence in High-Dimensional Datasets 
#' @details
#' \code{RDnp_Test} function tests the complete independence in high-dimensional
#' data sets without being affected by outliers. 
#' 
#' @importFrom cellWise estLocScale wrap
#' @importFrom MASS mvrnorm
#' @importFrom stats cor pnorm
#' @param X the data. It must be matrix.
#' @param alpha numeric parameter. It gives the rate of uncontaminated observations.
#'              Allowed values are between 0.5 and 1 and the default is 0.75.  
#' @export
#'
#' @return a list with 2 elements:
#' \item{TestValue}{The value of test statistic}
#' \item{pval}{The p value}
#' \item{robust}{Logical. Indicates whether the results are based on robust
#'  statistic. Here, it returns \code{robust=TRUE}}
#' @references Bulut, H (Unpublished). A Robust Test Statistic for Independence in High Dimensional Data
#' @author Hasan BULUT <hasan.bulut@omu.edu.tr>
#' @examples
#' 
#' # Under H0
#' library(MASS)
#' data_H0<-mvrnorm(n = 20,mu = rep(0,30),Sigma = diag(30))
#' RDnp_Test(data_H0)
#' 
#' # Under H1
#' library(MASS)
#' data_H1<-mvrnorm(n = 20,mu = rep(0,30),Sigma = (diag(30)+1))
#' RDnp_Test(data_H1)


RDnp_Test<-function(X,alpha=0.75){
  n=nrow(X)
  p=ncol(X)
  
  estX    <- cellWise::estLocScale(X,alpha=alpha)
  Xw      <- cellWise::wrap(X, estX$loc, estX$scale)$Xw
  Rw=cor(Xw)
  
  m=ncol(Rw)    
  rDnp=0
  for (i in 2:m){
    k=i-1
    for(j in 1:k){
      rDnp=rDnp+abs(Rw[i,j])
    }
  }
  
  C=(p*(p-1))/2
  GAMMA=(gamma(1)*gamma((n-1)/2))/(sqrt(pi)*gamma(n/2))
  MUnp=C*GAMMA
  TAUKarenp=C*(1/(n-1)-GAMMA^2)
  
  rDnp_star=(rDnp-MUnp)/sqrt(TAUKarenp)
  p_value=pnorm(rDnp_star,lower.tail=FALSE)
  
  return(list(TestValue=rDnp_star, pval=p_value,robust=TRUE))
}