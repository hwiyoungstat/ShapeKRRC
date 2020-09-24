#' Naive Ridge Regression Classifier
#' 
#' @param sub_X Training data for building subspaces
#' @param Test_X Test data
#' @param Test_Y Label of Test set
#' @param labmda Regularization parameter for Ridge Reg
#' 
#' @return Est_RRC_label Estimated label
#' @export
RRC1 <- function(sub_X,Test_X,Test_Y,lambda){
  
  G       <- length(sub_X)
  K       <- dim(sub_X[[1]])[1]
  N_sub   <- dim(sub_X[[1]])[2]
  N_Test  <- length(Test_Y)
  
  Hat     <- array(0,c(K,K,G))
  
  for(i in 1:G){
    X         <- sub_X[[i]]
    HX        <- Conj(t(X))
    Hat[,,i]  <- X%*%solve(HX%*%X + diag(lambda,N_sub))%*%HX 
  }
  
  Est_RRC_Mat <- matrix(0,N_Test,G)
  
  for(i in 1:G){
    Estimated_Test  <- Hat[,,i]%*%Test_X
    Est_RRC_Mat[,i] <- colSums(Mod(Test_X-Estimated_Test)^2)
  }
  
  Est_RRC_label <- apply(Est_RRC_Mat,1,which.min)
  Miss_rate  <- sum(Est_RRC_label==Test_Y)/N_Test
  
  return(Est_RRC_label)
}  


