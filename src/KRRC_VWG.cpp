#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
double VWdist(arma::cx_mat x,arma::cx_mat y){
  arma::cx_mat A = x-y;
  double res = trace(real(A.t() * A));
  return(res);
}




// [[Rcpp::export]]
double EVWG(arma::cx_mat x, arma::cx_mat y, double Sig){
  double res = exp(-VWdist(x,y)/Sig);
  return(res);
}







//' Kernel Ridge Regreggion Classifier with Veronese Whitney Gaussian
//'
//' @param x Training Array for building subspaces
//' @param xtest Test set (Array)
//' @param y Test label
//' @param Sig the regularization parameter for kernel
//' @param lambda the regularization parameter for Ridge Reg
//' @return Est_Mat distance matrix
//' @export
// [[Rcpp::export]]
arma::mat  KRRC_VWG(Rcpp::List x, arma::cx_cube xtest, arma::vec y, double Sig, double lambda){
  int G = x.size();
  int N_Test; N_Test = xtest.n_slices;  
  arma::vec N_Sub(G);
  Rcpp::List Kmat(G);
  Rcpp::List Kvec(G);
  
  arma::mat Est_Mat(N_Test, G, arma::fill::zeros);
  
  for(int i=0; i<G; i++){
    
    arma::cx_cube sub_temp = x[i];
    N_Sub(i) = sub_temp.n_slices;
    
    arma::mat Temp_mat(N_Sub(i),N_Sub(i)); Temp_mat.fill(0);
    arma::mat Temp_vec(N_Sub(i),N_Test); Temp_vec.fill(0);
    
    for(int j=0; j<N_Sub(i); j++){
     
      for(int k=0; k<N_Sub(i); k++){
        Temp_mat(j,k) = EVWG(sub_temp.slice(j),sub_temp.slice(k), Sig);
      }
      
      for(int l=0; l<N_Test; l++){
        Temp_vec(j,l) = EVWG(sub_temp.slice(j),xtest.slice(l), Sig);
      }
    }
    
    Kmat[i] = Temp_mat;
    Kvec[i] = Temp_vec;
    
    arma::mat I(N_Sub(i),N_Sub(i)); I.eye();
    arma::mat Lambda = lambda*I;
    
    arma::mat MAT1(N_Sub(i),N_Sub(i),arma::fill::zeros);
    MAT1 = inv(Temp_mat +Lambda);
    
    arma::mat MAT2(N_Sub(i), N_Sub(i),arma::fill::zeros); 
    MAT2 = Temp_mat + 2*Lambda;   
  

    for(int s=0; s<N_Test; s++){
      arma::mat ki(N_Sub(i),1); ki.fill(0);
      arma::mat kt(1,N_Sub(i)); kt.fill(0);
      
      ki = Temp_vec.col(s);
      kt = -ki.t();
      Est_Mat(s,i) = arma::as_scalar(kt * MAT1 * MAT2 * MAT1 * ki); 
    }
  }
  
  
  return(Est_Mat);
}