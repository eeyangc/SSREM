#include "calcorr.hpp"
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;
#include "time.h"

void vec_sub(vec& v0, double m_value, double squaresum){
  v0 -= m_value;
  if(squaresum > 0)
    v0 /= squaresum;
}

arma::Col<int> get_interval(arma::Col<int> index, uword i){
  Col<int> u_i;
  if(i == 0){
    u_i = linspace< Col<int> >(0, index(0) - 1, index(0));
  }else{
    uword total = index(i) - index(i-1);
    u_i = linspace< Col<int> >(index(i-1), index(i) - 1, total);
  }
  return u_i;
}



arma::fmat cal_blockcor(Mat<unsigned>& X){
  uword size1 = X.n_rows;
  uword size2 = X.n_cols;
  vec meanX(size2);
  vec sqrtsum(size2);
  fmat Xnorm = conv_to<fmat>::from(X);
  for (int i = 0; i < size2; i++) { //calculate the mean of the vector and sqrt sum
    meanX[i] = sum(Xnorm.col(i))*1.0/size1;
    Xnorm.col(i) -=  meanX[i];
    fvec v_i = Xnorm.col(i);
    fmat pd = v_i.t() * v_i;
    sqrtsum[i] = sqrt(pd.at(0));
    if(sqrtsum[i] > 0){
      Xnorm.col(i) /= sqrtsum[i];
    }
  }
  arma::fmat corr(size2, size2);
  arma::fmat eyeI(size2, size2);
  eyeI.eye();
  fmat cor_ij(1,1);
  corr.eye();
  for(int i = 0; i < size2; i++){
    for (int j = i + 1; j < size2; j++) {
      cor_ij = ((Xnorm.col(i)).t() * (Xnorm.col(j)));
      double value = cor_ij.at(0);
      corr(i,j) = value;
      corr(j,i) = value;
    }
  }

  corr *= 0.9;
  corr += 0.1*eyeI;

  return corr;

}

