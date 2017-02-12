#ifndef calcorr_hpp
#define calcorr_hpp
#include <armadillo>
#include <stdio.h>
using namespace arma;
using namespace std;


arma::fmat cal_blockcor(Mat<unsigned>& X);

//to standarize the vector
void vec_sub(vec& v0, double m_value, double squaresum);

//intervals between adjacent  chromsomes
Col<int> get_interval(arma::Col<int> index, uword i);


#endif /* correlationcal_hpp */
