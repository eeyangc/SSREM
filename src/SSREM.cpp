#include <Rcpp.h>
#include "ssremdata.hpp"
#include <iostream>
#include<armadillo>
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(BH)]]


NumericMatrix fmat2NMatrix(arma::fmat X){
     NumericMatrix nm(X.n_rows, X.n_cols);
     for(int i = 0; i < X.n_rows; i++){
        for(int j = 0; j < X.n_cols; j++){
          nm(i,j) = X(i,j);
        }
     }
     return nm;
}

template<typename T>
NumericMatrix array2matrix(T* X, int n_row, int n_col){
  NumericMatrix nm(n_row, n_col);
  for(int i = 0; i < n_row * n_col; i++){
      nm[i] = X[i];
  }
  return nm;
}

IntegerMatrix imat2NMatrix(arma::Mat<int> X){
  IntegerMatrix nm(X.n_rows, X.n_cols);
  for(int i = 0; i < X.n_rows; i++){
    for(int j = 0; j < X.n_cols; j++){
      nm(i,j) = X(i,j);
    }
  }
  return nm;
}

// [[Rcpp::export]]
RcppExport SEXP  load_category(Rcpp::String anno_start, Rcpp::String anno_end,
                               Rcpp::String category_file,Rcpp::String SNP="SNP",Rcpp::String CHR="CHR",Rcpp::String BP="BP") {

    vector<string> category_identifiers;
    category_identifiers.push_back(anno_start);
    category_identifiers.push_back(anno_end);
    Annos ans = load_category(category_file, category_identifiers, SNP, CHR, BP);
    List ret;
    ret["snp_names"] = ans.snp_names;
    ret["annos"] = imat2NMatrix(ans.annos);
    return ret;
}

// [[Rcpp::export]]
RcppExport SEXP  load_summary(Rcpp::String summary_file, Rcpp::String beta="b",Rcpp::String sd = "SE", Rcpp::String N = "N",
                              Rcpp::String Allele1="Allele1", Rcpp::String Allele2="Allele2",
                              Rcpp::String SNP="SNP",Rcpp::String CHR="CHR",Rcpp::String BP="BP") {

  vector<string> identifiers;
  identifiers.push_back(beta);
  identifiers.push_back(sd);
  identifiers.push_back(N);
  vector<string> allele_identifiers;
  allele_identifiers.push_back(Allele1);
  allele_identifiers.push_back(Allele2);
  string snpIdentifier = SNP;

  Summary summary(summary_file, identifiers,allele_identifiers, 0, snpIdentifier);
  vector<string> snp_names;
  vector<int> A1;
  vector<int> A2;
  for(int i = 0; i < summary.P; i++){
      snp_names.push_back(summary.chroms->snps[i].name);
      A1.push_back(summary.chroms->A1Array[i]);
      A2.push_back(summary.chroms->A2Array[i]);
  }

  List ret;
  ret["snp_names"] = snp_names;
  ret["A1"] = A1;
  ret["A2"] = A2;
  ret["summary"] = fmat2NMatrix(*summary.lpsummary);

  return ret;
}

//' @title
//' SSREM_More
//' @description
//' Continue SSREM Gibbs Training after previous gibbs procedure
//' @param more number of iterations
//' @param num_thread number of threads
//'
//' @param cache file name of the cache for the gibbs results to recover the model
//'
//' @return List of gibbs results, the same as the output of SSREM function
// [[Rcpp::export]]
RcppExport SEXP SSREM_More(int more, int num_thread = -1,  Rcpp::String cache = "ssrem.cache", Rcpp::String logfile="screen"){
  streambuf* coutBuf = NULL;
  ofstream of(logfile);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  cout << "Another circle of traing for SSREM gibbs......" << endl;
    GenoInfo geno;
    std::ifstream inStream(cache, std::ios::in | std::ios::binary);
    boost::archive::binary_iarchive iar(inStream);
    iar >> (geno);
    inStream.close();
    geno.set_numit(more);
    geno.set_burnin(0);
    if(num_thread > 0){
      geno.set_thread_num(num_thread);
    }
    geno.loadX();
    geno.cal_blockscorr();
    GibbsOut gibbs = geno.gibbs_blocks();

    std::ofstream outStream(cache, std::ios::out | std::ios::binary | std::ios::trunc);
    boost::archive::binary_oarchive oar(outStream);
    oar << (geno);
    outStream.flush();
    outStream.close();

    int n_region = gibbs.num_in_region.size();
    int ngibbs = geno.get_ngibbs();
    List ret;
    ret["beta"] = array2matrix(geno.get_beta_est(), geno.getP(), 1);
    ret["regions"] = array2matrix(geno.get_regions(), geno.getP(), 1);
    ret["gibbs_enrich"] = array2matrix(gibbs.enrich, n_region, ngibbs);
    ret["gibbs_sigma2beta"] = array2matrix(gibbs.sigma_beta, n_region, ngibbs);
    ret["gibbs_heritability"] = array2matrix(gibbs.heiritability, n_region, ngibbs);
    ret["num_in_region"] = gibbs.num_in_region;

    if(logfile != "screen"){
      cout.rdbuf(coutBuf);
      cout << "Write Personal Information over..." << endl;
    }
    of.flush();
    of.close();
    return ret;
}

//' @title
//' SSREM
//' @description
//' SSREM Gibbs Training
//'
//' @param geno_file  genotype data file of plink data format
//'
//' @param block_file  block file for the block boundary of snps
//'
//' @param summaries  input of summary statistics, data frames with column names "snp,A1,A2,b,se,N"
//'
//' @param annotations input of category information, data frames with snp column and category informations
//'
//' @param beta identifier for beta, default value is b
//'
//' @param sd identifier for se, default value is se
//'
//' @param N identifier for N(sample size), default value is N
//'
//' @param Allele1 identifier for allele 1, default value is A1
//'
//' @param Allele2 identifier for allele 2, default value is A2
//'
//' @param SNP identifier for SNP, default value is SNP
//'
//' @param CHR identifier for Chromsome, default value is CHR
//'
//' @param BP identifier for Base Pair, default value is BP
//'
//' @param num_thread number of threads, default value is 1
//'
//' @param burnin number of burnin for gibbs, default value is 100
//'
//' @param numit number of iteration for gibbs after burnin, default value is 2000
//'
//' @param thin value of interval for recording gibbs results, default value is 10
//'
//' @param logfile indicate the log file name, default value is "screen" for displaying the output information on the screen
//'
//' @param cache file name of the cache for the gibbs results to recover the model
//'
//' @return List of gibbs results, refer to the end of the function for the details
// [[Rcpp::export]]
RcppExport SEXP SSREM(Rcpp::String geno_file,
                    Rcpp::String block_file, DataFrame summaries, DataFrame annotations,
                    Rcpp::String beta="b",Rcpp::String sd = "se", Rcpp::String N = "N",
                    Rcpp::String Allele1="A1", Rcpp::String Allele2="A2",
                    Rcpp::String SNP="SNP",Rcpp::String CHR="CHR",Rcpp::String BP="BP",
                    int num_thread = 1,int burnin = 100, int numit = 2000, int thin = 10,
                    Rcpp::String logfile="screen",Rcpp::String cache = "ssrem.cache") {
  streambuf* coutBuf = NULL;
  ofstream of(logfile);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }

  CharacterVector snps_from_ss = summaries[0];
  IntegerVector A1 = summaries[Allele1];
  IntegerVector A2 = summaries[Allele2];
  NumericVector b = summaries[beta];
  NumericVector se = summaries[sd];
  NumericVector n = summaries[N];



  CharacterVector snps = annotations[0];

  IntegerMatrix annos(annotations.nrows(),annotations.size() - 1);
  for(int i = 1; i < annotations.size(); i++){
    NumericVector col_i = annotations[i];
    annos.column(i-1) = col_i;
  }
  int* lp_annos = &annos[0];

  Mat<int> annos_mat(lp_annos,annos.nrow(), annos.ncol(), false);

  vector<string> identifiers;
  identifiers.push_back(beta);
  identifiers.push_back(sd);
  identifiers.push_back(N);
  vector<string> allele_identifiers;
  allele_identifiers.push_back(Allele1);
  allele_identifiers.push_back(Allele2);
  string snpIdentifier = SNP;

  vector<string> snps_vector = Rcpp::as< vector<string> >(snps_from_ss);
  vector<int> A1_vector = Rcpp::as< vector<int> >(A1);
  vector<int> A2_vector = Rcpp::as< vector<int> >(A2);
  vector<float> b_vector = Rcpp::as< vector<float> >(b);
  vector<float> se_vector = Rcpp::as< vector<float> >(se);
  vector<float> n_vector = Rcpp::as< vector<float> >(n);

  Summary summary(snps_vector, A1_vector, A2_vector, b_vector, se_vector, n_vector);

  GenoInfo geno(geno_file);

  summary.cal_overlap(geno);
  geno.cal_blocks(block_file);
  std::vector<string> xs =  Rcpp::as< vector<string> >(snps);
  geno.map_category(xs, annos_mat);

  geno.set_thread_num(num_thread);
  geno.set_numit(numit);
  geno.set_thin(thin);
  geno.set_burnin(burnin);

  geno.cal_blockscorr();

  GibbsOut gibbs = geno.gibbs_blocks();
  int n_region = gibbs.num_in_region.size();
  int ngibbs = geno.get_ngibbs();

  std::ofstream outStream(cache, std::ios::out | std::ios::binary | std::ios::trunc);
  boost::archive::binary_oarchive oar(outStream);
  oar << (geno);
  outStream.flush();
  outStream.close();

  List ret;
  ret["beta"] = array2matrix(geno.get_beta_est(), geno.getP(), 1);
  ret["regions"] = array2matrix(geno.get_regions(), geno.getP(), 1);
  ret["gibbs_enrich"] = array2matrix(gibbs.enrich, n_region, ngibbs);
  ret["gibbs_sigma2beta"] = array2matrix(gibbs.sigma_beta, n_region, ngibbs);
  ret["gibbs_heritability"] = array2matrix(gibbs.heiritability, n_region, ngibbs);
  ret["num_in_region"] = gibbs.num_in_region;


  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  return ret;

}

