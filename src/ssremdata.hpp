#ifndef ssremdata_hpp
#define ssremdata_hpp


#include <armadillo>
#include<iostream>
#include <stdio.h>
#include "calcorr.hpp"
#include "plinkfun.hpp"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <thread>
#include <mutex>
#include <map>
#include <random>
#include "datablock.hpp"


#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>

#include <fstream>
//#pragma warning(disable: 4244)
#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>




using namespace std;
using namespace arma;
using namespace boost;

class GibbsOut{
public:
  float* sigma_beta;
  float* heiritability;
  float* enrich;
  double* beta_est;
  int* regions;
  int gibbs_size;
  int P;
  Col<float> num_in_region;
  GibbsOut(double* beta_est,int* regions, int P, float* sigma_beta, float* heiritability,float* enrich, int gibbs_size);
  int get_num_gibbs(){
    return gibbs_size / num_in_region.size();
  }


};

class Chroms;

class SNP{
public:
  SNP(){

  }
  SNP(string name, int idx, int type, int BP){
    this -> name = name;
    this -> idx = idx;
    this -> type = type;
    this -> BP = BP;
    this -> region = -1;
  }
  SNP(const SNP& snp){
    this -> name = snp.name;
    this -> idx = snp.idx;
    this -> type = snp.type;
    this -> BP = snp.BP;
    this -> region = snp.region;
  }
  string name;
  int idx;
  int type;
  int BP;
  int region;
  bool operator<(const SNP&) const;
  bool operator>(const SNP&) const;
  bool operator != (const SNP&) const;
  bool operator == (const SNP&) const;

};


const int const_chrom_no = 22;
Chroms read_snpnames(string filename, int P);

class Chroms{
public:
  //      chromsomes[const_chrom_no];
  Col<int> chromsome;
  Col<int> chrom_list;
  Col<int> A1Array;
  Col<int> A2Array;
  int n_region;
  vector<SNP> snps;
  int chrom_no;
  int P;
  arma::Col<int> index;
  //to allocate the chromsomes for the chroms, chr_idx = -1 for one chromsome, chr_idx != -1 for const_chrom_no 22
  Chroms(int P){
    chromsome.resize(P);
    A1Array.resize(P);
    A2Array.resize(P);
    this -> P = P;
  }

  Chroms(string bimfile, int P){
    *this = read_snpnames(bimfile, P);
  }

  void summary(){
    chrom_list = arma::unique(chromsome);
    uvec indices = sort_index(chrom_list);
    chrom_list = chrom_list(indices);
    this -> chrom_no = (int)chrom_list.size();
    index.resize(this -> chrom_no);
    for(int i = 0; i < this -> chrom_no; i++ ){
      this -> index[i] = (int)sum(this -> chromsome == chrom_list[i]);
    }
  }


  void clear(){
    chromsome.reset();
    A1Array.reset();
    A2Array.reset();
    index.reset();
  }

  ~Chroms(){
    //     delete[] snps;
  }
  Chroms(){
  }

};

#define pvalue_v 0
#define zvalue_v 1

#define beta_ss 0
#define pvalue_ss 1

#define from_ss 0
#define from_x 1
#define from_category 2

static const double zero = 0.0;
static const double one = 1.0;

/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers);

double take_sum(vector<double> values);

fvec read_phenotypes(string filename, int N);

float normalCFD(float value);

class GenoInfo;

class Summary{
public:
  int P;
  int chrom_no;
  int type;
  int N_threshold;
  Chroms* chroms;
  fmat* lpsummary;
  fmat* convert_lpsummary;
  vector<string> gwasidentifiers;
  void convert(){
    clock_t t1 = clock();
    if(type == zvalue_v){
      convert_lpsummary = new Mat<float>(lpsummary ->n_rows, lpsummary -> n_cols);
      float zvalue = 0;
      for(uword i = 0; i < lpsummary ->n_rows; i++){
        for(uword j = 0; j < lpsummary -> n_cols; j++){
          zvalue = (*lpsummary).at(i,j);
          convert_lpsummary->at(i,j) = 2*(1 - normalCFD(zvalue));
        }
      }
    }
    cout << "Convert z-value to p-value, time elapsed " << (clock() - t1) / CLOCKS_PER_SEC << endl;
  }
  Summary();
  Summary(string summaryfile, vector<string> identifiers,vector<string> alleleIdentifiers,int N_threshold, string snpidentifier = "SNP",
          string chromidentifier = "CHR",  int type = beta_ss);

  Summary(vector<string>& snps, vector<int>& A1, vector<int>& A2, vector<float>& b, vector<float>& se, vector<float>& N,int type = beta_ss);
  void cal_overlap(GenoInfo& genoinfo);
  void clear(){
    if(lpsummary != NULL)
      delete lpsummary;
    if(convert_lpsummary != NULL)
      delete convert_lpsummary;
    this -> chroms -> clear();
  }

};

static std::default_random_engine generator;
static std::normal_distribution<double> distribution(0,1);

class GenoInfo{
private:
  int n_thread = 1;  //to be serialized
  int numit = 1000;  //to be serialized
  int thin = 10;     //to be serialized
  int burnin = 100;  //to be serialized
  int ngibbs = 0;
  int P = 0;  //snp number    //to be serialized
  int n_region = 1;

  float* lp_heritability_gibbs;
  float* lp_sigma2beta_gibbs;
  float* lp_enrich_gibbs;
  int* regions;         //region array for each SNP
  double* beta_est; //beta to be estimated  //to be serialized
  double* lp_sigma2beta; //to record the sigma2beta for each SNP
  float* betah;// = v_beta.memptr();
  float* seh;// = v_se.memptr();
  float* nh;// = v_N.memptr();
  double* ht;//to be serialized
  arma::Mat<unsigned> X;
  arma::fvec y;

  u64* xindex;
  string genofile;
  Chroms* chroms;
  arma::Mat<double> co_lambda;
  vector<DataBlock> datablocks;  //to be serialized
  uword n_block;
  gsl_rng** random_gen;
  int bandwidth;
  int N;  //sample size   //to be serialized
  Col<float>  a_beta_region;
  double* sigma2beta_by_region;

  void initial_random_generator();
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, unsigned) {
    ar & genofile;
    ar & n_thread;
    ar & numit;
    ar & thin;
    ar & burnin;
    ar & datablocks;
    ar & P;
    ar & n_region;
    ar & ngibbs;

    if (Archive::is_loading::value)
    {
      assert(regions == nullptr);
      regions = new int[P];
    }
    ar & boost::serialization::make_array<int>(regions, P);

    if (Archive::is_loading::value)
    {
      assert(beta_est == nullptr);
      beta_est = new double[P];
    }
    ar & boost::serialization::make_array<double>(beta_est, P);

    if (Archive::is_loading::value)
    {
      assert(lp_sigma2beta == nullptr);
      lp_sigma2beta = new double[P];
    }
    ar & boost::serialization::make_array<double>(lp_sigma2beta, P);

    if (Archive::is_loading::value)
    {
      assert(betah == nullptr);
      betah = new float[P];
    }
    ar & boost::serialization::make_array<float>(betah, P);

    if (Archive::is_loading::value)
    {
      assert(seh == nullptr);
      seh = new float[P];
    }
    ar & boost::serialization::make_array<float>(seh, P);

    if (Archive::is_loading::value)
    {
      assert(nh == nullptr);
      nh = new float[P];
    }
    ar & boost::serialization::make_array<float>(nh, P);

    if (Archive::is_loading::value)
    {
      assert(xindex == nullptr);
      xindex = new uword[P];
    }
    ar & boost::serialization::make_array<uword>(xindex, P);

    if (Archive::is_loading::value)
    {
      assert(lp_heritability_gibbs == nullptr);
      lp_heritability_gibbs = new float[ngibbs * n_region];
    }
    ar & boost::serialization::make_array<float>(lp_heritability_gibbs, ngibbs * n_region);

    if (Archive::is_loading::value)
    {
      assert(lp_sigma2beta_gibbs == nullptr);
      lp_sigma2beta_gibbs = new float[ngibbs * n_region];
    }
    ar & boost::serialization::make_array<float>(lp_sigma2beta_gibbs, ngibbs * n_region);

    if (Archive::is_loading::value)
    {
      assert(lp_enrich_gibbs == nullptr);
      lp_enrich_gibbs = new float[ngibbs * n_region];
    }
    ar & boost::serialization::make_array<float>(lp_enrich_gibbs, ngibbs * n_region);

  }



  void resize_gibbs_results(int P0, int P);
  double gibbs_update_by_threads(std::thread* threads);
  double gibbs_heritability_by_threads(std::thread* threads);
  void sample_beta_by_region(double* lp_beta2sum_region, int P);
  void print_result_in(int iter, int gibbsample);
  double set_gibbs_result(int gibbsample, double* lp_heritability_region);

  void cal_num_in_region(){
    a_beta_region.resize(this -> n_region);
    Col<int> tmp_region(this -> regions, this -> P, false);
    for(int i = 0; i < n_region; i++){
      a_beta_region[i] = sum(tmp_region == i);
    }
  }

public:

  int getP()
  {
    return this->P;
  }

  void set_thread_num(int n_thread){
    this -> n_thread = n_thread;
  }

  void set_numit(int numit){
    this -> numit = numit;
  }

  void set_thin(int thin){
    this -> thin = thin;
  }

  void set_burnin(int burnin){
    this -> burnin = burnin;
  }

  int get_ngibbs(){
    return this -> ngibbs;
  }

  double* get_beta_est(){
    return beta_est;
  }

  int* get_regions(){
    return regions;
  }


  Chroms* getChroms(){
    return this -> chroms;
  }

  void filter(Col<uword> xindex, Col<uword> sindex, Summary* summary);

  void copy_summary_info(const arma::Mat<float>& lpSummary){
    this -> betah = new float[P];
    this -> seh = new float[P];
    this -> nh = new float[P];
    memcpy(this -> betah, lpSummary.col(0).colmem, P * sizeof(float)); //vector of beta effect size
    memcpy(this -> seh, lpSummary.col(1).colmem, P * sizeof(float));  //vector of stand err
    memcpy(this -> nh, lpSummary.col(2).colmem, P * sizeof(float)); //vector of stand err

  }

  GenoInfo(){
    bandwidth = 200;
    n_thread = 1;  //to be serialized
    numit = 1000;  //to be serialized
    thin = 10;     //to be serialized
    burnin = 100;  //to be serialized
    ngibbs = 0;
    P = 0;  //snp number    //to be serialized
    n_region = 1;
    lp_heritability_gibbs = nullptr;
    lp_sigma2beta_gibbs = nullptr;
    lp_enrich_gibbs = nullptr;
    regions = nullptr;         //region array for each SNP
    beta_est = nullptr;  //beta to be estimated  //to be serialized
    lp_sigma2beta = nullptr; //to record the sigma2beta for each SNP
    betah = nullptr; // = v_beta.memptr();
    seh = nullptr;// = v_se.memptr();
    nh = nullptr;// = v_N.memptr();
    ht = nullptr;
  }

  int loadX();
  Col<int> chrom_vec;
  GenoInfo(string stringname);
  vector<umat> load_block_file(string block_file);
  ~GenoInfo(){
    if(regions != nullptr){
      delete[] regions;
    }
  }

  DataBlock* next(){
    std::lock_guard<std::mutex> lock(_mtx);
    if (_current_idx >= datablocks.size()) {
      return nullptr;
    }
    _current_idx++;
    return &datablocks[(_current_idx - 1)];
  }

  int init_access()
  {
    std::lock_guard<std::mutex> lock(_mtx);
    _current_idx = 0;
    return 0;
  }

  arma::fmat calCorr();
  vector< vector<int> > blocks_by_thread;
  void allocate_blocks();

  void cal_blockscorr();
  void cal_blocks(string block_file);
  void map_category(string category_file, vector<string> region_identifiers, string snpidentifier,
                    string chromidentifier , string BPidentifier);
  void sse_gibbs_block(DataBlock* block, int thread_id);
  void SSREM_By_Thread(int thread_id);
  void SSREM_CalCorr_By_Thread( );

  void sse_gibbs_block_per_heritability(DataBlock* block);
  void SSREM_Ht_By_Thread(int thread_id);

  void map_category(vector<string>& snps,const Mat<int>& annos);

  void setbandwidth(int bandwidth){
    this -> bandwidth = bandwidth;
  }

  void gibbs_coheritability(DataBlock& datablock, const float* beta0, const float* s0,const float* N0, double* sigma2beta_vec, double* beta_est);

  void cal_overlap(Summary summary, Col<uword>& uindex, Col<uword>& sindex);
  GibbsOut gibbs_blocks();
  std::mutex _mtx;
  int _current_idx;

};




int chroms_overlap(Chroms& chrom_x, Chroms& chrom_s,Col<uword>& xindex, Col<uword>& sindex);
int snps_overlap(vector<SNP> chrom_x_i, vector<SNP> chrom_s_i, Col<uword>& xindex, Col<uword>& sindex, int reset = 0);
void adjust_allele(Chroms& chrom, int i);
void adjust_allele(Col<int>& wholeAllele);

void sse_gibbs_block(DataBlock& datablock, const float* beta0, const float* s0,const float* N0 , double* sigma2beta_vec, double* beta_est);
void sse_gibbs_update(DataBlock& datablock, double* beta_est, int* lp_regions, double* beta2sum_in_region);
void sse_gibbs_block_per_heritability(DataBlock& datablock, const float* beta0, const float* s0,const float* N0, double* beta_est, double* ht);
void sse_gibbs_block_merge_heritability(DataBlock& datablock,double* ht,int* lp_regions, double* lp_heritability_region);
Annos load_category(string category_file, vector<string> identifiers, string snpidentifier, string chromidentifier, string BPidentifier);
void test_threads(int n_thread);
int init_access();
double sum_vec();

template<typename T>
void initial_array(T*& a, int P, T v0){
  if(a == nullptr){
    a = new T[P];
    for(int i = 0; i < P; i++){
      a[i] = v0;
    }
  }
}

template<typename T>
void reset_array(T*& a, int P, T v0){
  if(a == nullptr){
    a = new T[P];
    for(int i = 0; i < P; i++){
      a[i] = v0;
    }
  }else{
    for(int i = 0; i < P; i++){
      a[i] = v0;
    }
  }
}

template<typename T>
void resize_array(T*& a, int P0, int P){
  if(a != nullptr && P > P0){
    T* a_n = new T[P];
    for(int i = 0; i < P0; i++){
      a_n[i] = a[i];
    }
    for(int i = P0; i < P; i++){
      a_n[i] = 0;
    }
    delete[] a;
    a = a_n;
  }else if(a == nullptr){
    a = new T[P];
  }
}

#endif
