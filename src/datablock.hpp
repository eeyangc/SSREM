#ifndef datablock_hpp
#define datablock_hpp

#include <stdio.h>
#include <string>
#include <armadillo>

#include <boost/algorithm/string.hpp>

#include <boost/algorithm/string.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>

using namespace arma;
using namespace std;

class DataBlock{
public:
  int chrom;
  int block_no;
  int start;
  int end;
  DataBlock(int chrom, int block_no);
  DataBlock();
  DataBlock(const DataBlock& block);
  fmat corr;
  double heritability_beta;
  double heritability;
private:
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, unsigned) {
    ar & chrom;
    ar & block_no;
    ar & start;
    ar & end;
  }
};


class Annos{
public:
    vector<string> snp_names;
    Mat<int> annos;
};

#endif /* DataBlock_hpp */
