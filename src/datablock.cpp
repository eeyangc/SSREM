#include "datablock.hpp"

DataBlock::DataBlock(int chrom, int block_no){
  this -> chrom = chrom;
  this -> block_no = block_no;
  this -> heritability_beta = 0;
}

DataBlock::DataBlock(){

}

DataBlock::DataBlock(const DataBlock& block){
  this -> chrom = block.chrom;
  this -> block_no = block.block_no;
  this -> start = block.start;
  this -> end = block.end;
}
