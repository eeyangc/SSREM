
#include "ssremdata.hpp"
#include <iomanip>
#include <boost/algorithm/string.hpp>



float normalCFD(float value)
{
  return 0.5 * erfc(-value * M_SQRT1_2);
}


/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers){

  uword size = identifiers.size();
  Col<int> pos(size);
  for(uword i = 0; i < size; i++){
    pos[i] = -1;
    for(uword j = 0; j < fields.size(); j++){
      if(boost::iequals(fields[j], identifiers[i])){
        pos[i] = (int)j ;
        break;
      }
    }

  }
  return pos;

}


GenoInfo::GenoInfo(string stringname) {
  this -> genofile = stringname;
  cout << endl <<"Loading GEnotype Data......" << endl;
  string famfile = stringname;
  famfile += ".fam";
  string bimfile = stringname;
  bimfile += ".bim";

  this -> N =  getLineNum(famfile);
  this -> P =  getLineNum(bimfile);

  fvec y = read_phenotypes(famfile, N);

  unsigned* X = new unsigned[ N * P];
  clock_t t1 = clock();
  readPlink(stringname,N, P, X);
  cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
  cout <<"Sample Size =" <<  N << " SNP Number:" << P << endl;
  arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
  // arma::fmat corMat = calCorr(Xdata, index, avbIndex, bandwidth);
  this->bandwidth = 200;
  this -> X = *Xdata;
  this -> chroms = new Chroms(bimfile, P);//&chroms;
  this -> chroms -> summary();
  this -> y = y;
  this -> ngibbs = 0;
  beta_est = nullptr; //beta to be estimated
  lp_sigma2beta = nullptr; //to record the sigma2beta for each SNP
  lp_heritability_gibbs = nullptr;
  lp_sigma2beta_gibbs = nullptr;
  lp_enrich_gibbs = nullptr;
  delete[] X;

}

void adjust_allele(Chroms& chrom, int i){
  chrom.A1Array.replace((int)'T', (int)'A');
  chrom.A1Array.replace((int)'G', (int)'C');
  chrom.A2Array.replace((int)'T', (int)'A');
  chrom.A2Array.replace((int)'G', (int)'C');
}

void adjust_allele(Col<int>& wholeAllele){
  wholeAllele.replace((int)'T', (int)'A');
  wholeAllele.replace((int)'G', (int)'C');
}


int snps_overlap(vector<SNP> chrom_x_i, vector<SNP> chrom_s_i, Col<uword>& xindex, Col<uword>& sindex, int reset){

  vector<SNP> common_snp_in_x;
  vector<SNP> common_snp_in_ss;

  if(reset == 1){
    for(int i = 0; i < chrom_x_i.size(); i++){
      chrom_x_i[i].idx = i;
    }
  }

  sort(chrom_s_i.begin(), chrom_s_i.end());
  sort(chrom_x_i.begin(), chrom_x_i.end());

  set_intersection(chrom_x_i.begin(),chrom_x_i.end(),chrom_s_i.begin(),chrom_s_i.end(),back_inserter(common_snp_in_x));

  set_intersection(chrom_s_i.begin(),chrom_s_i.end(), chrom_x_i.begin(),chrom_x_i.end(),back_inserter(common_snp_in_ss));

  Mat<uword> indexPair(common_snp_in_x.size(), 2);
  for(int i = 0; i < common_snp_in_x.size(); i++){
    indexPair(i,0) = common_snp_in_x[i].idx;
    indexPair(i,1) = common_snp_in_ss[i].idx;
  }

  uvec sorted_index = sort_index(indexPair.col(0));
  indexPair = indexPair.rows(sorted_index);
  xindex = indexPair.col(0);
  sindex = indexPair.col(1);

  return (int)common_snp_in_x.size();

}

fvec read_phenotypes(string filename, int N){
  std::ifstream ifs(filename.c_str());

  std::string line;
  fvec y(N);
  string snpname;
  float pheno = 0;
  string tmp_str;
  int idx = 0;
  int phenopos = 5;
  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    for(int i = 0; i < phenopos; i++){
      iss >> tmp_str;
    }
    iss >> pheno;
    y[idx] = pheno;

  }
  ifs.close();
  return y;
}

Chroms read_snpnames(string filename, int P){
  std::ifstream ifs(filename.c_str());
  std::string line;
  Chroms chroms(P);//snps(P);
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  int BP = 0;
  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split( fields, line, boost::is_any_of(" \t *"));
    chromsome = (int)atoi(fields[0].c_str());
    snpname = fields[1];
    BP = atoi(fields[3].c_str());
    int a1 = *fields[4].c_str();
    int a2 = *fields[5].c_str();
    chroms.chromsome[i] = chromsome;
    SNP snp(snpname, (int)i, from_x, BP);
    chroms.snps.push_back(snp);
    chroms.A1Array[i] = a1;
    chroms.A2Array[i] = a2;
    i++;
  }
  ifs.close();
  return chroms;
}

int chroms_overlap(Chroms* chrom_x, Chroms* chrom_s,Col<uword>& xindex, Col<uword>& sindex){
  int total_overlap = 0;
  adjust_allele(chrom_s -> A1Array);
  adjust_allele(chrom_s -> A2Array);
  adjust_allele(chrom_x -> A1Array);
  adjust_allele(chrom_x -> A2Array);
  total_overlap = snps_overlap(chrom_x -> snps, chrom_s -> snps,  xindex, sindex);
  return total_overlap;
}

void Summary::cal_overlap(GenoInfo& genoinfo)
{
  cout << endl << "Calculating overlap between summary and genotype data....." << endl;
  Col<uword> xindex;
  Col<uword> sindex;
  chroms_overlap(genoinfo.getChroms(), this -> chroms, xindex, sindex);
  genoinfo.filter(xindex, sindex, this);
}

void GenoInfo::sse_gibbs_block(DataBlock* block, int thread_id){
  fmat* corr = &block ->  corr;
  float s = 0;
  float s2 = 0;
  float betak = 0;
  float* lp_corr = corr -> memptr();

  int P = block -> end - block -> start + 1;
  for(int k = 0 ; k < P; k++){
    int idx = k + block -> start;//lp_index[k];
    s = seh[idx];
    betak = betah[idx];
    s2 = s*s;
    double sigma2e = 1 / (1/s2 + lp_sigma2beta[idx]);
    double corr_sum = 0;
    for(int j = 0; j  < P; j++){
      int idx2 = j + block -> start;
      if(j == k) continue;
      corr_sum += lp_corr[k * P + j] * beta_est[idx2] / seh[idx2];

    }
    double mu_k = sigma2e * (betak/s2 - corr_sum / s);
    double betasample = mu_k + randn() * sqrt(sigma2e);
    beta_est[idx] = betasample;
  }
}


int num_total_blocks = 0;
void GenoInfo::SSREM_By_Thread(int thread_id){
  while(true){
    DataBlock* current_block = this -> next();
    if(current_block == nullptr){
      break;
    }
    sse_gibbs_block(current_block, thread_id);
  }
}


void GenoInfo::SSREM_Ht_By_Thread(int thread_id){
  while(true){
    DataBlock* current_block = this -> next();
    if(current_block == nullptr){
      break;
    }
    sse_gibbs_block_per_heritability(current_block);
  }

}

void sse_gibbs_update(DataBlock& datablock, double* beta_est, int* lp_regions, double* beta2sum_in_region){
  datablock.heritability_beta = 0;
  for(uword idx = datablock.start; idx <= datablock.end; idx++){
    beta2sum_in_region[lp_regions[idx]] += beta_est[idx] * beta_est[idx];
  }

}

void GenoInfo::sse_gibbs_block_per_heritability(DataBlock* block){
  fmat* corr = &block -> corr;
  float betak = 0;
  float* lp_corr = corr -> memptr();
  int P = block -> end - block -> start + 1;
  for(uword k = 0; k < P; k++){
    uword idx = k + block -> start;
    double n_sk_2 = nh[idx] * seh[idx] *  seh[idx];
    betak = betah[idx];
    double var_sum = n_sk_2 + betak * betak;
    for(uword j = 0; j < P; j++){
      uword idx2 = j + block -> start;
      double betaj = betah[idx2];
      double n_sk_2_hat = nh[idx2] * seh[idx2] * seh[idx2];
      double var_sum_hat = n_sk_2_hat + betaj * betaj;
      ht[idx] += (*(lp_corr + k * P + j)) * beta_est[idx] * beta_est[idx2] / sqrt(var_sum * var_sum_hat);
    }

  }
}

void sse_gibbs_block_merge_heritability(DataBlock& datablock,double* ht,int* lp_regions, double* lp_heritability_region){
  uword P =  datablock.end - datablock.start + 1;
  for(uword k = 0; k < P; k++){
    uword idx = k + datablock.start;
    lp_heritability_region[lp_regions[idx]] += ht[idx];
  }
}

void GenoInfo::filter(Col<uword> xindex, Col<uword> sindex, Summary* summary){
  Col<int> allele1X = chroms -> A1Array.rows(xindex);
  Col<int> allele2X = chroms -> A2Array.rows(xindex);

  Col<int> allele1S = summary -> chroms -> A1Array.rows(sindex);
  Col<int> allele2S = summary -> chroms -> A2Array.rows(sindex);

  Col<int> sum4diff = (allele1X + allele2X) - (allele1S + allele2S);

  uvec not_equal = find(sum4diff != 0);
  uvec equal = find(sum4diff == 0);
  cout << "Discard inconsistent snps in Genotype Data and Summary Data File" << endl;
  for(int i = 0; i < not_equal.size(); i++){
    uword idx = not_equal[i];
    cout <<"idx:" << idx <<" "<< chroms -> snps[xindex[idx]].name<<":" << xindex[idx] <<"  "<< (char)allele1X[idx] <<" "<< (char)allele2X[idx] <<"   "<< (*summary -> chroms).snps[sindex[idx]].name <<":" << sindex[idx] <<"  " <<
      (char)allele1S[idx] <<" "<< (char)allele2S[idx]<< endl;
  }

  xindex = xindex.rows(equal);
  sindex = sindex.rows(equal);

  allele1X = chroms -> A1Array.rows(xindex);
  allele2X = chroms -> A2Array.rows(xindex);

  allele1S = summary -> chroms -> A1Array.rows(sindex);
  allele2S = summary -> chroms -> A2Array.rows(sindex);

  X = X.cols(xindex);

  arma::Mat<float> lpSummary = summary -> lpsummary -> rows(sindex);

  Col<int> direction;
  direction.reset();

  fvec N_vec = lpSummary.col(2);

  if(summary -> type == beta_ss){
    direction = conv_to< Col<int> >::from(allele1X == allele1S);
    Col<float> beta_hat = (lpSummary).col(0);
    uvec adjustidx = find(direction == 0);
    cout <<"Adjust directions for " <<adjustidx.size() << " SNPs" << endl;
    beta_hat(adjustidx) = - beta_hat(adjustidx);
    (lpSummary).col(0) = beta_hat;
  }
  copy_summary_info(lpSummary);
  chroms -> P = (int)xindex.size();
  chroms -> chromsome = chroms -> chromsome.elem(xindex);
  chroms -> A1Array = chroms -> A1Array.elem(xindex);
  chroms -> A2Array = chroms -> A2Array.elem(xindex);

  vector<SNP> snps;
  for(int i = 0; i < chroms -> P; i++)
  {
    snps.push_back(chroms -> snps[xindex[i]]);
  }
  chroms -> snps = snps;
  chroms -> summary();

  this -> P = chroms -> P;
  this -> xindex = new uword[this -> P];
  memcpy(this -> xindex, xindex.memptr(), sizeof(uword) * xindex.size());

}

void GenoInfo::cal_blocks(string block_file){
  vector<arma::fmat> corr_blocks;
  vector<umat>  blocks = load_block_file(block_file.c_str());
  if(chroms -> index.size() == 0){
    chroms -> summary();
  }
  Col<int> pos(chroms -> index.size());
  Col<int> tmpsum = cumsum(chroms-> index);
  uword n_chrom = chroms -> index.size();
  vector<int> block_vector;
  int last_block_no = -1;
  for (uword i = 0; i < n_chrom; i++) {
    umat block_i = blocks[i];
    Col<int> u_i = get_interval(tmpsum, i);
    uword n_block = block_i.n_rows;
    for(int kk = 0; kk < u_i.size(); kk++){
      bool in = false;
      int block_no = -1;
      int k = u_i[kk];
      for(int j = 0; j < n_block; j++){
        in = (chroms -> snps[k].BP <= block_i(j,1))&&(chroms -> snps[k].BP >= block_i(j,0));
        if(in){
          block_no = j;
          break;
        }
      }
      if(block_no != -1){
        if(block_no < last_block_no){
          int last_chrom = this -> datablocks[this -> datablocks.size()-1].chrom;
          if(last_chrom == i){
            cout <<"Error Order" << endl;
          }
        }
        if(block_no != last_block_no && last_block_no != -1){
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[ 0 ];
          block.end = block_vector[block_vector.size() - 1];
          this -> datablocks.push_back(block);
          block_vector.clear();
        }
        block_vector.push_back(k);
        last_block_no = block_no;
      }
    }
  }
  if(block_vector.size() > 0){
    DataBlock block((int)(n_chrom - 1), last_block_no);
    block.start = block_vector[ 0 ];
    block.end = block_vector[block_vector.size() - 1];
    this -> datablocks.push_back(block);
  }

  chrom_vec.resize(datablocks.size());
  chrom_vec.zeros();
  for(int i = 0; i < this -> datablocks.size(); i++){
    int chrom_idx = this -> datablocks[i].chrom;
    uvec idx_u = find(this -> chroms -> chrom_list == chrom_idx + 1);
    uword idx = idx_u[0];
    chrom_vec[i] = (int)idx;
  }
}

void GenoInfo::SSREM_CalCorr_By_Thread( ){
  while(true){
    DataBlock* current_block = this -> next();
    if(current_block == nullptr){
      break;
    }
    Mat<unsigned> subX = this -> X.cols(current_block -> start,
                                        current_block -> end);
    arma::fmat  cor_sub = cal_blockcor(subX);
    current_block -> corr = cor_sub;
  }

}

vector<umat> GenoInfo::load_block_file(string block_file){

  vector<umat> blocks_mat;
  vector<int *> blocks;
  std::ifstream ifs(block_file.c_str());

  std::string line;
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  std::getline(ifs, line);

  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    i++;
  }

  ifs.close();
  ifs.open(block_file.c_str());
  std::getline(ifs, line);
  chrom_vec.resize(i);
  umat block_mat(i,2);
  i = 0;
  while(std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::split( fields, line, boost::is_any_of("\t"));
    chromsome = atoi(fields[0].substr(3).c_str());
    chrom_vec[i] = chromsome - 1;
    block_mat(i,0) = atof(fields[1].c_str());
    block_mat(i,1) = atof(fields[2].c_str());
    i++;
  }

  for(int k = 0; k < 22; k++){
    blocks_mat.push_back(block_mat.rows(find(chrom_vec == k)));
  }
  ifs.close();
  return blocks_mat;
}


int GenoInfo::loadX(){
  string filename = this -> genofile;
  string famfile = filename;
  famfile += ".fam";
  this -> N =  getLineNum(famfile);

  string bimfile = filename;
  bimfile += ".bim";
  int P0 =  getLineNum(bimfile);


  if(this -> y.size() == 0){
    y = read_phenotypes(filename + ".fam", N);
  }

  unsigned* X = new unsigned[ N * P0];
  clock_t t1 = clock();
  readPlink(filename, N, P0, X);
  cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
  cout <<"Sample Size =" <<  N << " SNP Number:" << P << endl;

  arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P0, false,false);
  if(xindex != NULL){
    cout << this -> P << endl;
    arma::Col<uword> xindex_vec(xindex, this -> P, false,false);
    this -> X = *Xdata;
    this -> X = this -> X.cols(xindex_vec);
  }
  delete[] X;

  return 1;
}


void GenoInfo::cal_blockscorr(){
  time_t t1_o = time(NULL);
  this -> init_access();
  std::thread* threads = new (std::nothrow) std::thread[n_thread];
  for(int thread_id = 0; thread_id < n_thread; thread_id++){
    threads[ thread_id ] = std::thread(&GenoInfo::SSREM_CalCorr_By_Thread,this);
  }
  for(uword k = 0; k < n_thread; k++){
    threads[k].join();
  }
  delete[] threads;
  threads = nullptr;
  cout << "Elapsed time  for calculating "<<datablocks.size() << " blocks is " << (time(NULL) - t1_o)*1.0 << endl;
  this -> X.clear();
}

void  GenoInfo::map_category(vector<string>& snp_names,const Mat<int>& annos){
  int P = (int)snp_names.size();
  if(annos.n_rows != P){
    cout << "Number of SNPs not Match" << endl;
    exit(0);
  }
  Chroms* chroms = new Chroms(P);
  int n_region = (int)annos.n_cols + 1;
  Row<int> row_i(annos.n_cols);
  int region = 0;
  for(int i = 0; i < P; i++){
    SNP snp(snp_names[i], i, from_category, 0);
    row_i = annos.row(i);
    region = 0;
    uvec none_zero = find(row_i == 1);
    uword nz = none_zero.size();
    if(nz > 0 ){
      region = (int)none_zero.at(nz - 1) + 1;
    }
    snp.region = region;
    chroms -> snps.push_back(snp);
  }

  Col<uword> xindex;
  Col<uword> sindex;
  int reset_idx = 1;
  int total_overlap = snps_overlap(this->chroms -> snps, chroms -> snps, xindex, sindex, reset_idx);
  cout <<"overlapped number is "<< total_overlap << endl;
  uword i1 = 0;
  uword i2 = 0;

  this-> n_region = (int)n_region;

  Col<int> regions(this->chroms -> P);
  regions.zeros();
  regions -= 1;

  for(int i = 0; i < xindex.size(); i++){
    i1 = xindex[i];
    i2 = sindex[i];
    region = chroms -> snps[i2].region;
    this-> chroms -> snps[i1].region = region;
    regions[i1] = region;
  }
  this -> regions = new int[this->chroms -> P];
  memcpy(this -> regions, regions.memptr(), this->chroms -> P * sizeof(int));
  cout <<"end of map category P=" << this->chroms-> P << endl;
  this -> P =  this->chroms -> P;

}


int idx = 0;
std::mutex _mtx1;

int next(){
  //   std::lock_guard<std::mutex> lock(_mtx1);
  if (idx >= 32) {
    //     std::lock_guard<std::mutex> unlock(_mtx);

    return -1;
  }
  idx++;
  //    std::lock_guard<std::mutex> unlock(_mtx);
  return (idx - 1);
}

int init_access()
{
  std::lock_guard<std::mutex> lock(_mtx1);
  idx = 0;
  return 0;
}


Annos load_category(string category_file, vector<string> identifiers, string snpidentifier, string chromidentifier, string BPidentifier){

  cout << endl << "loading category data......" << endl;
  int P = getLineNum(category_file) - 1;
  std::ifstream ifs(category_file.c_str());
  std::string line;
  std::getline(ifs, line);
  vector <string> fields;
  boost::split( fields, line, boost::is_any_of(" \t *"));


  identifiers.push_back(BPidentifier);
  identifiers.push_back(chromidentifier);
  identifiers.push_back(snpidentifier);

  Col<int> pos = getPositions(fields, identifiers);

  uword snp_index = pos[pos.size()-1];
  int chr_index = pos[pos.size()-2];
  int bp_index = pos[pos.size()-3];

  pos.set_size(pos.size() - 3);

  int chrom_id = 1;
  int BP = 0;

  int region_start = -1;
  int region_end = -1;

  if(pos.size() == 1){
    region_start = pos(pos.size() - 1);
    region_end = region_start;
  }else if(pos.size() == 2){
    region_start = pos(pos.size() - 2);
    region_end = pos(pos.size() - 1);
  }

  uword n_region = pos.size();
  if(pos.size() <= 2){
    n_region = (region_end - region_start + 1);
  }
  Mat<int>* lpRegions = new Mat<int>(P, n_region, fill::zeros);

  vector<string> snp_names;
  for(uword k = 0; k < P; k++){
    std::getline(ifs, line);
    boost::split( fields, line, boost::is_any_of(" \t *"));
    if(chr_index != -1){
      chrom_id = (int)atoi(fields[chr_index].c_str());
    }
    if(bp_index != -1){
      BP = (int)atoi(fields[bp_index].c_str());
    }

    snp_names.push_back(fields[snp_index]);

    if(region_start == -1){
      for(uword j = 0; j < pos.size(); j++){
        if(pos[j] == -1) continue;
        string value = fields[pos[j]];
        if(value != "NA"){
          int v = (int)atoi(value.c_str());
          (*lpRegions)(k,j) = v;
        }
      }
    }else{
      for(int i = region_start; i <= region_end; i++){
        string value = fields[i];
        if(value != "NA"){
          int idx = i - region_start;
          int v = (int)atoi(value.c_str());
          (*lpRegions)(k,idx) = v;
        }
      }
    }

  }

  Annos annos;
  annos.snp_names = snp_names;
  annos.annos = (*lpRegions);

  return annos;


}

void GenoInfo::map_category(string category_file, vector<string> identifiers, string snpidentifier, string chromidentifier, string BPidentifier){


  clock_t t1 = clock();
  cout << endl << "loading category data......" << endl;
  int P = getLineNum(category_file) - 1;
  std::ifstream ifs(category_file.c_str());
  std::string line;
  std::getline(ifs, line);
  vector <string> fields;
  boost::split( fields, line, boost::is_any_of(" \t *"));


  identifiers.push_back(BPidentifier);
  identifiers.push_back(chromidentifier);
  identifiers.push_back(snpidentifier);

  Col<int> pos = getPositions(fields, identifiers);

  uword snp_index = pos[pos.size()-1];
  int chr_index = pos[pos.size()-2];
  int bp_index = pos[pos.size()-3];
  Chroms* chroms = new Chroms(P);

  pos.set_size(pos.size() - 3);

  int chrom_id = 1;
  int BP = 0;

  int region_start = -1;
  int region_end = -1;

  if(pos.size() == 1){
    region_start = pos(pos.size() - 1);
    region_end = region_start;
  }else if(pos.size() == 2){
    region_start = pos(pos.size() - 2);
    region_end = pos(pos.size() - 1);
  }

  uword n_region = pos.size();
  if(pos.size() <= 2){
    n_region = (region_end - region_start + 1);
  }
  Col<int> row_k(n_region);
  Mat<int>* lpRegions = new Mat<int>(P, n_region, fill::zeros);

  int region = 0;
  for(uword k = 0; k < P; k++){
    std::getline(ifs, line);
    boost::split( fields, line, boost::is_any_of(" \t *"));
    if(chr_index != -1){
      chrom_id = (int)atoi(fields[chr_index].c_str());
    }
    if(bp_index != -1){
      BP = (int)atoi(fields[bp_index].c_str());
    }

    SNP snp(fields[snp_index], (int)k, from_category, BP);

    if(region_start == -1){
      for(uword j = 0; j < pos.size(); j++){
        if(pos[j] == -1) continue;
        string value = fields[pos[j]];
        if(value != "NA"){
          int v = (int)atoi(value.c_str());
          (*lpRegions)(k,j) = v;
          row_k[j] = v;
        }
      }
    }else{
      for(int i = region_start; i <= region_end; i++){
        string value = fields[i];
        if(value != "NA"){
          int idx = i - region_start;
          int v = (int)atoi(value.c_str());
          (*lpRegions)(k,idx) = v;
          row_k[idx] = v;
        }
      }
    }


    region = 0;
    uvec none_zero = find(row_k == 1);
    uword nz = none_zero.size();
    if(nz > 0 ){
      region = (int)none_zero.at(nz - 1) + 1;
    }
    snp.region = region;
    chroms -> snps.push_back(snp);

  }

  Col<uword> xindex;
  Col<uword> sindex;

  int reset_idx = 1;
  int total_overlap = snps_overlap(this->chroms
                                     -> snps, chroms -> snps, xindex, sindex, reset_idx);
  cout <<"overlapped number is "<< total_overlap << endl;
  uword i1 = 0;
  uword i2 = 0;

  this-> n_region = (int)n_region + 1;

  Col<int> regions(this->chroms -> P);               //to be serialized

  regions.zeros();
  regions -= 1;

  for(int i = 0; i < xindex.size(); i++){
    i1 = xindex[i];
    i2 = sindex[i];
    region = chroms -> snps[i2].region;
    this->chroms -> snps[i1].region = region;
    regions[i1] = region;
  }

  cout <<"map category P=" << P <<"this->chroms.P="<<this->chroms -> P << endl;
  this -> regions = new int[this->chroms -> P];
  memcpy(this -> regions, regions.memptr(), this->chroms -> P * sizeof(int));

  cout << "Elapsed Time for Category Data is " <<(clock() - t1)*1.0/CLOCKS_PER_SEC << endl;


}

void GenoInfo::sample_beta_by_region(double* lp_beta2sum_region, int P){
  if(sigma2beta_by_region == nullptr)
    sigma2beta_by_region = new double[n_region];
  Col<int> tmp_region(this -> regions, this -> P, false);
  for(uword ir = 0; ir < this -> n_region; ir ++){
    if(a_beta_region[ir] <= 0) continue;
    double b_beta = 2.0 / lp_beta2sum_region[ir];
    double sigma2beta = as_scalar(randg( 1, distr_param(a_beta_region[ir] / 2, b_beta)));
    uvec region_ir = find(tmp_region == (int)ir);
    for(int i = 0; i < region_ir.size(); i++){
      this -> lp_sigma2beta[region_ir[i]] = sigma2beta;
    }
    sigma2beta_by_region[ir] = sigma2beta;
  }
}

double GenoInfo::gibbs_update_by_threads(std::thread* threads){
  time_t start_t = time(NULL);
  for(int thread_id = 0; thread_id < n_thread; thread_id++){
    threads[ thread_id ] = std::thread(&GenoInfo::SSREM_By_Thread, this, thread_id);
  }
  for(uword k = 0; k < n_thread; k++){
    threads[k].join();
  }
  time_t end_t = time(NULL);
  return (end_t - start_t)*1.0;
}

double GenoInfo::gibbs_heritability_by_threads(std::thread* threads){
  time_t start_t = time(NULL);
  for(int thread_id = 0; thread_id < n_thread; thread_id++){
    threads[ thread_id ] = std::thread(&GenoInfo::SSREM_Ht_By_Thread, this, thread_id);
  }

  for(uword thread_id = 0; thread_id < n_thread; thread_id++){
    threads[ thread_id ].join();
  }
  time_t end_t = time(NULL);
  return (end_t - start_t)*1.0;
}

void GenoInfo::resize_gibbs_results(int P0, int P){
  resize_array(this -> lp_heritability_gibbs, P0, P);
  resize_array(this -> lp_sigma2beta_gibbs, P0, P);
  resize_array(this -> lp_enrich_gibbs, P0, P);
}

GibbsOut GenoInfo::gibbs_blocks( ){
  cout << endl << "Starting Gibbs Sampling ......" << endl;
  this -> random_gen = new gsl_rng*[n_thread];
  for(int i = 0; i < n_thread; i++){
    this -> random_gen[i] =  gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(this -> random_gen[i], time(NULL));
  }
  int ngibbs_more =  numit / thin;

  n_block = this -> datablocks.size();
  initial_array(this -> beta_est, this -> P, zero);
  initial_array(this -> lp_sigma2beta, this -> P, one);
  int gibbsample = this -> ngibbs;
  time_t t1_o = time(NULL);
  int n_gibbs = this -> ngibbs + ngibbs_more;
  resize_gibbs_results(this -> ngibbs * n_region, n_gibbs * n_region);
  cal_num_in_region();
  double* lp_beta2sum_region = new double[n_region];
  cout << "number of variance components=" << n_region << endl;
  this -> ht = nullptr;
  reset_array(this -> ht, this -> P, zero);
  sigma2beta_by_region = nullptr;
  std::thread* threads = new (std::nothrow) std::thread[n_thread];
  if (threads == nullptr) {
    fprintf(stderr, "Fail to allocate thread\n");
  }
  double time4parallel = 0;
  double* lp_heritability_region = new double[n_region];
  reset_array(lp_heritability_region, this -> n_region, zero);
  for (uword iter = 0; iter < (numit + burnin); iter++) {

    this -> init_access();
    reset_array(lp_beta2sum_region, n_region, zero);

    time4parallel += gibbs_update_by_threads(threads);


    for(uword k = 0; k < n_block; k++){
      sse_gibbs_update(this -> datablocks[k], beta_est, this -> regions, lp_beta2sum_region);
    }

    sample_beta_by_region(lp_beta2sum_region, this -> P);

    if (iter >= burnin) {
      if ( (iter - burnin) % thin == 0) {
        reset_array(lp_heritability_region, this -> n_region, zero);
        reset_array(this -> ht, this -> P, zero);
        this -> init_access();
        time4parallel += gibbs_heritability_by_threads(threads);
        for(uword k = 0; k < n_block; k++){
          sse_gibbs_block_merge_heritability(this -> datablocks[k], ht, this -> regions, lp_heritability_region);
        }
        double total_heritability = set_gibbs_result(gibbsample, lp_heritability_region);
        print_result_in((int)iter, gibbsample);
        cout << " Total Heritability = " << total_heritability << endl << endl;
        gibbsample = gibbsample + 1;
      }
    }

  }
  GibbsOut gibbsout(this->beta_est, this -> regions, this-> P, lp_sigma2beta_gibbs, lp_heritability_gibbs, lp_enrich_gibbs, n_gibbs * n_region);
  gibbsout.num_in_region = a_beta_region;

  this -> ngibbs = n_gibbs;

  cout <<"Time for " <<numit + burnin << " iterations  is " << (time(NULL) - t1_o) <<" seconds by " <<n_thread <<" threads!"<< endl;
  cout <<"Time for parallel is " << time4parallel << " seconds!"<< endl;
  delete[] threads;
  delete[] this -> ht;
  delete[] lp_beta2sum_region;
  delete[] lp_heritability_region;
  threads = nullptr;
  return gibbsout;
}

void GenoInfo::print_result_in(int iter, int gibbsample){
  cout << "iter=" << iter << endl;
  for(int ir = 0; ir < n_region; ir ++){
    if(a_beta_region[ir] <= 0) continue;
    cout << ir <<":(sigma=" << lp_sigma2beta_gibbs[ir + gibbsample * n_region] <<", ht=" << lp_heritability_gibbs[ir + gibbsample * n_region]<<" enrich=" << this -> lp_enrich_gibbs[ir + gibbsample * n_region] <<" n=" << a_beta_region[ir] <<") " << endl;
  }
}


GibbsOut::GibbsOut(double* beta_est,int* regions, int P, float* sigma_beta, float* heiritability,float* enrich, int gibbs_size){
  this -> P = P;
  this -> gibbs_size = gibbs_size;
  this -> beta_est = new double[P];
  this -> regions = new int[P];
  this -> sigma_beta = new float[gibbs_size];
  this -> heiritability = new float[gibbs_size];
  this -> enrich = new float[gibbs_size];

  memcpy(this -> beta_est, beta_est, P * sizeof(double));
  memcpy(this -> regions, regions, P * sizeof(int));
  memcpy(this -> sigma_beta, sigma_beta, gibbs_size * sizeof(float));
  memcpy(this -> heiritability, heiritability, gibbs_size * sizeof(float));
  memcpy(this -> enrich, enrich, gibbs_size * sizeof(float));
}


double GenoInfo::set_gibbs_result(int gibbsample, double* lp_heritability_region){
  double total = 0;
  for(uword ir = 0; ir < n_region; ir ++)
  {
    total += lp_heritability_region[ir];
  }
  double average_ht =  total / this -> P;
  for(int ir = 0;ir < n_region; ir ++)
  {
    if(sigma2beta_by_region[ir] > 0 && a_beta_region[ir] > 0){
      this -> lp_heritability_gibbs[ir + gibbsample * n_region] = lp_heritability_region[ir];
      this -> lp_sigma2beta_gibbs[ir + gibbsample * n_region] =  1 / sigma2beta_by_region[ir];
      this -> lp_enrich_gibbs[ir + gibbsample * n_region] = (lp_heritability_region[ir] / a_beta_region[ir]) / average_ht;
    }
  }
  return total;
}

Summary::Summary(string summaryfile, vector<string> identifiers,vector<string> alleleIdentifiers, int sample_size_threshold, string snpidentifier,
                 string chromidentifier,  int type){
  this -> chrom_no = 1;
  this -> type = type;
  cout << "loading summary data...." << endl;
  this -> P = getLineNum(summaryfile) - 1;
  std::ifstream ifs(summaryfile.c_str());
  std::string line;
  std::getline(ifs, line);
  vector <string> fields;
  boost::split( fields, line, boost::is_any_of(" \t *"));

  identifiers.push_back(alleleIdentifiers[0]);
  identifiers.push_back(alleleIdentifiers[1]);
  identifiers.push_back(chromidentifier);
  identifiers.push_back(snpidentifier);


  Col<int> pos = getPositions(fields, identifiers);
  clock_t t1 = clock();
  uword snp_index = pos[pos.size()-1];
  int chr_index = pos[pos.size()-2];
  int allele1_idx = pos[pos.size()-4];
  int allele2_idx = pos[pos.size()-3];
  this -> P = this -> P - 1;
  this -> chroms = new Chroms(P);
  if(chr_index == -1){
    this -> chroms -> chrom_no = 1;
  }
  pos.set_size(pos.size()-4);
  gwasidentifiers = identifiers;
  lpsummary = new Mat<float>(P, pos.size(), fill::randn);
  int i = 0;
  float sample_N = 0;
  for(uword k = 0; k < P; k++)
  {
    std::getline(ifs, line);
    boost::split( fields, line, boost::is_any_of(" \t *"));
    char allele1 = *fields[allele1_idx].c_str();
    char allele2 = *fields[allele2_idx].c_str();
    int chrom_id = 1;
    if(chr_index != -1){
      chrom_id = (int)atoi(fields[chr_index].c_str());
    }

    //     this -> chroms -> push(fields[snp_index], chrom_id, (int)allele1, (int)allele2);
    sample_N = 0;
    int sample_idx = 2;
    string value = fields[pos[sample_idx]];
    sample_N = (float)atof(value.c_str());
    if(sample_N >= sample_size_threshold)  //for height
    {
      this -> chroms -> chromsome[i] = chrom_id;
      SNP snp(fields[snp_index], (int)i, from_ss, -1);
      this -> chroms -> snps.push_back(snp);
      this -> chroms -> A1Array[i] = (int)allele1;
      this -> chroms -> A2Array[i] = (int)allele2;
      for(uword j = 0; j < pos.size(); j++){
        if(pos[j] == -1) continue;
        string value = fields[pos[j]];
        if(value != "NA"){
          float v = (float)atof(value.c_str());
          (*lpsummary)(i,j) = v;
        }
      }
      i++;
    }

  }
  cout <<"Summary Statistics:" << i <<" snps kept" << endl;
  this -> chroms -> summary();
  this -> chroms -> chromsome.resize(i);
  this -> chroms -> A1Array.resize(i);
  this -> chroms -> A2Array.resize(i);
  lpsummary -> resize(i, pos.size());
  this -> P = i;

  if ( chr_index != -1){
    uvec indices = find_unique(this -> chroms -> chromsome);
    this -> chrom_no = (int)indices.size();
  }


  cout <<"Read summary time is " << (clock() - t1)/CLOCKS_PER_SEC << endl;
}

Summary::Summary(vector<string>& snps, vector<int>& A1, vector<int>& A2, vector<float>& b, vector<float>& se, vector<float>& N, int type){
  this -> P = snps.size();
  this -> chroms = new Chroms(P);
  lpsummary = new Mat<float>(P, 3);
  this -> type = type;
  for(uword k = 0; k < P; k++){
    SNP snp(snps[k], k, from_ss, -1);
    this -> chroms -> snps.push_back(snp);
    this -> chroms -> A1Array[k] = A1[k];
    this -> chroms -> A2Array[k] = A2[k];
    (*lpsummary)(k,0) = b[k];
    (*lpsummary)(k,1) = se[k];
    (*lpsummary)(k,2) = N[k];
  }

  this -> chrom_no = 1;

}




bool SNP::operator<(const SNP& obj) const{
  return this -> name < obj.name;
}
bool SNP::operator>(const SNP& obj) const{
  return this -> name > obj.name;
}

bool SNP::operator != (const SNP& obj) const{
  return this -> name.compare(obj.name) > 0;
}


bool SNP::operator == (const SNP& obj) const{
  return this -> name.compare(obj.name) == 0;
}

double take_sum(vector<double> values){
  double sum = 0;
  for(int i = 0; i < values.size(); i++){
    sum += values[i];
  }
  return sum;
}
