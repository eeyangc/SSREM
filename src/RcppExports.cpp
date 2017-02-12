// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// load_category
RcppExport SEXP load_category(Rcpp::String anno_start, Rcpp::String anno_end, Rcpp::String category_file, Rcpp::String SNP, Rcpp::String CHR, Rcpp::String BP);
RcppExport SEXP SSREM_load_category(SEXP anno_startSEXP, SEXP anno_endSEXP, SEXP category_fileSEXP, SEXP SNPSEXP, SEXP CHRSEXP, SEXP BPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type anno_start(anno_startSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type anno_end(anno_endSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type category_file(category_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type SNP(SNPSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type BP(BPSEXP);
    rcpp_result_gen = Rcpp::wrap(load_category(anno_start, anno_end, category_file, SNP, CHR, BP));
    return rcpp_result_gen;
END_RCPP
}
// load_summary
RcppExport SEXP load_summary(Rcpp::String summary_file, Rcpp::String beta, Rcpp::String sd, Rcpp::String N, Rcpp::String Allele1, Rcpp::String Allele2, Rcpp::String SNP, Rcpp::String CHR, Rcpp::String BP);
RcppExport SEXP SSREM_load_summary(SEXP summary_fileSEXP, SEXP betaSEXP, SEXP sdSEXP, SEXP NSEXP, SEXP Allele1SEXP, SEXP Allele2SEXP, SEXP SNPSEXP, SEXP CHRSEXP, SEXP BPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type summary_file(summary_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type Allele1(Allele1SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type Allele2(Allele2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type SNP(SNPSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type BP(BPSEXP);
    rcpp_result_gen = Rcpp::wrap(load_summary(summary_file, beta, sd, N, Allele1, Allele2, SNP, CHR, BP));
    return rcpp_result_gen;
END_RCPP
}
// SSREM_More
RcppExport SEXP SSREM_More(int more, int num_thread, Rcpp::String cache, Rcpp::String logfile);
RcppExport SEXP SSREM_SSREM_More(SEXP moreSEXP, SEXP num_threadSEXP, SEXP cacheSEXP, SEXP logfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type more(moreSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type cache(cacheSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type logfile(logfileSEXP);
    rcpp_result_gen = Rcpp::wrap(SSREM_More(more, num_thread, cache, logfile));
    return rcpp_result_gen;
END_RCPP
}
// SSREM
RcppExport SEXP SSREM(Rcpp::String geno_file, Rcpp::String block_file, DataFrame summaries, DataFrame annotations, Rcpp::String beta, Rcpp::String sd, Rcpp::String N, Rcpp::String Allele1, Rcpp::String Allele2, Rcpp::String SNP, Rcpp::String CHR, Rcpp::String BP, int num_thread, int burnin, int numit, int thin, Rcpp::String logfile, Rcpp::String cache);
RcppExport SEXP SSREM_SSREM(SEXP geno_fileSEXP, SEXP block_fileSEXP, SEXP summariesSEXP, SEXP annotationsSEXP, SEXP betaSEXP, SEXP sdSEXP, SEXP NSEXP, SEXP Allele1SEXP, SEXP Allele2SEXP, SEXP SNPSEXP, SEXP CHRSEXP, SEXP BPSEXP, SEXP num_threadSEXP, SEXP burninSEXP, SEXP numitSEXP, SEXP thinSEXP, SEXP logfileSEXP, SEXP cacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type geno_file(geno_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type block_file(block_fileSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type summaries(summariesSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type annotations(annotationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type Allele1(Allele1SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type Allele2(Allele2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type SNP(SNPSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type BP(BPSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type numit(numitSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type logfile(logfileSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type cache(cacheSEXP);
    rcpp_result_gen = Rcpp::wrap(SSREM(geno_file, block_file, summaries, annotations, beta, sd, N, Allele1, Allele2, SNP, CHR, BP, num_thread, burnin, numit, thin, logfile, cache));
    return rcpp_result_gen;
END_RCPP
}