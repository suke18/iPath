// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix CumSum(NumericMatrix m) {
    int ncol = m.ncol();
    int nrow = m.nrow();
    for (int j = 0; j < ncol; ++j) {
        for (int i = 1; i < nrow; ++i) {
            m(i, j) += m(i - 1, j);
        }
    }
    return m;
}


// [[Rcpp::export]]
double caliES2(const std::vector<double> &ranks, const std::vector<int> &pos) {
    // pos is the integer vector indicating the matched postions with the pathway
    int Nh = pos.size();
    double iES = 0.0;

    if (Nh == 0){
        iES = -1.0;
    }else{
        double NS = 0.0;
        int N = ranks.size();
        double Nm = N - Nh;
        double Norm_Nontag = 1/ Nm;

        std::vector<int> tag(N, 0);
        std::vector<int> nonTag(N, 1);

        // normalized number for sum of tagged hits
        for (int j = 0; j < Nh; j++){
            int tmpPos = pos[j];
            tag[tmpPos] = 1;
            nonTag[tmpPos] = 0;
            NS += ranks[tmpPos];
        }
        double Norm_tag = 1/ NS;

        // norm tag = 1/NS
        // norm no tag = 1/Nm
        vector<double> interCumVec(N, 0);
        vector<double> CumVec(N);

        for (int i = 0; i < N; i++){
            interCumVec[i] = tag[i]*ranks[i]*Norm_tag - nonTag[i] * Norm_Nontag;
            // cummulative sum
            for (int j = 0; j <= i; j++){
                CumVec[i] += interCumVec[j];
            }
        }

        // determine the ES
        double maxES = *max_element(CumVec.begin(), CumVec.end());
        double minES = *min_element(CumVec.begin(), CumVec.end());

        if (maxES > -minES){
            iES = maxES;
        }else{
            iES = minES;
        }
    }
    return iES;
}



