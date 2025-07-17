/*===========================================================================*/
/* eRah C functions -                                                        */
/* Distributed under GNU General Public License version 3                    */
/*===========================================================================*/

#include <Rcpp.h>
#include <string>
#include <cmath>

//For the RUNMIN/RUNMAX functions:
#include <vector>
#include <cfloat>
#include <Rmath.h>  // Include this for R_FINITE


//#include <boost/lexical_cast.hpp>   

using namespace Rcpp;

//using boost::lexical_cast;

//How to compile a Rcpp function:
//Execute: compileAttributes(pckg_dir), from R

double correlation_dot_product(NumericVector x, NumericVector y) {
        int n = x.size();
    if (n != y.size() || n == 0) {
        stop("Vectors must be of the same size and non-empty");
    }

    // Compute dot product
    double dot_product = sum(x * y);

    // Compute the norms (magnitudes) of the vectors
    double norm_x = sqrt(sum(x * x));
    double norm_y = sqrt(sum(y * y));

    // Compute and return the cosine similarity (Pearson correlation)
    return (dot_product) / (norm_x * norm_y);
}

// [[Rcpp::export]]
NumericMatrix getRTSpectMatrix(NumericMatrix mat, NumericVector rt_vect,  NumericVector class_vect, double rt_threshold, double corr_threshold) {
    int p = mat.ncol();  // Number of columns (variables)
    //int n = mat.nrow();   // Number of rows (observations)
    std::vector<std::vector<double>> results;  // Store indices and distance above threshold

    // Compute pairwise correlations
    for (int i = 0; i < p; ++i) {
        for (int j = i + 1; j < p; ++j) {  // Avoid diagonal and redundant pairs
            if(class_vect[i]!=class_vect[j]){
                double rtdistance = std::abs(rt_vect[i] - rt_vect[j]);  // Absolute distance between rt_vect[i] and rt_vect[j]
                if(rtdistance <= rt_threshold){
                    double corr = correlation_dot_product(mat(_, i), mat(_, j));
                    if (corr >= corr_threshold) {  // Store only if above threshold
                        results.push_back({(double)i + 1, (double)j + 1, corr, rtdistance});  // Store 1-based indices and distance
                    }
                } 
            }     
        }
    }
    // Convert vector to NumericMatrix for output
    int nx = results.size();
    NumericMatrix output(nx, 4);  
    for (int k = 0; k < nx; ++k) {
        output(k, 0) = results[k][0];  // Variable 1 index
        output(k, 1) = results[k][1];  // Variable 2 index
        output(k, 2) = results[k][2];  // Dot-product similarity
        output(k, 3) = results[k][3];  // RT absolute distance
    }

    return output;
}

/*===========================================================================*/
/* runfunc - running window functions                                        */
/* Adapted (originally in C) to C++ from Jarek Tuszynski functions           */
/* in caTools (archived on CRAN)                                             */
/*===========================================================================*/

inline double SumErr(double a, double b, double ab) {
    return (((a > b) == (a > -b)) ? (b - (ab - a)) : (a - (ab - b)));
}

inline void SUM_1(double x, int n, double &Sum, double &Err, int &Num) {
    if (R_FINITE(x)) {
        double y = Sum;
        Err += x;
        Sum += Err;
        Num += n;
        Err = SumErr(y, Err, Sum);
    }
}

/*==================================================================================*/
/* Mean function applied to (running) window. All additions performed using         */
/* addition algorithm which tracks and corrects addition round-off errors (see      */  
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps)*/
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/

// [[Rcpp::export]]
NumericVector runmean(NumericVector In, int nWin) {
    int n = In.size();
    NumericVector Out(n, NA_REAL);
    int k2 = nWin / 2;
    double Sum = 0, Err = 0;
    int Num = 0;
    
    for (int i = 0; i < k2; i++) {
        SUM_1(In[i], 1, Sum, Err, Num);
    }
    
    for (int i = k2; i < nWin; i++) {
        SUM_1(In[i], 1, Sum, Err, Num);
        Out[i - k2] = (Num ? (Sum + Err) / Num : NA_REAL);
    }
    
    for (int i = nWin; i < n; i++) {
        SUM_1(In[i], 1, Sum, Err, Num);
        SUM_1(-In[i - nWin], -1, Sum, Err, Num);
        Out[i - k2] = (Num ? (Sum + Err) / Num : NA_REAL);
    }
    
    for (int i = n - k2; i < n; i++) {
        SUM_1(-In[i - nWin], -1, Sum, Err, Num);
        Out[i] = (Num ? (Sum + Err) / Num : NA_REAL);
    }
    
    return Out;
}
/*==================================================================*/
/* minimum function applied to moving (running) window              */ 
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/

// [[Rcpp::export]]
NumericVector runmin(NumericVector In, int nWin) {
    int n = In.size();
    NumericVector Out(n, NA_REAL);
    int k2 = nWin / 2;
    
    double Min = DBL_MAX;
    for (int i = 0; i < k2; i++) {
        Min = std::min(Min, In[i]);
    }
    
    for (int i = k2; i < nWin - 1; i++) {
        Min = std::min(Min, In[i]);
        Out[i - k2] = (Min == DBL_MAX ? NA_REAL : Min);
    }
    
    double ptOut = DBL_MAX;
    for (int i = nWin - 1; i < n; i++) {
        if (ptOut == Min) {
            Min = DBL_MAX;
            for (int j = i - nWin + 1; j <= i; j++) {
                Min = std::min(Min, In[j]);
            }
        } else {
            Min = std::min(Min, In[i]);
        }
        ptOut = In[i - nWin + 1];
        Out[i - k2] = (Min == DBL_MAX ? NA_REAL : Min);
    }
    
    for (int i = n - k2; i < n; i++) {
        if (ptOut == Min) {
            Min = DBL_MAX;
            for (int j = i - nWin + 1; j < n; j++) {
                Min = std::min(Min, In[j]);
            }
        }
        ptOut = In[i - nWin + 1];
        Out[i] = (Min == DBL_MAX ? NA_REAL : Min);
    }
    
    return Out;
}





