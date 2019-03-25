// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//' @title TTfast
//'
//' @description
//' Compiled fast t-tests on matrices. Takes a binary matrix
//' X with zero and non-zero values, and a matrix Y of
//' continuous values. Computes the t-test on each Y column
//' using the respective X column to define the two groups.
//' If Y is a matrix with one column, that column is used
//' to test with grouping derived from every column in X.
//' This function is used in LESYMAP with a binarized X
//' matrix derived from lesioned voxels in the brain.
//'
//' @param X binary matrix of voxels (columns) for all
//' subjects (rows).
//' @param Y matrix of behavioral scores of same size as X
//' or a matrix with a single column.
//' @param computeDOF (default=true) chooses whether to compute
//' degrees of freedom. Set to false to save time during
//' permutations.
//' @param varEqual (default=true) chooses whether to compute
//' Student t-scores (true) or Welch d-scores (false). The
//' only difference is the assumption on variance which for
//' t-scores must be satisfied. This assumption is often
//' violated in some voxels, and the use of Welch
//' (varEqual=false) is recommended for more accurate results.
//'
//' @return List with two vectors:
//' \itemize{
//' \item\code{statistic} - Student T or Welch D
//' \item\code{df} - degrees of freedom
//' }
//'
//' @examples
//' set.seed(1234)
//' lesmat = matrix(rbinom(60,1,0.2), ncol=2)
//' set.seed(12345)
//' behavior = cbind( rnorm(30) )
//' set.seed(123456)
//' behavior = cbind ( behavior, rnorm(30) )
//' test = LESYMAP::TTfast(lesmat, behavior)
//' test$statistic[,1] # -2.359317  1.040766
//'
//' @author Dorian Pustina
//'
//' @export
// [[Rcpp::export]]
List TTfast(const arma::mat& X, const arma::mat& Y,
            bool computeDOF = true, bool varEqual = true) {

  // Student and Welch formulas used here can be found at
  // https://www.statsdirect.co.uk/help/parametric_methods/utt.htm

  // make sure the two matrixes have same nrows
  if (X.n_rows != Y.n_rows) {
    // correct way to throw error from Rcpp
    throw Rcpp::exception("Input matrices have different number of rows.");
  }

  bool dualMatrix = false;
  if (X.n_rows != Y.n_rows) {
    dualMatrix = true;
  }

  if (dualMatrix && X.n_cols != Y.n_cols) {
    throw Rcpp::exception("Matrix Y must be single column or same column number as X.");
  }

  // initialize y vector with first column
  colvec y = Y.col(0);



  // initialize output vectors
  int k = X.n_cols;
  int N = X.n_rows;
  vec statistic = zeros<vec>(k);
  vec df = zeros<vec>(k);

  // fill already df if possible
  // for Welch we will fill voxel by voxel
  if (computeDOF && varEqual) {
    // Student t degrees of freedom
    df.fill(N - 2.0);
  }

  // set up variables
  double varpooled = 0.0;
  double n0 = 0.0; // adding 0.0 to int -> double
  double n1 = 0.0;
  double mean0 = 0.0;
  double mean1 = 0.0;
  double var0 = 0.0;
  double var1 = 0.0;



  // loop through voxels
  for (int vox=0 ; vox < k ; ++vox) {

    // check for user interruption
    if (vox % 50000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    if (dualMatrix) {
      y = Y.col(vox);
    }

    vec thisvox = X.col(vox);

    uvec indx0 = find(thisvox == 0);
    uvec indx1 = find(thisvox != 0);

    n0 = indx0.size() + 0.0; // adding 0.0 to int -> double
    n1 = indx1.size() + 0.0;

    mean0 = mean(y.elem(indx0));
    mean1 = mean(y.elem(indx1));

    var0 = var(y.elem(indx0), 0); // type=0 means divided by n-1
    var1 = var(y.elem(indx1), 0);

    if (varEqual) {
      // using a contorted calculation to avoid problems of variable conversion
      // we take variance type 1 and multiply by n to get sum of squared error
      varpooled = ( (var(y.elem(indx0), 1) * n0) + (var(y.elem(indx1), 1) * n1) ) / (N - 2.0);
    }

    // compute t-score
    if (varEqual) {
      // regular student t-score for variance equal
      statistic[vox] = (mean0 - mean1) / sqrt( varpooled * (1.0/n0 + 1.0/n1)  );
    } else {
      // Welch d-score for variance unequal
      statistic[vox] = (mean0 - mean1) / sqrt( (var0 / n0) + (var1 / n1)  );
    }


    if (computeDOF && ! varEqual) {
      // Welch d degrees of freedom
      df[vox] = pow( (var0/n0 + var1/n1) , 2) /
        (
          ( pow( (var0/n0) , 2) / (n0 - 1.0) ) +
          ( pow( (var1/n1) , 2) / (n1 - 1.0) )
        );
    }


  } // end for loop on voxels

  // Return output list
  return List::create( Named("statistic") = statistic,
                       Named("df") = df);

}
