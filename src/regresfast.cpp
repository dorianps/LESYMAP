// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//' @title Fast linear regressions
//'
//' @description
//' Takes a matrix of voxels and a vector of behavior
//' and runs fast regressions for each voxel. Covariates
//' can be defined (i.e. age) to find the effect of each
//' voxel on behavior within the context of other predictive
//' factors.
//'
//' @param X matrix of voxlels (columns) for all
//' subjects (rows).
//' @param y vector of behavioral scores.
//' @param covariates matrix with one or more columns.
//' Must be of same length as behavior. This variable
//' should always be set, and the next argument can tell
//' if covariates should be used or not.
//' @param hascovar logical to tell whether covariates
//' should be used.
//'
//' @return List with (1) statistic, (2) n = number of
//' subjects, (3) kxfm = degrees of freedom.
//'
//' @author Dorian Pustina
//'
//' @export
// [[Rcpp::export]]
List regresfast(const arma::mat& X, const arma::colvec& y, const arma::mat& covariates, bool hascovar = false) {

  int k = X.n_cols;
  int n = X.n_rows;
  vec statistic = zeros<vec>(k);
  vec pvals = zeros<vec>(k);

  mat xmat = mat(X.n_rows, 2, fill::ones); // initialize a small matrix

  // add covariates if requested
  if (hascovar) {
    xmat = join_horiz(xmat, covariates);
  }

  int kxmat = xmat.n_cols; // number of predictors in xmat

  // compute statistic for each voxel
  for (int vox=0 ; vox < k ; ++vox) {

    // check for user interruption
    if (vox % 50000 == 0) {
      Rcpp::checkUserInterrupt();
    }

    xmat.col(0) = X.col(vox);
    colvec coef = solve(xmat, y);
    colvec resid = y - xmat*coef;
    double sig2 = as_scalar(trans(resid)*resid/(n-kxmat));

    colvec stderrest = sqrt(sig2 * diagvec( inv(trans(xmat)*xmat)) );
    colvec tval = coef / stderrest;
    statistic[vox] = tval[0];
  }

  // return the statistics
  return List::create(
    Named("statistic") = statistic,
    Named("n") = n,
    Named("kxmat") = kxmat
    );
}
