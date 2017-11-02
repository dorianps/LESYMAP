// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//' @title Fast Brunner-Munzel tests (v2) with permutations
//'
//' @description
//' Takes a binary matrix of voxels and a vector of behavior
//' and runs Brunner-Munzel tests on each voxel.
//' This is a fast function that corrects for infinite values
//' with a similar approach as the nparcomp package.
//' It calculates p-values by running permutations of each
//' voxel and using the ratio of times the real BM score
//' exceeds the permuted BM score.
//'
//'
//' @param X binary matrix ov voxlels (columns) for all
//' subjects (rows)
//' @param y vector of behavioral scores.
//' @param computeDOF (true) chooses whether to compute degrees
//' of freedom. Set to false to save time during permutations.
//'
//' @return List with these objects:
//' - statistic - BM values
//' - dfbm - degrees of freedom
//' - pvalue - permutation-based probability value
//'
//' @author Dorian Pustina
//'
//' /@export
// [[Rcpp::export]]
List BMperm(const arma::mat& X, const arma::colvec& y,
                     int npermBM = 20000, bool computeDOF = true) {
  // initialize output vectors
  int k = X.n_cols;
  int N = X.n_rows;
  vec statistic = zeros<vec>(k);
  vec dfbm = zeros<vec>(k);
  vec pvalue = ones<vec>(k);
  vec permvec = zeros<vec>(npermBM);

  // find duplicated indices
  uvec yuniq = find_unique(y);
  uvec dupli = ones<uvec>(y.size()); // zeros<uvec>(y.size());
  dupli.elem(yuniq).zeros();
  bool hasdupli = false;
  if( dupli.max() > 0) {
    hasdupli = true;
  }

  // rank of behavior vector, computed once
  colvec r = conv_to<colvec>::from(sort_index( sort_index(y) ) + 1);
  if (hasdupli) {
    for (int v=0; v<dupli.size(); ++v) {
      if (dupli[v]) {
        uvec sameindx = find(y == y[v]);
        r.elem(sameindx).fill( mean(r.elem(sameindx)) );
      }
    }
  }

  // loop through voxels
  for (int vox=0 ; vox < k ; ++vox) {
    for (int perm=1 ; perm <= npermBM ; ++perm) {

      // check for user interruption
      if (vox % 50000 == 0) {
        Rcpp::checkUserInterrupt();
      }

      vec thisvox = X.col(vox);

      // permute the voxel
      if (perm < npermBM) {
        while (all(X.col(vox) == thisvox)) {
          thisvox = shuffle(X.col(vox));
        }
      }

      uvec indx0 = find(thisvox == 0);
      uvec indx1 = find(thisvox != 0);

      int n1 = indx0.size();
      int n2 = indx1.size();


      // ranks within each group
      colvec r1 = conv_to<colvec>::from(sort_index( sort_index( y.elem(indx0) ) ) + 1);
      colvec r2 = conv_to<colvec>::from(sort_index( sort_index( y.elem(indx1) ) ) + 1);

      if (hasdupli) {
        // compute averages for ties
        for (int v=0; v<dupli.size(); ++v) {
          if (dupli[v]) {
            uvec sameindx = find(y.elem(indx0) == y[v]);
            if (sameindx.size() > 0) {
              r1.elem(sameindx).fill( mean(r1.elem(sameindx)) );
            }
            sameindx = find(y.elem(indx1) == y[v]);
            if (sameindx.size() > 0) {
              r2.elem(sameindx).fill( mean(r2.elem(sameindx)) );
            }
          }
        }
      }

      double p = 1.0/n1*(mean(r.elem(indx1)) - (n2+1.0)/2.0);
      if (p == 0) { p = pow(10, -5); }
      if (p == 1) { p = 1-pow(10, -5); }


      double S1 = 1.0/(n1-1.0)*(sum(square(r.elem(indx0)-r1))-n1*pow(mean(r.elem(indx0)-(n1+1.0)/2.0),2));
      double S2 = 1.0/(n2-1.0)*(sum(square(r.elem(indx1)-r2))-n2*pow(mean(r.elem(indx1)-(n2+1.0)/2.0),2));

      if (S1 == 0) { S1 = 1/(4.0*n1); }
      if (S2 == 0) { S2 = 1/(4.0*n2); }

      double s = N/( (n1+0.0) * (n2+0.0))*(S1/ (n2+0.0) + S2/ (n1+0.0));

      // if last loop, its the original order, save it to statistic
      // otherwise is one of the permutations
      if (perm == npermBM) {
        statistic[vox] =  -(p-0.5)*sqrt( (N+0.0) / s ); // notice the inversion '-', to assume smaller behavior for more lesion
      } else {
        permvec[perm] = -(p-0.5)*sqrt( (N+0.0) / s );
      }

      // // enable this if needed
      // if (computeDOF) {
      //   double m1 = mean(r.elem(indx0));
      //   double m2 = mean(r.elem(indx1));
      //   double v1 = sum(square(r.elem(indx0) - r1 - m1 + (n1 + 1)/2.0))/(n1 - 1.0);
      //   double v2 = sum(square(r.elem(indx1) - r2 - m2 + (n2 + 1)/2.0))/(n2 - 1.0);
      //   dfbm[vox] = (pow(n1 * v1 + n2 * v2, 2))/((pow(n1 * v1, 2))/(n1 - 1.0) + (pow(n2 * v2, 2))/(n2 - 1.0));
      // }


    } // end perm loop

    // now we have all permutation values, compute pvalue
    // if positive calculate higher, if negative calculate lower
    double excess = 0.0;
    if (statistic[vox] > 0) {
      excess = sum(permvec >= statistic[vox]);
    } else {
      excess = sum(permvec <= statistic[vox]);
    }

    pvalue[vox] = (excess + 1.0) / (npermBM + 1.0);

  } // end voxel loop

  // Return output list
  return List::create( Named("statistic") = statistic,
                       Named("dfbm") = dfbm,
                       Named("pvalue") = pvalue);

}
