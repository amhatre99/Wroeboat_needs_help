// to enable parallel computing with "#pragma omp"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <iostream>

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// This function accepts a cylinder parameterization and a set of points; it returns a boolean vector with the points inside as true
// [[Rcpp::export]]
vec get_points_inside(
    vec indiv,
    mat p,
    uword cyl_type // cyl_type is 0 for "straight cut" cylinder (nothing beyond the end discs is included) and 1 for rounded off
)
{
  vec center1(3);
  center1[0] = indiv[0] - indiv[5]/2*sin(indiv[3])*cos(-datum::pi/2+indiv[4]);
  center1[1] =   indiv[1] - indiv[5]/2*cos(indiv[3])*cos(-datum::pi/2+indiv[4]);
  center1[2] =   indiv[2] - sin(-datum::pi/2+indiv[4])*indiv[5]/2;
  vec center2(3);
  center2[0] = indiv[0] + indiv[5]/2*sin(indiv[3])*cos(-datum::pi/2+indiv[4]);
  center2[1] = indiv[1] + indiv[5]/2*cos(indiv[3])*cos(-datum::pi/2+indiv[4]);
  center2[2] = indiv[2] + sin(-datum::pi/2+indiv[4])*indiv[5]/2;
  
  vec u(3);
  u[0] = center2[0] - center1[0];
  u[1] = center2[1] - center1[1];
  u[2] = center2[2] - center1[2];
  u = normalise(u);
  
  vec xds = (p.col(0) - indiv[0])*u[0] + // UNSAFE
    (p.col(1) - indiv[1])*u[1] + // UNSAFE
    (p.col(2) - indiv[2])*u[2]; // UNSAFE
  
  // HERE: we can optimize a little bit by considering only points potentially belonging to the cylinder
  
  mat pros(p.n_rows, 3);
  pros.col(0) = xds*u[0]+indiv[0];
  pros.col(1) = xds*u[1]+indiv[1];
  pros.col(2) = xds*u[2]+indiv[2];
  
  vec low_z, high_z;
  
  uword diff_dim = as_scalar(find(abs(center1 - center2) > 2.0*datum::eps, 1)); // CORRECTED A BUG HERE TO AVOID THE CASE Z[1] == Z[3]
  
  if (center1[diff_dim] > center2[diff_dim]) {
    low_z = center2;
    high_z = center1;
  } else {
    low_z = center1;
    high_z = center2;
  }
  if (cyl_type == 0) {
    // open-ended cylinder without anything protruding
    pros.rows(find(pros.col(diff_dim)<low_z[diff_dim], 0)).fill(datum::nan); // UNSAFE
    pros.rows(find(pros.col(diff_dim)>high_z[diff_dim], 0)).fill(datum::nan); // UNSAFE
  } else if (cyl_type == 1) {
    // rounded off cylinder
    uvec lows = find(pros.col(diff_dim)<low_z[diff_dim], 0); // UNSAFE
    uvec highs = find(pros.col(diff_dim)>high_z[diff_dim], 0); // UNSAFE
    for (uword iii=0; iii<lows.n_elem; iii++) {
      pros.row(lows(iii)) = low_z;
    }
    for (uword iii=0; iii<highs.n_elem; iii++) {
      pros.row(highs(iii)) = high_z;
    }
  }
  
  vec yds = sqrt(sum(square(pros-p),1));
  
  vec belongs_to_cyl(p.n_rows);
  belongs_to_cyl.fill(0);
  belongs_to_cyl(find(yds <= indiv[6])).fill(1); // thresh width
  // belongs_to_cyl = find(abs(yds-indiv[6]) < thresh, 0); // 2*thresh width
  // belongs_to_cyl = find((yds <= indiv[6]) && (yds >= indiv[6] - thresh), 0); // only outer points
  
  // return (sqrt(sum(square(center2-center1)))); // return the length of the cylinder
  return(belongs_to_cyl);
}
