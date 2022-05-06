// to enable parallel computing with "#pragma omp"
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <iostream>

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// rewrite of the R `cut` function, returning in the 1-based array convention the result
// values lt the first break are assigned to the "0" bin
// values gte the last break are assigned to the "n+1"

// [[Rcpp::export]]
uvec cut_nocheck(vec V, vec edges)
{
  uvec r(V.n_elem);
  uvec sol;
  for (uword i=0; i<V.n_elem; i++) {
    sol = find(edges > V[i], 1);
    if (sol.n_elem == 1) {
      r(i) = as_scalar(sol);
    } else {
      r(i) = 0; // = edges.n_elem
    }
  }
  return (r);
}

// [[Rcpp::export]]
uvec cut_check(vec V, vec edges)
{
  uvec r(V.n_elem);
  r.fill(datum::nan);
  uvec sol;
  uvec within = find( (V >= edges(0)) && (V <= edges(edges.n_elem-1)));
  for (uword i=0; i<within.n_elem; i++) {
    sol = find(edges > V[within[i]], 1);
    if (sol.n_elem == 1) {
      r(within[i]) = as_scalar(sol);
    } else {
      r(within[i]) = 0; // = edges.n_elem
    }
  }
  return (r);
}
// cut_nocheck(c(1:10), c(3:6)+0.5)
// cut_check(c(1:10), c(3:6)+0.5)


class  Scene {
  
public:
  // constructor with all global variables
  Scene(mat points_cpp_, double thresh_, uword maxPopSize) : 
  points_cpp(points_cpp_), thresh(thresh_)
  {
    varNo = 7; // the number of parameters in a solution
    hidDim = 4; // the number of "hidden dimensions"
    occupancy_cpp = umat(points_cpp_.n_rows, maxPopSize); // allocate the matrix containing the best solution per point (varNo params) with its score (1 param)
    occupancy_cpp.fill(0);
    stored_occupancy_cpp = umat(points_cpp_.n_rows, maxPopSize); // the stored copy allows to avoid re-computing occupancy for already seen solutions
    stored_occupancy_cpp.fill(0);
    surf_total = uvec(maxPopSize);
    surf_total.fill(datum::nan);
    scores = vec(maxPopSize);
    scores.fill(0);
    sol_type = uvec(maxPopSize);
    sol_type.fill(0); // no solution is initially computed
    Realized = uvec(maxPopSize);
    solutions = mat(maxPopSize, varNo);
    solutions.fill(datum::nan);
    sol_hidDim = mat(maxPopSize, hidDim);
    sol_hidDim.fill(datum::nan);
  }
  
  void iterate()
  {
    uvec sols_to_compute = find(sol_type == 1); // get new solutions
    uvec sols_to_retrieve = find(sol_type == 2); // get already computed solutions
    uword n_sols_to_compute = sols_to_compute.n_elem;
    uword n_sols_to_retrieve = sols_to_retrieve.n_elem;
    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel for schedule(auto)
    for (uword indiv_i=0; indiv_i < n_sols_to_compute; indiv_i++) {
      // compute the occupancy
      update_score(sols_to_compute(indiv_i));
      sol_type(sols_to_compute(indiv_i)) = 2;
    } 
    #pragma omp parallel for schedule(auto)
    for (uword indiv_i=0; indiv_i < n_sols_to_retrieve; indiv_i++) {
      // just retrieve the stored copy of the occupancy
      occupancy_cpp.col(sols_to_retrieve(indiv_i)) = stored_occupancy_cpp.col(sols_to_retrieve(indiv_i));
    }
    compute_realized_span();
  }
  
  void change_indiv(rowvec indiv, uword indiv_i)
  {
    solutions.row(indiv_i) = indiv;
    sol_hidDim.row(indiv_i).fill(datum::nan);
    scores(indiv_i) = 0;
    surf_total(indiv_i) = datum::nan; 
    sol_type(indiv_i) = 1; // not already computed
  }
  
  void remove_indiv(uword indiv_i)
  {
    solutions.row(indiv_i).fill(datum::nan);
    sol_hidDim.row(indiv_i).fill(datum::nan);
    scores(indiv_i) = 0;
    surf_total(indiv_i) = datum::nan;
    sol_type(indiv_i) = 0; // absent solution
  }
  
  double compute_cyl(
      const vec& indiv,
      vec& xds,
      mat& pros,
      vec& yds,
      uvec& belongs_to_cyl,
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
    
    xds = (points_cpp.col(0) - indiv[0])*u[0] + // UNSAFE
      (points_cpp.col(1) - indiv[1])*u[1] + // UNSAFE
      (points_cpp.col(2) - indiv[2])*u[2]; // UNSAFE
    
    // HERE: we can optimize a little bit by considering only points potentially belonging to the cylinder
    
    pros.resize(points_cpp.n_rows, 3);
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
    
    yds = sqrt(sum(square(pros-points_cpp),1)); // can optimize here by not squarerooting
    
    belongs_to_cyl = find(abs(yds-indiv[6]) < thresh/2.0, 0); // thresh width
    // belongs_to_cyl = find(abs(yds-indiv[6]) < thresh, 0); // 2*thresh width
    // belongs_to_cyl = find((yds <= indiv[6]) && (yds >= indiv[6] - thresh), 0); // only outer points
    
    // return (sqrt(sum(square(center2-center1)))); // return the length of the cylinder
    return(indiv[5]);
  }
  
  void update_score(uword indiv_i) {
    // declare and compute the variables describing the cylinder
    double h; vec xds; mat pros; vec yds; uvec belongs_to_cyl;
    h = compute_cyl(vectorise(solutions.row(indiv_i)), xds, pros, yds, belongs_to_cyl, 0);
    
    vec a_breaks;
    vec h_breaks;
    
    // to compute the best fit orientation only once
    bool return_orient = true;
    bool desired_orientation_computed = false;
    
    vec xds_belongs = xds(belongs_to_cyl);
    mat points_cpp_belongs = points_cpp.rows(belongs_to_cyl);
    mat pros_belongs = pros.rows(belongs_to_cyl);
    
    //double score_normals = datum::nan;
    vec angles;
    mat desired_orientation;
    
    vec orient(3);
    double proj_length;
    
    if (return_orient) {
      if (belongs_to_cyl.n_elem == 0) {
        orient.fill(datum::nan);
        proj_length = datum::nan;
      } else if (belongs_to_cyl.n_elem == 1) {
        orient = vectorise(points_cpp_belongs - pros_belongs);
        proj_length = xds_belongs(0);
      } else {
        desired_orientation = points_cpp_belongs - pros_belongs;
        desired_orientation_computed = true;
        orient = vectorise(mean(desired_orientation, 0));
        proj_length = mean(xds_belongs);
      }
      orient = normalise(orient);
      sol_hidDim(indiv_i, 0) = proj_length;
      sol_hidDim(indiv_i, 1) = orient(0);
      sol_hidDim(indiv_i, 2) = orient(1);
      sol_hidDim(indiv_i, 3) = orient(2);
    }
    
    // double surface = 2 * datum::pi * solutions(indiv_i,6) * h;
    
    if (belongs_to_cyl.n_elem <= 1) {
      occupancy_cpp.col(indiv_i).fill(datum::nan);
    } else {
      if (!desired_orientation_computed) {
        desired_orientation = points_cpp_belongs - pros_belongs;
      }
      vec normalized_vec = normalise(vectorise(desired_orientation.row(0)));
      
      mat base_orientation = repmat(normalized_vec.t(), desired_orientation.n_rows, 1);
      vec cyl_angles = acos( sum(desired_orientation % base_orientation, 1) /
                               ( sqrt(sum(desired_orientation % desired_orientation, 1)) ) );
      // ^note: would be more efficient to normalize the desired orientation right away.
      
      // prevent a bug downstream with datum::nan inserted in the vector
      uvec bad_rounding = find_nonfinite(cyl_angles);
      if (bad_rounding.n_elem > 0) {
        cyl_angles(bad_rounding).fill(0);
      }
      
      // Now orient the result: NOT OPTIMIZED
      vec center1(3); // this paragraph is an ugly copy-and-paste
      vec indiv = vectorise(solutions.row(indiv_i));
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
      
      vec side = (normalized_vec(1)*desired_orientation.unsafe_col(2) - normalized_vec(2)*desired_orientation.unsafe_col(1)) * u(0) \
        + (normalized_vec(2)*desired_orientation.unsafe_col(0) - normalized_vec(0)*desired_orientation.unsafe_col(2)) * u(1) \
        + (normalized_vec(0)*desired_orientation.unsafe_col(1) - normalized_vec(1)*desired_orientation.unsafe_col(0)) * u(2);
      
      uvec to_switch = find(side < 0);
      cyl_angles(to_switch) = 2.0*datum::pi - cyl_angles(to_switch);
      
      // Below: breaks with only full-bins (extended at the bounds)
      double n_h = ceill(h / thresh)+1; // should be uword instead of double
      double n_a = ceill(2.0*datum::pi * solutions(indiv_i,6) / thresh)+1;
      if ((n_h > 999) | (n_a > 999) | (n_h < 1) | (n_a < 1) ) {
        std::cout << "Rounding error with" << std::endl;
        std::cout << "  h=" << h << std::endl;
        std::cout << "  r=" << solutions(indiv_i,6) << std::endl;
        std::cout << "  thresh=" << thresh << std::endl;
        occupancy_cpp.col(indiv_i).fill(datum::nan);
      } else {
        // std::cout << "n_h=" << n_h << "    n_a=" << n_a << std::endl;
        h_breaks = vec(n_h);
        for (uword i=0; i<n_h; i++) {
          h_breaks(i) = -h/2.0 + ((double)i)*thresh;
        }
        h_breaks(0) -= thresh/10.0; //!/2.0
        h_breaks(h_breaks.n_elem-1) += thresh/10.0; //!/2.0
        a_breaks = vec(n_a);
        for (uword i=0; i<n_a; i++) {
          a_breaks(i) = ((double)i)*thresh/solutions(indiv_i,6);
        }
        a_breaks(0) -= thresh/solutions(indiv_i,6)/10.0; // /2.0
        a_breaks(a_breaks.n_elem-1) += thresh/solutions(indiv_i,6)/10.0; // /2.0
      
        uvec binned_h = cut_nocheck(xds_belongs, h_breaks);
        uvec binned_a = cut_nocheck(cyl_angles, a_breaks);
        
        uvec binned_surface = ((binned_h) + 1000 * (binned_a)) * (binned_h>0) * (binned_a>0) ;
        // re-initialize the column by filling with 0
        // occupancy_cpp.unsafe_col(indiv_i).fill(0); // UNSAFE
        // occupancy_cpp.unsafe_col(indiv_i)(belongs_to_cyl) = binned_surface; // QAD // UNSAFE
        uvec newcol(occupancy_cpp.n_rows);
        newcol.fill(0);
        newcol(belongs_to_cyl) = binned_surface;
        occupancy_cpp.col(indiv_i) = newcol;
        surf_total(indiv_i) = (h_breaks.n_elem-1) * (a_breaks.n_elem-1); // was (x-1)*(y-1) +1
      }
    }
    stored_occupancy_cpp.col(indiv_i) = occupancy_cpp.col(indiv_i); // store a copy
  }
  
  void compute_realized_span()
  {
    uword n_not_realized;
    Realized.fill(0);
    Realized(find(sol_type == 0)).fill(1); // we won't consider absent solutions
    Realized(find(surf_total == 0)).fill(1); // removing solutions not touching a single point
    scores(find(surf_total == 0)).fill(0); // (^)
    uvec not_realized = find(Realized == 0);
    while (!prod(Realized)) {
      n_not_realized = not_realized.n_elem;
#pragma omp parallel for schedule(auto)
      for (uword i=0; i<n_not_realized; i++) {
        // re-compute the best approximation of the realized score for all remaining solutions
        uvec uniq = find_unique(occupancy_cpp.unsafe_col(not_realized(i))); // UNSAFE
        scores(not_realized(i)) = ((double)uniq.n_elem) / surf_total(not_realized(i));
      }
      uvec best_not_realized = find(scores(not_realized) == scores(not_realized).max(), 1);
      if (best_not_realized.n_elem == 0) {
        break;
      }
      
      // updates the score of the best and mark it as realized
      uword best_nr = as_scalar(not_realized(best_not_realized));
      uvec to_update = find(occupancy_cpp.col(best_nr) != 0, 0); // UNSAFE

      Realized(best_nr) = 1;
      
      // shrink the potential span of un-realized solutions
      not_realized = find(Realized == 0);
      occupancy_cpp(to_update, not_realized).fill(0);
    }
  }
  
  // This function is intended for use for debug or for singleton population (where potential==realized)
  double compute_potential_span(uword i)
  {
    uvec uniq = find_unique(occupancy_cpp.unsafe_col(i)); // UNSAFE
    if ((sol_type(i) == 0) || (surf_total(i) == 0)) {
      return(0);
    } else {
      return(((double)uniq.n_elem) / surf_total(i));
    }
  }
  
  double tst()
  {
    vec a(5);
    a(0) = 0;
    a(1) = 1;
    a(2) = datum::nan;
    // a(2) = 1;
    
    a(3) = 200;
    a(4) = 1;
    
    // a = a(find_finite(a));
    
    return(median(a));
  }
  
  mat get_points() {return(points_cpp);}
  umat get_occupancy() {return(occupancy_cpp);}
  vec get_scores() {return(scores);}
  uvec get_st() {return(surf_total);}
  uvec get_type() {return(sol_type);}
  mat get_indivs() {return(solutions);}
  mat get_hid_indivs() {return(sol_hidDim);}
  
private:
  mat solutions;
  mat sol_hidDim;
  mat points_cpp;
  umat occupancy_cpp;
  umat stored_occupancy_cpp;
  uvec surf_total;
  uvec Realized;
  uvec sol_type;
  vec scores;
  double thresh;
  uword varNo;
  uword hidDim;
};

RCPP_MODULE(mod_scene) {
  class_<Scene>( "Scene" )
  .constructor<mat, double, int>()
  .method("update_score", &Scene::update_score)
  .method("compute_realized_span", &Scene::compute_realized_span)
  .method("compute_potential_span", &Scene::compute_potential_span)
  .method("change_indiv", &Scene::change_indiv)
  .method("remove_indiv", &Scene::remove_indiv)
  .method("iterate", &Scene::iterate)
  // .method("tst", &Scene::tst)
     .method("get_points", &Scene::get_points)
     .method("get_occupancy", &Scene::get_occupancy)
     .method("get_scores", &Scene::get_scores)
     .method("get_st", &Scene::get_st)
     .method("get_type", &Scene::get_type)
     .method("get_indivs", &Scene::get_indivs)
     .method("get_hid_indivs", &Scene::get_hid_indivs)
  ;
}
