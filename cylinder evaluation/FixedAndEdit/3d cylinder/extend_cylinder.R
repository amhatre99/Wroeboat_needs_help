##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# fix the population size
n_sol = nrow(Solutions)

# Initialization of the C++ model
require('Rcpp')
require('RcppArmadillo')
#setwd('/media/jean/ext4/wooddebris/')
# Enable OpenMP
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp('cost_fct_cpp.cpp', rebuild=F)
Scene_cpp = new(Scene, Points, thresh, n_sol*2)

# Initialization of the population as the best solutions
for (sol_i in 1:n_sol) {
  Scene_cpp$change_indiv(Solutions[sol_i, 1:7], sol_i-1)
}
Scene_cpp$iterate()

# 
# last_iteration_did_extend = F
# for (extension in 1:(round(min_length/thresh))) {
#   # side_to_change_first = sample.int(2, size=n_sol, replace=T)
#   # for (side in 1:2) {
#     scores = Scene_cpp$get_scores()
#     Extended_Solutions = Solutions
#     Extended_Solutions[,6] = Extended_Solutions[,6] + thresh*2 # NEEDS TO CHANGE WHEN SWITCHING TO SIDES
#       for (sol_i in 1:n_sol) {
#         # NEEDS TO INSERT CODE HERE WHEN SWITCHING TO SIDES
#         Scene_cpp$change_indiv(Extended_Solutions[sol_i, 1:7], sol_i-1)
#     }
#     Scene_cpp$iterate()
#     new_scores = Scene_cpp$get_scores()
#     circular_patches_number = ceiling(2*pi*Solutions[,7]/thresh)
#     axis_patches_number = ceiling((Extended_Solutions[,6]-thresh*2) / thresh) # NEEDS TO CHANGE WHEN SWITCHING TO SIDES
#     new_filled_patches = new_scores * (circular_patches_number*(axis_patches_number+2)) - # NEEDS TO CHANGE WHEN SWITCHING TO SIDES
#       scores * (circular_patches_number*axis_patches_number)
#     extension_score = new_filled_patches / circular_patches_number / 2 # NEEDS TO CHANGE WHEN SWITCHING TO SIDES
#     # to_extend = which(extension_score >= acceptance_ratio)
#     to_extend = which(extension_score >= 0)
#     Solutions[to_extend, 1:7] = Extended_Solutions[to_extend, 1:7]
#     # Initialization of the population as the best solutions (with extension)
#     for (sol_i in 1:n_sol) {
#       Scene_cpp$change_indiv(Solutions[sol_i, 1:7], sol_i-1)
#     }
#     Scene_cpp$iterate()
#     
#     # }
# }


# last_iteration_did_extend = F
# for (extension in 1:(round(min_length/thresh))) {
#   # side_to_change_first = sample.int(2, size=n_sol, replace=T)
#   for (side in c(-1,1)) {
#     scores = Scene_cpp$get_scores()
#     Extended_Solutions = Solutions
#     for (sol_i in 1:n_sol) {
#       Extended_Solutions[sol_i, 1:3] = Extended_Solutions[sol_i, 1:3] + side*get_cylinder_axis(Extended_Solutions[sol_i,]) * thresh/2
#       Extended_Solutions[sol_i, 6] = Extended_Solutions[sol_i, 6] + thresh
#       Scene_cpp$change_indiv(Extended_Solutions[sol_i, 1:7], sol_i-1)
#     }
#     Scene_cpp$iterate()
#     new_scores = Scene_cpp$get_scores()
#     circular_patches_number = ceiling(2*pi*Solutions[,7]/thresh)
#     axis_patches_number = ceiling(Solutions[,6] / thresh)
#     new_filled_patches = new_scores * (circular_patches_number*(axis_patches_number+1)) - 
#       scores * (circular_patches_number*axis_patches_number)
#     extension_score = new_filled_patches / circular_patches_number
# 
#         # to_extend = which(extension_score >= acceptance_ratio)
#     to_extend = which(extension_score >= 0)
#     # Initialization of the population as the best solutions (with extension)
#     for (sol_i in to_extend) {
#       Solutions[sol_i, 1:7] = Extended_Solutions[sol_i, 1:7]
#       Scene_cpp$change_indiv(Solutions[sol_i, 1:7], sol_i-1)
#       last_iteration_did_extend = T
#     }
#     Scene_cpp$iterate()
#   }
#   if (!last_iteration_did_extend) {
#     break
#   }
# }

if (!('extend_acceptance_ratio' %in% ls())) {
  extend_acceptance_ratio = acceptance_ratio/2
}

last_iteration_did_extend = F
for (extension in 1:(round(min_length/thresh))) {
  side_to_change_first = sample(c(-1,1), size=n_sol, replace=T)
  for (side in c(-1,1)) {
    scores = Scene_cpp$get_scores()[(n_sol+1):(2*n_sol)]
    Extended_Solutions = Solutions
    for (sol_i in 1:n_sol) {
      Extended_Solutions[sol_i, 1:3] = Extended_Solutions[sol_i, 1:3] + side_to_change_first[sol_i]*side*get_cylinder_axis(Extended_Solutions[sol_i,]) * ((thresh + Extended_Solutions[sol_i, 6])/2)
      Extended_Solutions[sol_i, 6] = thresh
      Scene_cpp$change_indiv(Extended_Solutions[sol_i, 1:7], sol_i-1+n_sol)
    }
    Scene_cpp$iterate()
    new_scores = Scene_cpp$get_scores()[(n_sol+1):(2*n_sol)]
    
    to_extend = which(new_scores >= extend_acceptance_ratio)
    # Initialization of the population as the best solutions (with extension)
    for (sol_i in to_extend) {
      Solutions[sol_i, 1:3] = Solutions[sol_i, 1:3] + side_to_change_first[sol_i]*side*get_cylinder_axis(Solutions[sol_i,]) * thresh/2
      Solutions[sol_i, 6] = Solutions[sol_i, 6] + thresh
      Scene_cpp$change_indiv(Solutions[sol_i, 1:7], sol_i-1)
      last_iteration_did_extend = T
    }
    Scene_cpp$iterate()
  }
  if (!last_iteration_did_extend) {
    break
  }
}
