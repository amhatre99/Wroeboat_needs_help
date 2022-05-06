##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

###############################
# SEARCH SPACE INITIALIZATION #
###############################

# Initialize bounding box
bbox1 = matrix(c(apply(Points, 2, min, na.rm=T),apply(Points, 2, max, na.rm=T)), ncol=3, byrow = T)
lower = c(bbox1[1,1:3], 0, 0, min_length, min_radius)
upper = c(bbox1[2,1:3], 2*pi, pi, max_length, max_radius)

#########################
# GLOBAL SEARCH OPTIONS #
#########################

# Genetic operators
varNo = 7
hidDim = 4
objDim = 1

# Diversity options
maxdist = NULL # don't use spatial diversity
# maxdist = min_length/2
randomOffsprings = 10 ## must be less than 50% of minPopSize

# Population size
minPopSize = 100
maxPopSize = 300

# Termination criterion
min_ite = 200 # constant_ite_termination will raise it
max_ite = 5000
constant_ite_termination = 200

# 3D plot options
ite_modulo_show3d = 50
max_3d_cylinders = 20
points_cutoff = 10000
show_color_model = T
show_bounding_box = T

# Histogram plot options
ite_modulo_showbp = 1

###########################
# VARIABLE INITIALIZATION #
###########################

# Store the full points cloud
full_Points = Points

# Initialization of the C++ model
require('Rcpp')
require('RcppArmadillo')
#setwd('/media/jean/ext4/wooddebris/')
# Enable OpenMP
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp('cost_fct_cpp.cpp', rebuild=F)


# Load genetic search functions: mutation/crossover/roulette, solution ordering and termination criteria
source('fn_nsga2.R')

# Perform the search on subsets
k=3
x_breaks = seq(bbox1[1,1], bbox1[2,1], by=k*max_length)
y_breaks = seq(bbox1[1,2], bbox1[2,2], by=k*max_length)
z_breaks = seq(bbox1[1,3], bbox1[2,3], by=k*max_length)
all_solutions = list()
for (x_i in seq_along(x_breaks)) {
  x = x_breaks[x_i]
  all_solutions[[x_i]] = list()
  for (y_i in seq_along(y_breaks)) {
    y = y_breaks[y_i]
    all_solutions[[x_i]][[y_i]] = list()
    for (z_i in seq_along(z_breaks)) {
      z = z_breaks[z_i]
      if (x_i<2)
        next
      if ((x_i==2)&(y_i<2))
        next
      if ((x_i==2)&(y_i==2)&(z_i==1))
        next
      
      # Initialization of the Points subset
      selected_points = which((full_Points[,1] > x-max_length/2) &
                                (full_Points[,1] < x+(k+0.5)*max_length) &
                                (full_Points[,2] > y-max_length/2) & 
                                (full_Points[,2] < y+(k+0.5)*max_length) &
                                (full_Points[,3] > z-max_length/2) &
                                (full_Points[,3] < z+(k+0.5)*max_length))
      Points = full_Points[selected_points,]
      if (nrow(Points) < 500) {
        next # abort if empty
      }
      cat('Doing x_i=', x_i, " ; y_i=", y_i, " ; z_i=", z_i, "\n Number of points: ",nrow(Points),"\n")
      lower[1:3] = c(x, y, z)
      upper[1:3] = c(x+k*max_length, y+k*max_length, z+k*max_length)
      Scene_cpp = new(Scene, Points, thresh, maxPopSize)
      
      # Initialization of the population
      Solutions = t(sapply(1:minPopSize, function(u) c(array(runif(length(lower), lower, upper)), rep(NA,objDim+hidDim))))
      for (i in 1:minPopSize) {
        Scene_cpp$change_indiv(Solutions[i, 1:7], i-1)
      }
      popSize_cpp = minPopSize
      
      # Initialize the color model
      model_light = model
      model_light$vb = model_light$vb[,selected_points]
      model_light$material$color = model_light$material$color[selected_points]
      model_light$normals = NULL
      # model_light = NULL
      
      # Perform the search
      source('main_loop.R')
      all_solutions[[x_i]][[y_i]][[z_i]] = Solutions[which(Solutions[,varNo+1] > 0),1:(varNo+objDim)]
    }
  }
}

cyl_nb = 0
display_acceptance_ratio = acceptance_ratio
display_acceptance_ratio = 0.4
displayed_sol = NULL
All_Sols = NULL
for (x_i in seq_along(x_breaks)) {
  for (y_i in seq_along(y_breaks)) {
    for (z_i in seq_along(z_breaks)) {
      for (s_i in 1:nrow(all_solutions[[x_i]][[y_i]][[z_i]]))
      {     
        try({
          All_Sols = rbind(All_Sols, all_solutions[[x_i]][[y_i]][[z_i]][s_i,1:8])
        })
      }
      for (s_i in which(all_solutions[[x_i]][[y_i]][[z_i]][,8]>display_acceptance_ratio)) {
        try({
          display_sol(all_solutions[[x_i]][[y_i]][[z_i]][s_i,1:7])
          displayed_sol = rbind(displayed_sol, all_solutions[[x_i]][[y_i]][[z_i]][s_i,1:7])
          cyl_nb = cyl_nb + 1
          cat('.')
        })
      }
    }
  }
}
plot3d(model, type='dots', alpha=0.1, add=T)



## Zooming on the area of interest, and going back from all the solutions
xmin=-5
xmax=5
ymin=0
ymax=10
zmin=13
zmax=23
selected_points = which((full_Points[,1] > xmin) &
                          (full_Points[,1] < xmax) &
                          (full_Points[,2] > ymin) & 
                          (full_Points[,2] < ymax) &
                          (full_Points[,3] > zmin) &
                          (full_Points[,3] < zmax))
Points = full_Points[selected_points,]
# Initialize the color model
model_light = model
model_light$vb = model_light$vb[,selected_points]
model_light$material$color = model_light$material$color[selected_points]
model_light$normals = NULL
display_acceptance_ratio = 0.2
plot3d(model_light, type='dots')
kept_display_sol = NULL
for (x_i in seq_along(x_breaks)) {
  for (y_i in seq_along(y_breaks)) {
    for (z_i in seq_along(z_breaks)) {
      try({
        for (s_i in which(all_solutions[[x_i]][[y_i]][[z_i]][,8]>display_acceptance_ratio)) {
          if ((all_solutions[[x_i]][[y_i]][[z_i]][s_i,1] > xmin) &
              (all_solutions[[x_i]][[y_i]][[z_i]][s_i,1] < xmax) &
              (all_solutions[[x_i]][[y_i]][[z_i]][s_i,2] > ymin) & 
              (all_solutions[[x_i]][[y_i]][[z_i]][s_i,2] < ymax) &
              (all_solutions[[x_i]][[y_i]][[z_i]][s_i,3] > zmin) &
              (all_solutions[[x_i]][[y_i]][[z_i]][s_i,3] < zmax)) {
            display_sol(all_solutions[[x_i]][[y_i]][[z_i]][s_i,1:7])
            cyl_nb = cyl_nb + 1
            cat('.')
            kept_display_sol = rbind(kept_display_sol, all_solutions[[x_i]][[y_i]][[z_i]][s_i,1:7])
          }
        }
      })
    }
  }
}
# saveRDS(kept_display_sol,'OHSU_restart_solutions.rds')
# saveRDS(Points,'OHSU_subset_points.rds')
# saveRDS(model_light,'OHSU_subset_points_color.rds')

#########################################
### WHEN I HAVE FREE CPU: RESTART THIS: #
#########################################
bbox1 = matrix(c(apply(Points, 2, min, na.rm=T),apply(Points, 2, max, na.rm=T)), ncol=3, byrow = T)
lower = c(bbox1[1,1:3], 0, 0, min_length, min_radius)
upper = c(bbox1[2,1:3], 2*pi, pi, max_length, max_radius)

Scene_cpp = new(Scene, Points, thresh, maxPopSize) # keeping the same pop size

# Initialization of the population
Solutions = t(sapply(1:minPopSize, function(u) c(array(runif(length(lower), lower, upper)), rep(NA,objDim+hidDim))))
Solutions[1:nrow(kept_display_sol), 1:7] = kept_display_sol
for (i in 1:minPopSize) {
  Scene_cpp$change_indiv(Solutions[i, 1:7], i-1)
}
popSize_cpp = minPopSize

# Perform the search with less restrictive bounds:
max_ite = 50000
constant_ite_termination = 1000
acceptance_ratio = 0.3
min_length = 1 # LAST ITERATION
source('main_loop.R')
