##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- parent.frame(2)$ofile
print(this.dir)
if (!is.null(this.dir)) {
this.dir <- dirname(this.dir)
setwd(this.dir)
}

## load the data from Tree 1
require('rgl')
require('Rcpp')
require('RcppArmadillo')

source('toy_examples_fct.R')
source('helpers.R')
source('fn_pointcloud.R')

# model = color.read.ply('Tree1_DG_1Hz_30kkey_3ktie_high_0_025_cropped.ply')
model = readRDS('Tree1_DG_1Hz_30kkey_3ktie_high_0_025.rds')
Points = cbind(t(model$vb[1:3, ]))
thresh = 0.03*0.1406769 # 
acceptance_ratio = 0.4 # less than half cover for the wood
min_radius = 0.2*0.1406769 # the branch has radius = 0.25 roughly
max_radius = 0.3*0.1406769 #
min_length = 0.5*0.1406769 #
max_length = 1*0.1406769 #

# Check constraints
if (('precision' %in% ls() == F) || (is.null(precision))) {
  check_constraints(thresh=thresh, min_radius=min_radius, min_length=min_length, max_length=max_length, Points_xor_precision=Points)
} else {
  check_constraints(thresh=thresh, min_radius=min_radius, min_length=min_length, max_length=max_length, Points_xor_precision=precision)
}

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
# Show the bounding boxes
for (x_i in seq_along(x_breaks)) {
  x = x_breaks[x_i]
  for (y_i in seq_along(y_breaks)) {
    y = y_breaks[y_i]
    for (z_i in seq_along(z_breaks)) {
      z = z_breaks[z_i]
      lower[1:3] = c(x, y, z)
      upper[1:3] = c(x+k*max_length, y+k*max_length, z+k*max_length)
      bounding_box = structure(list(x = c(1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0), y = c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1), z = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)), .Names = c("x", "y", "z"), row.names = c(NA, -24L), class = "data.frame")
      bounding_box = bounding_box * matrix(upper[1:3] - lower[1:3], ncol=3, nrow=24, byrow = T) + matrix(lower[1:3], ncol=3, nrow=24, byrow = T)
      segments3d(bounding_box, line_antialias = TRUE, col = "blue")
      lower_points = c(x-max_length/2, y-max_length/2, z-max_length/2)
      upper_points = c(x+(k+0.5)*max_length, y+(k+0.5)*max_length, z+(k+0.5)*max_length)
      bounding_box_points = structure(list(x = c(1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0), y = c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1), z = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)), .Names = c("x", "y", "z"), row.names = c(NA, -24L), class = "data.frame")
      bounding_box_points = bounding_box_points * matrix(upper_points[1:3] - lower_points[1:3], ncol=3, nrow=24, byrow = T) + matrix(lower_points[1:3], ncol=3, nrow=24, byrow = T)
      segments3d(bounding_box_points, line_antialias = TRUE, col = "green")
    }
  }
}
plot3d(model, type='dots', alpha=1, add=T)


all_solutions = list()
for (x_i in seq_along(x_breaks)) {
  x = x_breaks[x_i]
  all_solutions[[x_i]] = list()
  for (y_i in seq_along(y_breaks)) {
    y = y_breaks[y_i]
    all_solutions[[x_i]][[y_i]] = list()
    for (z_i in seq_along(z_breaks)) {
      z = z_breaks[z_i]
      ## 5,1,1 => the little branch
      # Initialization of the Points subset
      selected_points = which((full_Points[,1] > x-max_length/2) &
                                (full_Points[,1] < x+(k+0.5)*max_length) &
                                (full_Points[,2] > y-max_length/2) & 
                                (full_Points[,2] < y+(k+0.5)*max_length) &
                                (full_Points[,3] > z-max_length/2) &
                                (full_Points[,3] < z+(k+0.5)*max_length))
      Points = full_Points[selected_points,]
      if (nrow(Points) < 500) {
        cat('Skipping empty cube.\n')
        next # abort if empty
      }
      realized_x = range(Points[,1])
      realized_y = range(Points[,2])
      realized_z = range(Points[,3])
      cat('Doing x_i=', x_i, " ; y_i=", y_i, " ; z_i=", z_i, "\n Number of points: ",nrow(Points),"\n")
      lower[1:3] = c(max(c(x,realized_x[1])), max(c(y, realized_y[1])), max(c(z,realized_z[1])))
      upper[1:3] = c(min(c(x+k*max_length, realized_x[2])), min(c(y+k*max_length, realized_y[2])), min(c(z+k*max_length, realized_z[2])))
      Scene_cpp = new(Scene, Points, thresh, maxPopSize)
      
      # Initialization of the population
      Solutions = t(sapply(1:minPopSize, function(u) c(array(runif(length(lower), lower, upper)), rep(NA,objDim+hidDim))))
      for (i in 1:minPopSize) {
        Scene_cpp$change_indiv(Solutions[i, 1:7], i-1)
      }
      if ('All_Sols' %in% ls()) {
        j=1
        for (i in order(All_Sols[,8], decreasing = T)) {
          if (j < minPopSize) {
            if (All_Sols[i,8] > 0.1) {
              if (all(All_Sols[i, 1:7] >= lower) & all(All_Sols[i, 1:7] <= upper)) {
                Solutions[j, 1:7] = All_Sols[i, 1:7]
                Scene_cpp$change_indiv(Solutions[j, 1:7], j-1)
                j = j + 1
                cat('.')
              }
            }
          }
        }
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

# Plot the cylinders
# AND save them to displayed_sol
display_acceptance_ratio = 0.3 # Vary this threshold to change the cylinder acceptance level
cyl_nb = 0
displayed_sol = NULL
All_Sols = NULL
for (x_i in seq_along(x_breaks)) {
  for (y_i in seq_along(y_breaks)) {
    for (z_i in seq_along(z_breaks)) {
      try({
        
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
      })
    }
  }
}
plot3d(model, type='dots', alpha=0.1, add=T)


##########################
## EXTEND THE CYLINDERS ##
##########################

# Redo one last pass of optimizing with the full model

Solutions = displayed_sol # starts from the kept solutions

Points = full_Points
max_length = 1 
bbox1 = matrix(c(apply(Points, 2, min, na.rm=T),apply(Points, 2, max, na.rm=T)), ncol=3, byrow = T)
lower = c(bbox1[1,1:3], 0, 0, min_length, min_radius)
upper = c(bbox1[2,1:3], 2*pi, pi, max_length, max_radius)
thresh = 0.004220307*1.5 # larger threshold (probably useless)
acceptance_ratio = display_acceptance_ratio
extend_acceptance_ratio = display_acceptance_ratio

# Short search: 1000 iterations exactly
max_ite = 1000
constant_ite_termination = 1000
min_ite = 1000
Scene_cpp = new(Scene, Points, thresh, nrow(Solutions)*2)
for (i in 1:nrow(Solutions)) {
  Scene_cpp$change_indiv(Solutions[i, 1:7], i-1)
}

source('main_loop.R')

# Now extends the cylinders:

source('extend_cylinder.R')

# Plots the results:

for (i in 1:nrow(Solutions)) {
  if (Solutions[i, 8] > acceptance_ratio) {
    display_sol(Solutions[i, 1:7], showscene = F, col = 'purple', alpha = 0.5, sides=100)
  }
}
plot3d(Points, alpha=0.1, add=T)

# Compute the volume of the cylinders
sourceCpp('helper_fct_cpp.cpp', rebuild=F)
grid_step = 0.005
grid3d = expand.grid(x=seq(bbox1[1,1], bbox1[2,1], by=grid_step),
                     y=seq(bbox1[1,2], bbox1[2,2], by=grid_step),
                     z=seq(bbox1[1,3], bbox1[2,3], by=grid_step))
grid3d = as.matrix(grid3d)
inside = rep(0, nrow(grid3d))
for (i in 1:nrow(Solutions)) {
  if (Solutions[i, 8] > acceptance_ratio) {
    inside = inside + get_points_inside(Solutions[i, 1:7], grid3d, 0)
  }
}
plot3d(grid3d[which(inside>0),], col=inside[which(inside>0)], asp='iso')
volume = length(which(inside>0)) * grid_step^3

# Get the mesh from the intersection of cylinders
require('alphashape3d')
tst=ashape3d(grid3d[which(inside>0),], alpha = grid_step*2, pert=T, eps=1e-5)


# Save some views
plot3d(model, type='dots', alpha=0.1, xlab = '', ylab = '', zlab = '', axes = FALSE)
rgl.snapshot('tree1_original2.png', top=F)
plot(tst, col=rep('purple',3), transparency = 0.5, clear=F)
rgl.snapshot('tree1_overlayed2.png', top=F)
plot(tst, col=rep('purple',3), transparency = 1, clear=T)
rgl.snapshot('tree1_approx2.png', top=F)

