##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

create_cylinder = function(center=c(0,0,0), cyl_length=5, width=1, jitter=0.1, precision=0.1, doplot=F, orientation=c(pi/2,0))
{
  require('rgl')
  
  # precision = mean distance between two points
  lambda = seq(-cyl_length/2, cyl_length/2, by = precision)
  # lambda = seq(0, cyl_length, by = precision)
  
  
  cyl = cylinder3d(
    center = cbind(
      # center[1] + lambda, 
      # center[2], 
      # center[3]),
      center[1] - lambda*sin(orientation[2])*cos(-pi/2+orientation[1]),
      center[2] - lambda*cos(orientation[2])*cos(-pi/2+orientation[1]),
      center[3] - lambda*sin(-pi/2+orientation[1])),
    radius = width, 
    closed = T, sides=round(2*pi*width/precision))
  cyl$ib = NULL
  cyl$vb[1:3,] = cyl$vb[1:3,] + runif(-jitter, jitter, n = length(cyl$vb[1:3,]))
  if (doplot) {
    plot3d(cyl, type='dots')
  }
  return(cyl)
}
# create_cylinder(doplot = T)

create_cylinder_caps = function(center=c(0,0,0), cyl_length=5, width=1, jitter=0.1, precision=0.1, doplot=F, orientation=c(pi/2,0))
{
  require('rgl')
  
  # precision = mean distance between two points
  caps = NULL
  lambda = c(-cyl_length/2, cyl_length/2)
  for (r in seq(precision, width-precision, by=precision)) {
    cyl = cylinder3d(
      center = cbind(
        center[1] - lambda*sin(orientation[2])*cos(-pi/2+orientation[1]),
        center[2] - lambda*cos(orientation[2])*cos(-pi/2+orientation[1]),
        center[3] - lambda*sin(-pi/2+orientation[1])),
      radius = r, 
      closed = 0, sides=round(2*pi*r/precision))
    caps = rbind(caps, t(cyl$vb[1:3,]))
  }
  if (doplot) {
    plot3d(caps)
  }
  return(caps)
}
# create_cylinder_caps(doplot = T)

create_dupin_cyclide = function(b, c, d, center=c(0,0,0), doplot=T, precision=0.05)
{
  require('rgl')
  uu = seq(0, 2*pi, by = precision/2)
  vv = seq(0, 2*pi, by = precision)
  
  a = sqrt(b^2+c^2)  # also: a > b > 0 ; d >= 0
  
  dupin = NULL
  for (u in uu) {
    # for (v in vv) {
    # x = d*(c-a*cos(u)*cos(v)) + b^2*cos(u) / (a - c*cos(u)*cos(v))
    # y = b*sin(u)*(a-d*cos(v)) / (a - c*cos(u)*cos(v))
    # z = b*sin(v)*(c*cos(u)-d) / (a - c*cos(u)*cos(v))
    # dupin = rbind(dupin, c(x,y,z))
    # }
    x = d*(c-a*cos(u)*cos(vv)) + b^2*cos(u) / (a - c*cos(u)*cos(vv))
    y = b*sin(u)*(a-d*cos(vv)) / (a - c*cos(u)*cos(vv))
    z = b*sin(vv)*(c*cos(u)-d) / (a - c*cos(u)*cos(vv))
    dupin = rbind(dupin, cbind(x,y,z))
  }
  require('alphashape3d')
  
  # dupin = cbind(x, y, z)
  
  if (doplot) {
    clear3d()
    rr = range(dupin[])
    decorate3d(rr, rr, rr, box=F, axes=F, xlab='', ylab='', zlab='')
    dupin_shape = ashape3d(dupin, alpha=0.075, pert=T)
    plot(dupin_shape,  col=rep('purple',3), edges=F, vertices=F, transparency = 0.8)
  }
  return(dupin)
}
# dupin=create_dupin_cyclide(b=1,c=0.25,d=0.5, doplot = F)


show_dupin_cyclide = function(b, c, d, center=c(0,0,0), doplot=T, precision=0.05)
{
  require('rgl')
  uu = seq(0, 2*pi, by = precision)
  vv = seq(0, 2*pi, by = precision)
  
  a = sqrt(b^2+c^2)  # also: a > b > 0 ; d >= 0
  
  dupin = NULL
  open3d()
  for (u in uu) {
    x = d*(c-a*cos(u)*cos(vv)) + b^2*cos(u) / (a - c*cos(u)*cos(vv))
    y = b*sin(u)*(a-d*cos(vv)) / (a - c*cos(u)*cos(vv))
    z = b*sin(vv)*(c*cos(u)-d) / (a - c*cos(u)*cos(vv))
    lines3d(x,y,z)
  }
  for (v in vv) {
    x = d*(c-a*cos(uu)*cos(v)) + b^2*cos(uu) / (a - c*cos(uu)*cos(v))
    y = b*sin(uu)*(a-d*cos(v)) / (a - c*cos(uu)*cos(v))
    z = b*sin(v)*(c*cos(uu)-d) / (a - c*cos(uu)*cos(v))
    lines3d(x,y,z)
  }
  
}
# show_dupin_cyclide(b=1,c=0.25,d=0.5, doplot = F, precision = 0.25)



# create_torus = function(center=c(0,0,0), radius=5, width=1, jitter=0.1, precision=50, doplot=F)
# {
#   require('rgl')
#   theta <- seq(0, 2*pi, len = round(2*pi*precision))
#   cyl = cylinder3d(
#     center = cbind(
#       center[1] + radius*sin(theta), 
#       center[2] + radius*cos(theta), 
#       center[3]),
#     radius = width, 
#     closed = TRUE, sides=precision)
#   cyl$ib = NULL
#   cyl$vb[1:3,] = cyl$vb[1:3,] + runif(-jitter, jitter, n = length(cyl$vb[1:3,]))
#   if (doplot) {
#     plot3d(cyl, type='dots')
#   }
#   return(cyl)
# }
create_torus = function(center=c(0,0,0), radius=5, width=1, jitter=0.1, precision=0.1, doplot=F)
{
  require('rgl')
  # precision = mean distance between two points
  theta = seq(0, 2*pi, by = precision/radius)
  cyl = cylinder3d(
    center = cbind(
      center[1] + radius*sin(theta), 
      center[2] + radius*cos(theta), 
      center[3]),
    radius = width, 
    closed = TRUE, sides=round(2*pi*width/precision))
  cyl$ib = NULL
  cyl$vb[1:3,] = cyl$vb[1:3,] + runif(-jitter, jitter, n = length(cyl$vb[1:3,]))
  if (doplot) {
    plot3d(cyl, type='dots')
  }
  return(cyl)
}
# create_torus(doplot = T)

create_hybrid_rectangular_cuboid_cylinder = function(alpha=0.5, prop=1, center=c(0,0,0), shape_length=5, width=1, jitter=0, precision=0.1, doplot=F)
{
  shape_points = NULL
  perimeter = 2*width*4 *alpha+ 2*pi*width *(1-alpha)
  for (z in seq(0, shape_length, by=precision)) {
    for (angle in seq(0, 2*pi*prop, le=round(perimeter/precision)+1)[-1]) {
      pt = c(cos(angle)/max(abs(c(cos(angle),sin(angle)))), sin(angle)/max(abs(c(cos(angle),sin(angle)))), z) *alpha # cuboid part
      pt = pt + c(cos(angle), sin(angle), z) *(1-alpha) # cylindrical part
      shape_points = rbind(shape_points, center+c(0,0,-shape_length/2)+pt*c(width, width, 1))
    }
  }
  if (jitter > 0) {
    shape_points = shape_points + runif(-jitter, jitter, n = prod(dim(shape_points)))
  }
  if (doplot) {
    plot3d(shape_points, type='p', asp='iso')
  }
  return(shape_points)
}
# tst=create_hybrid_rectangular_cuboid_cylinder(alpha=1, doplot = T)


create_rectangular_cuboid = function(center=c(0,0,0), cube_length=5, width=1, jitter=0.1, precision=0.1, doplot=F)
{
  require('rgl')
  
  # precision = mean distance between two points
  lambda = seq(-cube_length/2, cube_length/2, by = precision)
  bigside = seq(-width/2, width/2, by = precision)
  smallside = seq(-width/2+precision, width/2-precision, by = precision)
  
  cube = NULL
  for (l in lambda) {
    for (s in bigside) {
      cube = rbind(cube,
                   c(center[1] + width/2,
                     center[2] + s,
                     center[3] + l))
      cube = rbind(cube,
                   c(center[1] - width/2,
                     center[2] + s,
                     center[3] + l))
    }
    for (s in smallside) {
      cube = rbind(cube,
                   c(center[1] + s,
                     center[2] + width/2,
                     center[3] + l))
      cube = rbind(cube,
                   c(center[1] + s,
                     center[2] - width/2,
                     center[3] + l))
    }
  }
  cube = cube + runif(-jitter, jitter, n = prod(dim(cube)))
  
  if (doplot) {
    plot3d(cube, type='p')
  }
  return(cube)
}
# unique(create_rectangular_cuboid(doplot = T))


# This functions checks the constraints and outputs warnings if they are badly defined
check_constraints = function(thresh, min_radius, min_length, max_length, Points_xor_precision=NULL)
{
  fail = F
  if (is.null(Points_xor_precision)) {
    stop('You need to supply "Points_xor_precision". Aborting.')
  }
  if (length(Points_xor_precision) > 1) {
    message('Will compute the median spacing between points in the mesh (using k=3)...')
    require('FNN')
    points_knn = get.knn(Points[complete.cases(Points[,1:3]), 1:3], k=3)
    precision = median(apply(points_knn$nn.dist, 1, median, na.rm=T), na.rm=T)
    message(paste('Found', precision))
  } else {
    precision = Points_xor_precision
  }
  if (thresh < precision * 1.5) {
    warning(paste0("the tolerance threshold ('thresh') is advised to be at least ", precision * 1.5, " -- ", precision * 1.75))
    thresh = precision * 1.5
    fail = T
  }
  if (min_radius < thresh * 4) {
    warning(paste0("the minimal radius ('min_radius') is advised to be at least ", thresh * 4, " -- ", thresh * 5))
    fail = T
  }
  if (min_length < thresh * 1.5) {
    warning(paste0("the minimal axis length ('min_length') is advised to be at least ", thresh * 1.5, " -- ", thresh * 2))
    fail = T
  }
  if (max_length > thresh * 999) {
    warning(paste0("in the current implementation, the maximal axis length ('max_length') should be at most ", thresh * 999))
    fail = T
  }
  if (fail) {
    stop('Please update the constraints. Aborting.')
  } else {
    message('Constraints are respected, moving forward.')
  }
}

# This functions returns the normalized cylinder axis
get_cylinder_axis = function(indiv)
{
  center1 = c(
    indiv[1] - indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] - indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] - sin(-pi/2+indiv[5])*indiv[6]/2
  )  
  center2 = c(
    indiv[1] + indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] + indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] + sin(-pi/2+indiv[5])*indiv[6]/2
  )
  centers = rbind(center1,center2)
  u = apply(centers, 2, diff)
  u = u / sqrt(sum(u^2))
  return(u)
}


update_points_score = function(indiv, bbox=NULL, fnn=F, first_obj='number', second_obj='', third_obj='', return_orient=F,
                               Points=Points)
{
  n_points = nrow(Points)-1
  # in this version, the phenotype is centered on the middle
  if (is.null(voxel)) {
    points_i = 1:n_points
  } else {
    points_i = which(label_voxels %in% voxels_isaboveground[which(segmented_voxels[,4]==voxel)])
  }
  
  # new phenotype
  center_middle = indiv[1:3]
  indiv[4] = max(c(indiv[4],0)) # enforce some boundaries
  indiv[4] = min(c(indiv[4],2*pi))
  indiv[5] = max(c(indiv[5],0))
  indiv[5] = min(c(indiv[5],pi))
  indiv[6] = max(indiv[6],0.05)
  
  center1 = c(
    indiv[1] - indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] - indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] - sin(-pi/2+indiv[5])*indiv[6]/2
  )  
  center2 = c(
    indiv[1] + indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] + indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] + sin(-pi/2+indiv[5])*indiv[6]/2
  )
  centers = rbind(center1,center2)
  
  r = indiv[7]
  
  candidates = rep(T, length(points_i))
  
  u = apply(centers, 2, diff)
  u = u / sqrt(sum(u^2))
  
  # # vectorized
  ps = Points[points_i[which(candidates)], 1:3]
  xds = (ps-matrix(center_middle, ncol=3, nrow=nrow(ps), byrow = T)) %*% u 
  pros = cbind(xds*u[1]+center_middle[1], xds*u[2]+center_middle[2], xds*u[3]+center_middle[3])
  
  
  # handle the out-of-main-axis points
  ## Rounded cylinder:
  #   if (center1[3] > center2[3]) {
  #     topcenter = matrix(center1,nrow=nrow(pros),ncol=3,byrow = T)
  #     bottomcenter = matrix(center2,nrow=nrow(pros),ncol=3,byrow = T)
  #   } else {
  #     topcenter = matrix(center2,nrow=nrow(pros),ncol=3,byrow = T)
  #     bottomcenter = matrix(center1,nrow=nrow(pros),ncol=3,byrow = T)
  #   }
  #   pros[pros[,3]>topcenter[,3],] = topcenter[pros[,3]>topcenter[,3],]
  #   pros[pros[,3]<bottomcenter[,3],] = bottomcenter[pros[,3]<bottomcenter[,3],]
  ## open-ended cylinder without anything protruding
  diff_dim = which(center1 != center2)[1] ## CORRECTED A BUG HERE TO AVOID THE CASE Z[1] == Z[3]
  if (center1[diff_dim] > center2[diff_dim]) {
    pros[pros[,diff_dim] > center1[diff_dim],] = NA
    pros[pros[,diff_dim] < center2[diff_dim],] = NA
  } else {
    pros[pros[,diff_dim] > center2[diff_dim],] = NA
    pros[pros[,diff_dim] < center1[diff_dim],] = NA
  }
  
  yds = sqrt(rowSums((pros-ps)^2))
  
  belongs_to_cyl = which(abs(yds-r) < thresh)
  
  if ((second_obj == 'normals') | (first_obj == c('number_normals'))) {
    if (fnn==T) {
      stop('dont use a second objective with kd-tree classification')
    }
    desired_orientation = ps[belongs_to_cyl,] - pros[belongs_to_cyl,]
    if (length(belongs_to_cyl) == 0) {
      score_normals=1e9
    } else if (length(belongs_to_cyl) == 1) {
      score_normals = acos( sum(desired_orientation*normals[points_i[belongs_to_cyl],]) / ( sqrt(sum(desired_orientation * desired_orientation) * sum(normals[points_i[belongs_to_cyl],] * normals[points_i[belongs_to_cyl],])) ) )  
      if (is.na(score_normals)) { # point exactly on the cylinder axis or undefined normal
        score_normals = 1e8
      }
    } else {
      angles = acos( rowSums(desired_orientation*normals[points_i[belongs_to_cyl],]) / ( sqrt(rowSums(desired_orientation * desired_orientation) * rowSums(normals[points_i[belongs_to_cyl],] * normals[points_i[belongs_to_cyl],])) ) )  
      score_normals = median(angles, na.rm=T)
    }
  }
  
  if (return_orient) {
    if (length(belongs_to_cyl) == 0) {
      ave_orient = c(NA, NA, NA)
      ave_dist = NA
    } else if (length(belongs_to_cyl) == 1) {
      ave_orient = ps[belongs_to_cyl,] - pros[belongs_to_cyl,]
      ave_dist = xds[belongs_to_cyl]
    } else {
      ave_orient = apply(ps[belongs_to_cyl,] - pros[belongs_to_cyl,], 2, function(x) mean(x,na.rm=T))
      ave_dist = mean(xds[belongs_to_cyl], na.rm=T)
    }
    ave_orient = ave_orient / sqrt(sum(ave_orient^2))
  }
  
  h = dist(centers)
  # surface = 2*pi*r*h + 2*pi*r^2 # surface of the closed cylinder
  surface = 2*pi*r*h # surface of the open cylinder
  # surface = h # length of the cylinder axis
  # surface = 2*pi*r^2 # surface of the disk
  # surface = 1 # no penalty
  
  
  Sys.sleep(0.01)
  cat('.')
  
  
  if (first_obj == 'number') {
    ret = - length(belongs_to_cyl)
  } else if (first_obj == 'number_surface') {
    ret = - length(belongs_to_cyl) / surface
  } else if (first_obj == 'MSE') {
    ret = mean((yds-r)^2, na.rm=T)
  } else if (first_obj == 'number_normals') {
    ret = - length(belongs_to_cyl) / score_normals
  } else if (first_obj == 'spanned_surface') {
    if (length(belongs_to_cyl) <= 1) {
      ret_vec = rep(0, n_points+1)
      if (return_orient) {
        return(list(potential=ret_vec, ave_orient=ave_orient, ave_dist=ave_dist))
      } else {
        return(ret_vec)
      }
    } else {
      cyl_orientation = ps[belongs_to_cyl,] - pros[belongs_to_cyl,]
      base_orientation = matrix(rep(cyl_orientation[1,]/sqrt(sum(cyl_orientation[1,]^2)), nrow(cyl_orientation)), ncol=3, byrow = T)
      # base_orientation = ps[belongs_to_cyl,] # WHY???
      # would be more optimized to normalize it right away
      # cyl_angles = acos( rowSums(cyl_orientation*base_orientation) / ( sqrt(rowSums(cyl_orientation * cyl_orientation)) * sqrt(rowSums(base_orientation * base_orientation)) ) )
      
      cyl_angles = acos( rowSums(cyl_orientation*base_orientation) / ( sqrt(rowSums(cyl_orientation * cyl_orientation) ) ) )
      
      # h_breaks = centered.seq(-h/2,h/2,thresh) # NOT GREAT
      # a_breaks = centered.seq(0,pi,thresh*2*pi) # WRONG: BUG
      h_breaks = full.seq(-h/2,h/2,thresh)
      a_breaks = full.seq(0, pi, pi / (2*pi*r/thresh)) # step = thresh / (2 r)
      # cat(h_breaks, '\n')
      binned_surface = matrix(NA, ncol=2, nrow=length(belongs_to_cyl))
      binned_surface[,1] = cut(xds[belongs_to_cyl], breaks = h_breaks, labels = F, include.lowest = T)
      binned_surface[,2] = cut(cyl_angles, breaks = a_breaks, labels = F, include.lowest = T)
      # ret = - nrow(unique(binned_surface))/((length(h_breaks)-1)*(length(a_breaks)-1))
      ret = binned_surface[,1] + 1000*binned_surface[,2]
    }
  }
  
  # if (length(belongs_to_cyl)) {
  #   not_cyl_ps = (1:nrow(ps))[-belongs_to_cyl]
  # } else {
  #   not_cyl_ps = (1:nrow(ps))
  # }
  # plot3d(ps[not_cyl_ps,], col = 'orange', asp=F)
  # plot3d(ps[belongs_to_cyl,], col = 'green', add=T)
  
  # if (second_obj=='normals') {
  #   ret = c(ret, score_normals)
  # } else if (second_obj == 'surface') {
  #   ret = c(ret, surface)
  # }
  # 
  # if (third_obj=='normals') {
  #   ret = c(ret, score_normals)
  # } else if (third_obj == 'surface') {
  #   ret = c(ret, surface)
  # }
  # 
  # if (return_orient) {
  #   ret = c(ret, orient, proj_length)
  # }
  
  ret_vec = rep(0, n_points+1)
  ret_vec[belongs_to_cyl] = ret
  ret_vec[n_points+1] = length(h_breaks) + 1000*length(a_breaks)
  
  if (return_orient) {
    return(list(potential=ret_vec, ave_orient=ave_orient, ave_dist=ave_dist))
  } else {
    return(ret_vec)
  }
}

display_sol = function(indiv, Points=NULL, voxel=NULL, showscene = F, col='orange', doshowvoxel=F, alpha=1, sides=25, box_xyz=NULL, colbox='blue')
{ 
  require('rgl')
  if (!is.null(Points)) {
    points = Points[1:(nrow(Points)-1), 1:3]
  } else {
    points = rbind(indiv[1:3], indiv[1:3]-3, indiv[1:3]+3)
  }
  center_middle = indiv[1:3]
  indiv[4] = max(c(indiv[4],0)) # enforce some boundaries
  indiv[4] = min(c(indiv[4],2*pi))
  indiv[5] = max(c(indiv[5],0))
  indiv[5] = min(c(indiv[5],pi))
  indiv[6] = max(indiv[6],0.05)
  
  center1 = c(
    indiv[1] - indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] - indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] - sin(-pi/2+indiv[5])*indiv[6]/2
  )  
  center2 = c(
    indiv[1] + indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] + indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] + sin(-pi/2+indiv[5])*indiv[6]/2
  )
  centers = rbind(center1,center2)
  
  r = indiv[7]
  h = dist(centers)
  surface = 2*pi*r*h + 2*pi*r^2
  
  if (showscene) {
    # rgl.clear(type='all', subscene=NA)
    plot3d(points, col = adjustcolor('black', alpha.f = 0.2) ,aspect=T, top=F)
  }
  
  tryCatch({
    cyl2 = cylinder3d(center=centers,radius=r,sides=sides)
    plot3d(cyl2, col = col, add=T, aspect=T, alpha=alpha, top=F)
  }, error = function(...){print("Solution rendering cancelled")})
  
  
  if (!is.null(box_xyz)) {
    axis_dir=get_cylinder_axis(indiv)
    radial_dir = vcrossp(rbind(axis_dir), cbind(0,0,1))
    radial_dir2 = vcrossp(rbind(axis_dir), rbind(radial_dir))
    axis_dir = axis_dir*box_xyz[1]
    radial_dir = radial_dir*box_xyz[2]
    radial_dir2 = radial_dir2*box_xyz[3]
    
    xyz0 = indiv[1:3] + (axis_dir+ radial_dir + radial_dir2)/2
    xyz1 = indiv[1:3] - (axis_dir+ radial_dir + radial_dir2)/2
    
    axis_dir = -axis_dir
    radial_dir = -radial_dir
    radial_dir2 = -radial_dir2
    
    segment_list = list(
      x=c(xyz0[1], xyz0[1]+axis_dir[1], xyz0[1], xyz0[1]+radial_dir[1], xyz0[1], xyz0[1]+radial_dir2[1],
          xyz1[1], xyz1[1]-axis_dir[1], xyz1[1], xyz1[1]-radial_dir[1], xyz1[1], xyz1[1]-radial_dir2[1],
          xyz0[1]+radial_dir[1], xyz0[1]+radial_dir[1]+radial_dir2[1],
          xyz0[1]+radial_dir2[1], xyz0[1]+radial_dir[1]+radial_dir2[1],
          xyz1[1]-radial_dir[1], xyz1[1]-radial_dir[1]-radial_dir2[1],
          xyz1[1]-radial_dir2[1], xyz1[1]-radial_dir[1]-radial_dir2[1],
          xyz1[1]-radial_dir[1], xyz1[1]-radial_dir[1]-axis_dir[1],
          xyz1[1]-radial_dir2[1], xyz1[1]-radial_dir2[1]-axis_dir[1]),
      y=c(xyz0[2], xyz0[2]+axis_dir[2], xyz0[2], xyz0[2]+radial_dir[2], xyz0[2], xyz0[2]+radial_dir2[2],
          xyz1[2], xyz1[2]-axis_dir[2], xyz1[2], xyz1[2]-radial_dir[2], xyz1[2], xyz1[2]-radial_dir2[2],
          xyz0[2]+radial_dir[2], xyz0[2]+radial_dir[2]+radial_dir2[2],
          xyz0[2]+radial_dir2[2], xyz0[2]+radial_dir[2]+radial_dir2[2],
          xyz1[2]-radial_dir[2], xyz1[2]-radial_dir[2]-radial_dir2[2],
          xyz1[2]-radial_dir2[2], xyz1[2]-radial_dir[2]-radial_dir2[2],
          xyz1[2]-radial_dir[2], xyz1[2]-radial_dir[2]-axis_dir[2],
          xyz1[2]-radial_dir2[2], xyz1[2]-radial_dir2[2]-axis_dir[2]),
      z=c(xyz0[3], xyz0[3]+axis_dir[3], xyz0[3], xyz0[3]+radial_dir[3], xyz0[3], xyz0[3]+radial_dir2[3],
          xyz1[3], xyz1[3]-axis_dir[3], xyz1[3], xyz1[3]-radial_dir[3], xyz1[3], xyz1[3]-radial_dir2[3],
          xyz0[3]+radial_dir[3], xyz0[3]+radial_dir[3]+radial_dir2[3],
          xyz0[3]+radial_dir2[3], xyz0[3]+radial_dir[3]+radial_dir2[3],
          xyz1[3]-radial_dir[3], xyz1[3]-radial_dir[3]-radial_dir2[3],
          xyz1[3]-radial_dir2[3], xyz1[3]-radial_dir[3]-radial_dir2[3],
          xyz1[3]-radial_dir[3], xyz1[3]-radial_dir[3]-axis_dir[3],
          xyz1[3]-radial_dir2[3], xyz1[3]-radial_dir2[3]-axis_dir[3]))
    
    segment_list = structure(segment_list, .Names = c("x", "y", "z"), row.names = c(NA, -length(segment_list[[1]])), class = "data.frame")
    segments3d(segment_list, line_antialias = TRUE, col = colbox, lwd=5)
    
  }
  
  aspect3d('iso')
  return(T)
}



display_sol_insetbox = function(indiv, box_xyz, Points, voxel=NULL, showscene3D = T, showscene2D = T, col='orange', colbox='blue', doshowvoxel=F, alpha=1, sides=25)
{ 
  require('rgl')
  
  center_middle = indiv[1:3]
  indiv[4] = max(c(indiv[4],0)) # enforce some boundaries
  indiv[4] = min(c(indiv[4],2*pi))
  indiv[5] = max(c(indiv[5],0))
  indiv[5] = min(c(indiv[5],pi))
  indiv[6] = max(indiv[6],0.05)
  
  center1 = c(
    indiv[1] - indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] - indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] - sin(-pi/2+indiv[5])*indiv[6]/2
  )  
  center2 = c(
    indiv[1] + indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] + indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] + sin(-pi/2+indiv[5])*indiv[6]/2
  )
  centers = rbind(center1,center2)
  
  r = indiv[7]
  h = dist(centers)
  surface = 2*pi*r*h + 2*pi*r^2
  
  # isolate points within the box
  axis_dir_norm = get_cylinder_axis(indiv)
  radial_dir_norm = vcrossp(rbind(axis_dir_norm), cbind(0,0,1))
  radial_dir_norm = radial_dir_norm / sqrt(sum(radial_dir_norm^2))
  radial_dir2_norm = vcrossp(rbind(axis_dir_norm), rbind(radial_dir_norm))

  axis_dir = axis_dir_norm*box_xyz[1]
  radial_dir = radial_dir_norm*box_xyz[2]
  radial_dir2 = radial_dir2_norm*box_xyz[3]
  
  xyz0 = indiv[1:3] + (axis_dir+ radial_dir + radial_dir2)/2
  xyz1 = indiv[1:3] - (axis_dir+ radial_dir + radial_dir2)/2
  
  axis_dir = -axis_dir
  radial_dir = -radial_dir
  radial_dir2 = -radial_dir2
  
  p1 = xyz0 + axis_dir
  p2 = xyz0 + radial_dir
  p3 = xyz0 + radial_dir2
  
  within = rep(T,nrow(Points))
  
  # Check that the dot product u.x is between u.P1 and u.P2  
  dp1 = rowSums(Points * matrix(axis_dir,ncol=3,nrow=nrow(Points),byrow=T))
  uO = sum(axis_dir * xyz0)
  uP1 = sum(axis_dir * p1)
  within[which((dp1 > uP1) | (dp1 < uO))] = F
  
  # Check that v.x is between v.P1 and v.P4
  dp2 = rowSums(Points * matrix(radial_dir,ncol=3,nrow=nrow(Points),byrow=T))
  vO = sum(radial_dir * xyz0)
  vP2 = sum(radial_dir * p2)
  within[which((dp2 < vO) | (dp2 > vP2))] = F

  # Check that w.x is between w.P1 and w.P5
  dp3 = rowSums(Points * matrix(radial_dir2,ncol=3,nrow=nrow(Points),byrow=T))
  wO = sum(radial_dir2 * xyz0)
  wP3 = sum(radial_dir2 * p3)
  within[which((dp3 < wO) | (dp3 > wP3))] = F
  
  if (showscene3D) {
    # show the points
    plot3d(Points[within,], col = adjustcolor('black', alpha.f = 0.2) ,aspect=T, top=F)
    # show the box and the solution
    display_sol(indiv, Points=NULL, showscene = F, col=col, alpha=alpha, sides=sides, box_xyz=box_xyz, colbox=colbox)
    aspect3d('iso')
  }

  if (showscene2D) {
    # compute the coordinates in 2D
    plan_Points = cbind(rowSums(Points[within,]*matrix(radial_dir_norm, ncol=3,nrow=sum(within),byrow=T)),
                        rowSums(Points[within,]*matrix(radial_dir2_norm,ncol=3,nrow=sum(within),byrow=T)))
    circle_coords = list(x = sum((center1+center2)/2*radial_dir_norm), y=sum((center1+center2)/2*radial_dir2_norm))
    bl_coords = list(x = sum(xyz1*radial_dir_norm), y=sum(xyz1*radial_dir2_norm))
    tr_coords = list(x = sum(xyz0*radial_dir_norm), y=sum(xyz0*radial_dir2_norm))
    r = indiv[7]
    # Initialize the plot
    plot(1,1,type='n', xlim=c(bl_coords$x, tr_coords$x), ylim=c(bl_coords$y,tr_coords$y),asp=1,xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
    # show the box
    rect(bl_coords$x, bl_coords$y, tr_coords$x, tr_coords$y, xpd=T, border = colbox, lwd=4)
    # show the points
    points(plan_Points, pch=20, col=grey(0.5))
    # show the solution
    require('plotrix')
    draw.circle(circle_coords$x, circle_coords$y, radius = r, border = col,lwd=4)
  }
  
  return(T)
}


box_display_sol = function(bbox, indiv, Points=NULL, voxel=NULL, showscene = F, col='orange', doshowvoxel=F, alpha=1, sides=25)
{ 
  require('rgl')
  if (!is.null(Points)) {
    points = Points[1:(nrow(Points)-1), 1:3]
  } else {
    points = rbind(indiv[1:3], indiv[1:3]-3, indiv[1:3]+3)
  }
  center_middle = indiv[1:3]
  indiv[4] = max(c(indiv[4],0)) # enforce some boundaries
  indiv[4] = min(c(indiv[4],2*pi))
  indiv[5] = max(c(indiv[5],0))
  indiv[5] = min(c(indiv[5],pi))
  indiv[6] = max(indiv[6],0.05)
  
  # reduce the length of the cylinder so that it fits within the bounding box:
  indiv[6] = sqrt(sum((bbox[1,]-bbox[2,])^2))/2
  uu=get_cylinder_axis(indiv)
  bbox_center = (bbox[1,]+bbox[2,])/2
  center_vec = bbox_center - center_middle
  center_vec = center_vec / norm (cbind(center_vec), type='2')
  indiv[1:3] = indiv[1:3] + uu*sum(uu*center_vec)
  
  center1 = c(
    indiv[1] - indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] - indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] - sin(-pi/2+indiv[5])*indiv[6]/2
  )  
  center2 = c(
    indiv[1] + indiv[6]/2*sin(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[2] + indiv[6]/2*cos(indiv[4])*cos(-pi/2+indiv[5]),
    indiv[3] + sin(-pi/2+indiv[5])*indiv[6]/2
  )
  centers = rbind(center1,center2)
  
  r = indiv[7]
  h = dist(centers)
  surface = 2*pi*r*h + 2*pi*r^2
  
  tryCatch({
    cyl2 = cylinder3d(center=centers,radius=r,sides=sides)
    plot3d(cyl2, col = col, add=F, aspect=F, alpha=alpha, top=F, xlim=c(bbox[,1]),ylim=c(bbox[,2]),zlim=c(bbox[,3]),forceClipregion=T)
    rgl.light(x=indiv[1],y=indiv[2],z=indiv[3])
    for (ii in 1:7)
      rgl.light(x=indiv[1]+rnorm(1,sd=indiv[7]),y=indiv[2]+rnorm(1,sd=indiv[7]),z=indiv[3]+rnorm(1,sd=indiv[7]))
  }, error = function(...){print("Solution rendering cancelled")})
  
  if (showscene) {
    # rgl.clear(type='all', subscene=NA)
    selected_points = which((Points[,1] > bbox[1,1]) &
                              (Points[,1] < bbox[2,1]) &
                              (Points[,2] > bbox[1,2]) & 
                              (Points[,2] < bbox[2,2]) &
                              (Points[,3] > bbox[1,3]) &
                              (Points[,3] < bbox[2,2]))
    plot3d(points[selected_points,], col = adjustcolor('black', alpha.f = 0.2) , top=F, add=T)
  }
  
  ref = uu
  r=sqrt(sum(ref[1:3]^2))
  theta = acos(ref[3]/r)
  phi = atan2(ref[2],ref[1])
  # view3d(userMatrix=rotationMatrix(theta/pi*360,(phi+pi)/2/pi*180,-1,-1))
  
  # aspect3d('iso')
  return(cyl2)
}




combine_cyl_mesh = function(Solutions, precision)
{
  require('alphashape3d')
  all_cyls = NULL
  for (i in 1:nrow(Solutions)) {
    all_cyls = rbind(all_cyls, t(create_cylinder(center=Solutions[i,1:3], cyl_length=Solutions[i,6], width=Solutions[i,7], jitter=0, precision=0.01, doplot=F, orientation=Solutions[i,c(5,4)])$vb[1:3,]))
    all_cyls = rbind(all_cyls, create_cylinder_caps(center=Solutions[i,1:3], cyl_length=Solutions[i,6], width=Solutions[i,7], jitter=0, precision=0.01, doplot=F, orientation=Solutions[i,c(5,4)]))
  }
  sh = ashape3d(all_cyls, alpha = precision/2, pert=T)
  return(sh)
}

generate_random_solution = function(lower, upper, objDim, hidDim, heuristic=NULL, Points=NULL)
{
  sol = NULL
  if (is.null(heuristic)) {
    sol = runif(length(lower), lower, upper)
  } else if (heuristic == 'find_segment') {
    if (is.null(Points)) {
      stop(paste0('Needs to have point for heuristic ',heuristic))
    } else {
      require('FNN')
      
      sol = rep(NA, 7)
      sol[7] = runif(1, lower[7], upper[7]) # draw a random radius
      sol[6] = runif(1, lower[6], upper[6]) # random length
      
      n = 1000
      points_list = get.knnx(Points[(1:(nrow(Points)-1)),1:3], rbind(Points[sample.int(nrow(Points)-1,1),1:3]), k=n)$nn.index
      leg1 = Points[points_list[1],1:3]
      leg2 = Points[points_list[n],1:3]
      
      # leg1 = runif(3, lower, upper)
      # leg2 = runif(3, lower, upper)
      # leg2 = leg1+c(1,1,1)
      
      
      leg_center = (leg1+leg2)/2
      
      a = leg2 - leg1
      a = a/sqrt(sum(a^2))
      a = a*sol[6]
      b = runif(3); b = b / sqrt(sum(b^2))
      
      perp = c(a[2] * b[3] - a[3] * b[2],
               a[3] * b[1] - a[1] * b[3],
               a[1] * b[2] - a[2] * b[1]
      )
      perp = perp/sqrt(sum(perp^2))
      
      sol[4] = atan2(a[1], a[2])
      if (sol[4] < 0) {
        sol[4] = 2*pi + sol[4]
      }
      sol[5] = pi - acos(a[3] / sol[6])
      sol[1:3] = leg_center + perp * sol[7]
      # cat('sol4 = ', sol[4], '      sol5 = ',sol[5], '\n')
      plot3d(sol[1], sol[2], sol[3], xlab='X')
      spheres3d(rbind(leg1), add=T, col = 'red', radius=0.05)
      spheres3d(rbind(leg2), add=T, col = 'purple', radius=0.05)
      spheres3d(rbind(leg_center), add=T, col = 'yellow', radius=0.025)
      spheres3d(rbind(sol[1:3]), add=T, col = 'green', radius=0.025)
      display_sol(sol, Points=NULL, showscene = F, alpha=0.5)
    }
  } else if (heuristic == 'CGAL_Surface_mesh_skeletonization') {
    # todo
  } else {
    stop('No heuristic given')
  }
  # for (i in 1:7) {
  #   if (sol[i] < lower[i]) {
  #     sol[i] = lower[i]
  #   } else if (sol[i] > upper[i]) {
  #     sol[i] = upper[i]
  #   }
  # }
  return (c(sol, rep(NA,objDim+hidDim)))
}

# sh = combine_cyl_mesh(Solutions[1:kept_solutions,1:7], thresh)
# all_scene = color.read.ply('/media/jean/ext4/wooddebris/RestorationLog_forRplot.ply', ShowSpecimen = F)
# all_scene = color.read.ply('/media/jean/ext4/wooddebris/RestorationLog_noground_forRanalysis.ply', ShowSpecimen = F)
# bp = color.read.ply('/media/jean/ext4/wooddebris/RestorationLog_noground_forRanalysis_BP.ply', ShowSpecimen = F)
# bp = color.read.ply('/media/jean/ext4/wooddebris/RestorationLog_noground_forRanalysis_Poisson.ply', ShowSpecimen = F)
# plot3d(all_scene, type='dots', xlab='', ylab='', zlab='', box=F, axes=F)
# plot.ashape3d(sh, col=rep('green',3), clear=F)
# 
# sh = combine_cyl_mesh(Solutions[1:kept_solutions,1:7], thresh*10)
# plot.ashape3d(sh, col=rep('green',3), clear=T)
# 
# # 
# # 
# saved_par3d = dput(par3d())
# saved_par3d2 = dput(par3d())
# # plot3d(all_scene$vb[1,], all_scene$vb[2,], all_scene$vb[3,], type='n',  xlab='', ylab='', zlab='', box=F, axes=F)
# par3d(saved_par3d)
# par3d(saved_par3d2)
# # plot.ashape3d(sh, col=rep('green',3), clear=F)
# # 
# # conv_factor = 0.3048 / 0.44
# 
# sh = combine_cyl_mesh(Solutions[1:kept_solutions,1:7], thresh*20)
# plot3d(tree, type='dots', xlab='', ylab='', zlab='', box=F, axes=F, alpha=0.1)
# plot.ashape3d(sh, col=rep('green',3), triangles=T, edges=F, vertices=F, clear = F)
# writeWebGL('/media/jean/ext4/wooddebris/FrontYardTreeLeafOff3_MODEL_webGL')
# 
# snapshot3d('~/Desktop/Poisson+cloud.png')
# 
# snapshot3d('~/Desktop/input_ground.png')
# 
# 
# plot3d(bp, type='wire', color='green', xlab='', ylab='', zlab='', box=F, axes=F, add=T)
# plot3d(bp, type='shade', color='green', alpha=0.4, xlab='', ylab='', zlab='', box=F, axes=F, add=T)
