##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

color.read.ply = function (file, ShowSpecimen = TRUE, addNormals = TRUE) 
{
  require('geomorph')
  plyfile <- scan(file = file, what = "char", sep = "\n", strip.white = TRUE, 
                  quiet = TRUE)
  is.ply <- grep("ply", plyfile)
  if ((length(is.ply) == 0)) 
    stop("File is not a PLY file")
  format <- unlist(strsplit(grep(c("format "), plyfile, value = TRUE), 
                            " "))
  if (format[2] != "ascii") 
    stop("PLY file is not ASCII format: ", "format = ", format[2:length(format)])
  poly <- NULL
  material <- NULL
  xline <- unlist(strsplit(grep(c("vertex "), plyfile, value = TRUE), 
                           " "))
  npoints <- as.numeric(xline[grep(c("vertex"), xline) + 1])
  yline <- unlist(strsplit(grep(c("element face"), plyfile, 
                                value = TRUE), " "))
  npoly <- as.numeric(yline[grep(c("face"), yline) + 1])
  headerend <- grep(c("end_header"), plyfile)
  ncolpts <- (length(grep(c("property"), plyfile)) - 1)
  cols <- grep(c("property"), plyfile, value = TRUE)
  x <- grep(c(" x"), cols)
  y <- grep(c(" y"), cols)
  z <- grep(c(" z"), cols)
  points <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            1):(headerend + npoints)], " "))))
  dim(points) <- c(ncolpts, npoints)
  xpts <- points[x, ]
  ypts <- points[y, ]
  zpts <- points[z, ]
  vertices <- rbind(xpts, ypts, zpts, 1)
  if (yline[3] == 0) 
    print("Object has zero faces")
  if (yline[3] != 0) {
    poly <- as.matrix(as.numeric(unlist(strsplit(plyfile[(headerend + 
                                                            npoints + 1):(headerend + npoints + npoly)], " "))))
    dim(poly) <- c((poly[1] + 1), npoly)
    poly <- poly[-1, ]
    poly = poly + 1
  }
  colinfo <- grep("property uchar red", plyfile)
  if (length(colinfo) != 0) {
    color <- rgb(points[4, ], points[5, ], points[6, ], maxColorValue = 255)
    material$color <- color
  }
  mesh <- list(vb = vertices, it = poly, primitivetype = "triangle", material=material)
  class(mesh) <- c("mesh3d", "shape3d")
  if (addNormals == TRUE) {
    mesh <- addNormals(mesh)
  }
  if (ShowSpecimen == TRUE) {
    clear3d()
    if (length(poly) == 0) {
      dot3d(mesh)
    }
    if (length(material) != 0) {
      shade3d(mesh)
    }
    shade3d(mesh, color = "gray")
  }
  return(mesh)
}


get_spacing = function(Points)
{
  require('FNN')
  n2 = get.knn(unique(Points), k=2)
  return(quantile(n2$nn.dist[,2], probs = c(0.05, 0.1, 0.25, 0.5), na.rm=T))
}


re_orient_normals = function(scene_full, plot=T, cutoff=NULL)
{
  require('Rvcg')
  # align the scene according to the normal median
  if (plot) {
    scene_full_back = scene_full
    scene_full_back$material$color = rep('grey',length(scene_full_back$material$color))
  }
  
  scene_full = vcgUpdateNormals(scene_full)
  normals = scene_full$normals[1:3,]
  if (is.null(cutoff)) {
    meanvec = apply(normals,1,median)
  } else {
    browser()
    #     dis = t(normals) - matrix(c(0,0,1),ncol=3,nrow=ncol(scene_full$normals),byrow = T)
    #     dis = sqrt(rowSums(dis*dis))
    #     normals[1:3,which(dis>1)] = - normals[1:3,which(dis>1)]
    #     dis = t(normals) - matrix(c(0,0,1),ncol=3,nrow=ncol(scene_full$normals),byrow = T)
    #     dis = sqrt(rowSums(dis*dis))
    #     meanvec = apply(normals[1:3,which(dis<cutoff)],1,mean)
  }
  
  current_orientation=rbind(meanvec)
  desired_orientation=rbind(c(0,0,1))
  
  scene_full = align_cloud(desired_orientation, current_orientation, scene_full)
  
  if (plot) {
    par(mfrow=c(1,2))
    dis = t(scene_full$normals[1:3,]) - matrix(c(0,0,1),ncol=3,nrow=ncol(scene_full$normals),byrow = T)
    hist(sqrt(rowSums(dis*dis)),
         main='before:')
    scene_full = vcgUpdateNormals(scene_full)
    dis = t(scene_full$normals[1:3,]) - matrix(c(0,0,1),ncol=3,nrow=ncol(scene_full$normals),byrow = T)
    hist(sqrt(rowSums(dis*dis)),
         main='after:')
    meanvec = apply(scene_full$normals[1:3,],1,median)
    plot3d(scene_full, type='dots')
    dot3d(scene_full_back)
    meancyl = cylinder3d(center=rbind(cbind(0,0,0),meanvec),radius=0.05,sides=25)
    meancyl2 = cylinder3d(center=rbind(cbind(0,0,0),c(0,0,1)),radius=0.05,sides=25)
    plot3d(meancyl,add=T,color='blue')
    plot3d(meancyl2,add=T,color='red')
    sqrt(sum(desired_orientation * desired_orientation))
    
  }
  return (scene_full)
}



vcrossp = function(a, b) {
  # vectorized cross-product
  result = matrix(NA, nrow(a), 3)
  result[,1] = a[,2] * b[,3] - a[,3] * b[,2]
  result[,2] = a[,3] * b[,1] - a[,1] * b[,3]
  result[,3] = a[,1] * b[,2] - a[,2] * b[,1]
  return (result)
}


plot_scene = function(PF,PR1,PR2=NULL)
{
  # debug function to make sure that rotation etc are bug-free
  plot3d(NULL,
         xlim=c(-1,2),
         ylim=c(-1,2),
         zlim=c(-1,2),
         xlab='x',
         ylab='y',
         zlab='z',
         add=F)
  spheres3d(0,0,0,radius=0.025,color='orange')
  text3d(0.5,0,0,'x')
  text3d(0,0.5,0,'y')
  text3d(0,0,0.5,'z')
  PFcyl = cylinder3d(PF,radius=0.05,sides=25)
  shade3d(PFcyl,color='#FF00FF',add=T)
  spheres3d(PF[2,1],
            PF[2,2],
            PF[2,3],
            radius=0.075,color='#FF00FF')
  PR1cyl = cylinder3d(PR1,radius=0.05,sides=25)
  shade3d(PR1cyl,color='pink',add=T)
  spheres3d(PR1[2,1],
            PR1[2,2],
            PR1[2,3],
            radius=0.075,color='pink')
  if (!is.null(PR2)) {
    PR2cyl = cylinder3d(PR2,radius=0.05,sides=25)
    shade3d(PR2cyl,color='red',add=T)
    spheres3d(PR2[2,1],
              PR2[2,2],
              PR2[2,3],
              radius=0.075,color='red')
  }
}

align_cloud = function(desired_orientation, current_orientation, mesh)
{
  # ensure the right dimensionality of the input
  if (is.null(dim(desired_orientation))) {
    desired_orientation = rbind(desired_orientation)
  }
  if (is.null(dim(current_orientation))) {
    current_orientation = rbind(current_orientation)
  }
  
  #   angle = 2 * atan(norm(desired_orientation*norm(current_orientation) - norm(desired_orientation)*current_orientation)
  #                    / norm(desired_orientation * norm(current_orientation) + norm(desired_orientation) * current_orientation))
  angle = acos( sum(desired_orientation*current_orientation) / ( sqrt(sum(desired_orientation * desired_orientation)) * sqrt(sum(current_orientation * current_orientation)) ) )  
  
  ref_vec = vcrossp(desired_orientation,current_orientation)
  return(rotate3d(mesh, angle, ref_vec[1], ref_vec[2], ref_vec[3]))
}

# PF = matrix(c(0,0,0,runif(3)), ncol=3,byrow = T) # fixed point, fushia  	
# PR = matrix(c(1,1,1,runif(3)), ncol=3,byrow = T) # rose, to be rotated
# PR2 = align_cloud(PF[2,], PR[2,]-PR[1,], PR)
# plot_scene(PF,PR,PR2)
# 
# 
# for (i in 1:9) {
#   PF = matrix(runif(6), ncol=3) # fixed point, fushia  	
#   PR = matrix(runif(6), ncol=3) # rose, to be rotated
#   
#   x = rbind(PF[2,] - PF[1,])
#   y = rbind(PR[2,] - PR[1,])
#   
#   angle = 2 * atan(norm(x*norm(y) - norm(x)*y)
#                    / norm(x * norm(y) + norm(x) * y))
#   
# 
#   ref_vec = vcrossp(x,y)
# #   PR2 = translate3d(rotate3d(
# #     translate3d(PR, -PR[1,1], -PR[1,2], -PR[1,3]), angle, ref_vec[1], ref_vec[2], ref_vec[3]),
# #     PR[1,1], PR[1,2], PR[1,3])
#   
#   PR2 = translate3d(rotate3d(
#          translate3d(PR, -PF[1,1], -PF[1,2], -PF[1,3]), angle, ref_vec[1], ref_vec[2], ref_vec[3]),
#          PF[1,1], PF[1,2], PF[1,3])
#   
#   open3d()
#   plot_scene(PF,PR,PR2); text3d(1,1,1,paste(round(angle/pi,digits = 2),'pi'))
#   
# }


score_rotation = function(a,ref_vec='x', angle=0)
{
  if (ref_vec=='x') {
    ref_vec = c(1,0,0)
  } else if (ref_vec=='y') {
    ref_vec = c(0,1,0)
  } else {
    cat('not a good value for ref_vec, should be "x" or "y"\n')
    return(-1)
  }
  cuttedv = cut((rotate3d(a, 
                          angle,
                          ref_vec[1], 
                          ref_vec[2], 
                          ref_vec[3]))$vb[3,],
                include.lowest = T, 
                labels = F,
                breaks=seq(min(a$vb[3,],na.rm = T), 
                           max(a$vb[3,],na.rm = T), 
                           le=101))
  return(max(table(cuttedv)))
}

rotate_best_angle = function(a, maxangle, le=21)
{
  best_score = 0
  best_ax = NA
  best_angle = NA
  pb = txtProgressBar(min = 0, max = le*2, style=3); pbi = 0
  for (ax in c('x','y')) {
    for (angle in seq(-maxangle,maxangle,le=le)) {
      score = score_rotation(a, ref_vec=ax, angle)
      if (score > best_score) {
        best_score = score
        best_ax = ax
        best_angle = angle
      }
      pbi = pbi + 1
      setTxtProgressBar(pb,pbi)
    }
  }
  cat('\nThe best angle found was:',best_angle,'around axis',best_ax,'\n')
  return(list(angle=best_angle,ax=best_ax))
}

optimize_horizontal_slow = function(a, factors = c(1,10,100))
{
  # might be stuck
  for (factor in factors) {
    cont = T
    cat('\nStarting with precision',factor,'...\n')
    while (cont) {
      best = rotate_best_angle(a, pi/4/factor)
      if (best$angle != 0) {
        a = rotate3d(a, best$angle, best$ax=='x', best$ax=='y', 0)
      } else {
        cont = F
      }
    }
  }  
  return (a)
}



optimize_horizontal = function(a, maxangle = pi/4, factors = c(1,10,100,1000), le=10,verbose=F)
{
  b = a
  a$normals = NULL
  a$material = NULL
  a$vb[1,] = a$vb[1,] - median(a$vb[1,])
  a$vb[2,] = a$vb[2,] - median(a$vb[2,])
  a$vb[3,] = a$vb[3,] - median(a$vb[3,])
  pb = txtProgressBar(min = 0, max = 2*le*length(factors)*4, style=3); pbi = 0
  for (factor in factors) {
    for (ax in c('x','y','x','y')) {
      for (direction in c(1, -1)) {
        lagg = -1 # we allow score regression only 2 times 
        best_score = max(table(cut(a$vb[3,],
                                   include.lowest = T, 
                                   labels = F,
                                   breaks=seq(min(a$vb[3,],na.rm = T), 
                                              max(a$vb[3,],na.rm = T), 
                                              le=10*factor+1))))
        best_angle = 0
        if (direction == 1) {
          seq_angle = seq(0,maxangle/factor,le=le+1)[2:(le+1)]
        } else {
          seq_angle = seq(0,-maxangle/factor,le=le+1)[2:(le+1)]
        }
        for (angle_i in seq_along(seq_angle)) {
          angle = seq_angle[angle_i]
          pbi = pbi + 1
          setTxtProgressBar(pb,pbi)
          # score = score_rotation(a, ref_vec=ax, angle)
          score = max(table(cut((rotate3d(a, 
                                          angle,
                                          ax=='x', 
                                          ax=='y', 
                                          1))$vb[3,],
                                include.lowest = T, 
                                labels = F,
                                breaks=seq(min(a$vb[3,],na.rm = T), 
                                           max(a$vb[3,],na.rm = T), 
                                           le=10*factor+1))))
          if (score > best_score) {
            best_score = score
            best_angle = angle
          } else {
            if (lagg >= 0) {
              # cat('should break!\n\n')
              break
            } else {
              lagg = lagg + 1
            }
          }
        }
        if (best_angle != 0) {
          if (verbose) {
            cat('performed a rotation of',best_angle,'around axis',ax,'\n\n')
          }
          a = rotate3d(a, best_angle, ax=='x', ax=='y', 0)
          a$vb[1,] = a$vb[1,] - median(a$vb[1,])
          a$vb[2,] = a$vb[2,] - median(a$vb[2,])
          a$vb[3,] = a$vb[3,] - median(a$vb[3,])
        }
        pbi = pbi + (le - angle_i)
        setTxtProgressBar(pb,pbi)
      }
    }
  }  
  b$vb = a$vb
  return (b)
}

# expandMesh = function(m)
# {
#   # function to repeat (expand) single-coded components of a mesh
#   ## for now, handles only "quad" primitive types
#   n = ncol(m$ib)
#   for (field in names(m$material)) {
#     if (length(m$material[[field]]) == 1) {
#       m$material[[field]] = rep(m$material[[field]], n)
#     }
#   }
#   return (m)
# }
# 
# concatenateMeshes = function(...)
# {
#   args = list(...)
#   if (length(args) == 1 && !inherits(args[[1]], "mesh3d")) 
#     args = unlist(args, recursive = FALSE)
#   argc = length(args)
#   if (argc < 2) 
#     stop("at least two arguments needed")
#   outmesh = expandMesh(args[[1]])
#   for (mi in 2:argc) {
#     m = expandMesh(args[[mi]])
#     outmesh$ib = cbind(outmesh$ib, m$ib+ncol(outmesh$vb))
#     outmesh$vb = cbind(outmesh$vb, m$vb)
#     for (field in names(m$material)) {
#       outmesh$material[[field]] = c(outmesh$material[[field]], m$material[[field]])
#     }
#   }
#   return (outmesh)
# }

