### These functions were taken from package 'nsga2R' and then modified to include the 
### custom mutation/cross-over operators

mutate = function(y, yl, yu, mum)
{
  if (y > yl) {
    ## modification to avoid boundary-stickiness
    # if ((y - yl) < (yu - y)) {
    #   delta = (y - yl)/(yu - yl)
    # }
    # else {
    #   delta = (yu - y)/(yu - yl)
    # }
    rnd = runif(1)
    mut_pow = 1/(mum + 1)
    if (rnd < 0.5) {
      delta = (y - yl)/(yu - yl) # modif
      xy = 1 - delta
      val = 2 * rnd + (1 - 2 * rnd) * (xy^(mum + 1))
      deltaq = val^mut_pow - 1
    }
    else {
      delta = (yu - y)/(yu - yl) # modif
      xy = 1 - delta
      val = 2 * (1 - rnd) + 2 * (rnd - 0.5) * (xy^(mum + 1))
      deltaq = 1 - val^mut_pow
    }
    y = y + deltaq * (yu - yl)
    if (y > yu) {
      y = yu
    }
    else if (y < yl) {
      y = yl
    }
    return (y)
  } else {
    xy = runif(1)
    return (yl + xy * (yu - yl))
  }
}

boundedPolyMutation = function (parent_chromosome, varNo, objDim, hidDim,
                                lowerBounds, upperBounds, mprob, mum) 
{
  popSize = nrow(parent_chromosome)
  child <- parent_chromosome
  for (i in 1:popSize) {
    for (j in 1:varNo) {
      if (runif(1) < mprob) {
        child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
      }
    }
  }
  return(child)
}

# New May 16
forced_mutation3D = function (parent_chromosome, varNo, objDim, hidDim, lowerBounds, upperBounds, mprob, mum, allowed_operators, verbose_return=F)
{
  parent_chromosome = rbind(parent_chromosome)
  child = parent_chromosome
  stuck = 0
  while (identical(child[1:varNo], parent_chromosome[1:varNo])) {
    stuck = stuck + 1
    if (stuck > 100) {break} # to solve cases where it is not possible to generate a new indiv given the allowed_operators
    if (verbose_return) {
      child_n_operator = mutation3D(parent_chromosome, varNo, objDim, hidDim,lowerBounds, upperBounds, mprob, mum, allowed_operators, verbose_return=T)
      child = child_n_operator$child
    } else {
      child = mutation3D(parent_chromosome, varNo, objDim, hidDim,lowerBounds, upperBounds, mprob, mum, allowed_operators)
    }
  }
  if (verbose_return) {
    return(child_n_operator)
  } else {
    return(child)
  }
}

# New May 16
fill_russian_roulette = function (n_to_fill, fitness, parent_chromosome, varNo, objDim, hidDim, lowerBounds, upperBounds, mprob, cprob, mum, muc, m_c_ratio, allowed_operators = 1:8, verbose=F)
{
  if (is.null(ncol(parent_chromosome))) {
    # forces parent_population to be a matrix even if it contains only one solution
    parent_population = which(complete.cases(rbind(parent_chromosome)))
  } else {
    parent_population = which(complete.cases(parent_chromosome[,1:(varNo+objDim+hidDim)]))
  }
  if (length(parent_population) == 0) {
    # we arrive here if no complete solution exists in the parent population -> bug or wrong initialization
    warning('Bad model initialization\n')
    parent_chromosome = NULL
  } else if (length(parent_population) == 1) {
    parent_chromosome = matrix(parent_chromosome[parent_population,], nrow=1)
  } else {
    parent_chromosome = parent_chromosome[parent_population, ]
  }
  if (verbose) {
    cat("Parent population of size ", dim(parent_chromosome), "\n")
    for (i in 1:nrow(parent_chromosome)) {
      for (j in 1:ncol(parent_chromosome)) {
        cat(parent_chromosome[i,j], " ")
      }
      cat("\n")
    }
  }
  fitness = fitness[parent_population] # handles NA in the parent chromosomes
  fitness = fitness / sum(fitness, na.rm=T) ## BEWARE!!! THIS MAXIMIZES THE FITNESS
  if (any(is.na(fitness))) {
    fitness[] = 1/length(fitness)
  }
  childs = matrix(NA, ncol=(varNo+objDim+hidDim), nrow=n_to_fill)
  all_selected = NULL
  for (i in 1:n_to_fill) {
    if (all(is.na(fitness))) {
      # falls back on random generation if all fitness are NA (or length(fitness) == 0) 
      cat('Generating random offspring because population contains no valid solution...\n')
      new_child = generate_random_solution(lowerBounds, upperBounds, objDim, hidDim, heuristic=NULL, Points=NULL)
    } else if (runif(1) <= m_c_ratio) {
      # mutation
      selected = which(cumsum(fitness) >= runif(1))[1]
      all_selected = c(all_selected, parent_population[selected])
      new_child = forced_mutation3D(parent_chromosome[selected,], varNo, objDim, hidDim, lowerBounds, upperBounds, mprob, mum, allowed_operators = allowed_operators)
      if (verbose) {
        cat("Mutated current solution:\n")
        for (j in 1:ncol(parent_chromosome)) {
          cat(parent_chromosome[selected,j], " ")
        }
        cat("\nInto:\n")
        for (j in 1:length(new_child)) {
          cat(new_child[j], " ")
        }
        cat("\n")
      }
    } else {
      # cross-over
      selected1 = which(cumsum(fitness) >= runif(1))[1]
      for (redo in 1:100) {
        # tries hard to select another parent
        selected2 = which(cumsum(fitness) >= runif(1))[1]
        if (!identical(selected1, selected2))
          break
      }
      if (!identical(selected1, selected2)) {
        # 2 different parents => do the cross-over
        new_child = Xover3D(parent_chromosome[c(selected1, selected2),], varNo, objDim, hidDim, lowerBounds, upperBounds, cprob, muc)[sample(c(1,2),size=1),]
        while (identical(new_child, parent_chromosome[selected1,]) | identical(new_child, parent_chromosome[selected2,])) {
          new_child = Xover3D(parent_chromosome[c(selected1, selected2),], varNo, objDim, hidDim, lowerBounds, upperBounds, cprob, muc)[sample(c(1,2),size=1),]
        }
        all_selected = c(all_selected, parent_population[selected1], parent_population[selected2])
        if (verbose) {
          cat("X-over of current solutions:\n")
          for (j in 1:ncol(parent_chromosome)) {
            cat(parent_chromosome[selected1,j], " ")
          }
          cat("\n")
          for (j in 1:ncol(parent_chromosome)) {
            cat(parent_chromosome[selected2,j], " ")
          }
          cat("\nInto:\n")
          for (j in 1:length(new_child)) {
            cat(new_child[j], " ")
          }
          cat("\n")
        }
      } else {
        # could not find 2 different parents => resorts to random generation
        cat('Generating random offspring because two solutions were not found for cross-over...\n')
        new_child = generate_random_solution(lowerBounds, upperBounds, objDim, hidDim, heuristic=NULL, Points=NULL)
        if (verbose) {
          cat("X-over failed! Generating random solution:\n")
          for (j in 1:length(new_child)) {
            cat(new_child[j], " ")
          }
          cat("\n")
        }
      }
    }
    childs[i, 1:varNo] = new_child[1:varNo]
  }
  # childs[, (varNo+1):(varNo+objDim+hidDim)] = NA
  return(list(childs=childs, parents=unique(all_selected)))
}


mutation3D = function (parent_chromosome, varNo, objDim, hidDim,
                       lowerBounds, upperBounds, mprob, mum, allowed_operators = 1:8, verbose_return=F)
{
  debug = NULL
  debug_plotnew = F
  debug_time = 5
  popSize = nrow(parent_chromosome)
  child = parent_chromosome
  selected_operators = rep(F, 7)
  for (i in 1:popSize) {
    ## TODO: symmetric rotation along the orient vector
    #     child = indiv #garbage
    #     u = apply(centers,2,diff)
    #     u = u / sqrt(sum(u^2))
    #     center_middle = apply(centers,2,mean)
    #     rot_axis1 = center_middle + u * proj_length
    #     rot_axis2 = center_middle + u * proj_length + orient * r
    #     rot_axis = rot_axis2 - rot_axis1
    #     
    #     xy_over_z = sqrt(sum((r*rot_axis[1:2])^2)) / abs(r*rot_axis[3])
    #     angle = pi/2
    #     
    #     rotindiv = c(rot_axis1+rotate3d(-indiv[1:3]+rot_axis1, angle, rot_axis[1], rot_axis[2], rot_axis[3]),
    #                  (indiv[4]+angle*(1-xy_over_z)) %% (2*pi),
    #                  (indiv[5]+angle*xy_over_z) %% pi,
    #                  indiv[6:7])
    #     display_sol(indiv,showscene=T,voxel = voxel, col='blue')
    #     display_sol(rotindiv,showscene=F,voxel = voxel, col='black')
    # move the center toward the best contact area (only if using hidDim)
    if (hidDim > 0) {
      if (1 %in% allowed_operators) {
        if ((runif(1) < mprob) & (!is.na(child[i,7+objDim+4]))) {
          selected_operators[1]=T
          if (1 %in% debug) {
            cat('\nRecenter Translation\n')
            display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
            spheres3d(child[i, 1:3],radius=0.01, col='red')
          }
          u = c(
            child[i,6]/2*sin(child[i,4])*cos(-pi/2+child[i,5]),
            child[i,6]/2*cos(child[i,4])*cos(-pi/2+child[i,5]),
            sin(-pi/2+child[i,5])*child[i,6]/2
          )
          u = u / sqrt(sum(u^2))
          for (j in 1:3) {
            child[i, j] = child[i, j] + u[j]*child[i,7+objDim+4]
            # ensure that we still are within the right boundaries
            if (child[i,j] < lowerBounds[j]) {
              child[i,j] = lowerBounds[j]
            } else if (child[i,j] > upperBounds[j]) {
              child[i,j] = upperBounds[j]
            }
          }
          # for (j in 4:5) {
          #   child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
          # }
          if (1 %in% debug) {
            display_sol(child[i,1:7], voxel, showscene = F, col='green')
            spheres3d(child[i, 1:3],radius=0.01, col='green')
            Sys.sleep(debug_time)
          }
        }
      }
    }
    # symmetric flip along orient (only if using hidDim)
    if (hidDim > 0) {
      if (2 %in% allowed_operators) {
        if ((runif(1) < mprob) & (!is.na(child[i,7+objDim+1]))) {
          selected_operators[2]=T
          if (2 %in% debug) {
            cat('\nSymmetric Flip\n')
            display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
          }
          r = child[i, 7]
          for (j in 1:3) {
            child[i, j] = child[i, j] + 2*r*child[i, 7+objDim+j]
            # ensure that we still are within the right boundaries
            if (child[i,j] < lowerBounds[j]) {
              child[i,j] = lowerBounds[j]
            } else if (child[i,j] > upperBounds[j]) {
              child[i,j] = upperBounds[j]
            }
          }   
          if (2 %in% debug) {
            display_sol(child[i,1:7], voxel, showscene = F, col='green')
            Sys.sleep(debug_time)
          }
        }
      }
    }
    # orient-based dilation (only if using hidDim)
    if (hidDim > 0) {
      if (3 %in% allowed_operators) {
        if ((runif(1) < mprob) & (!is.na(child[i,7+objDim+1]))) {
          selected_operators[3]=T
          if (3 %in% debug) {
            cat('\nOrient-based dilation\n')
            display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
          }
          r = child[i, 7]
          r_new = mutate(r, lowerBounds[7], upperBounds[7], 1) # 1 instead of mum, to really explore far away
          child[i, 7] = r_new
          for (j in 1:3) {
            child[i, j] = child[i, j] + (r-r_new)*child[i, 7+objDim+j]
            # ensure that we still are within the right boundaries
            if (child[i,j] < lowerBounds[j]) {
              child[i,j] = lowerBounds[j]
            } else if (child[i,j] > upperBounds[j]) {
              child[i,j] = upperBounds[j]
            }
          }
          if (3 %in% debug) {
            display_sol(child[i,1:7], voxel, showscene = F, col='green')
            Sys.sleep(debug_time)
          }
        }
      }
    }
    # random translation
    if (4 %in% allowed_operators) {
      if (runif(1) < mprob) {
        selected_operators[4]=T
        if (4 %in% debug) {
          cat('\nRandom translation\n')
          display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
        }
        for (j in 1:3) {
          child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
        }
        if (4 %in% debug) {
          display_sol(child[i,1:7], voxel, showscene = F, col='green')
          Sys.sleep(debug_time)
        }
      }
    }
    # random rotation
    if (5 %in% allowed_operators) {
      if (runif(1) < mprob) {
        selected_operators[5]=T
        if (5 %in% debug) {
          cat('\nRandom rotation\n')
          display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
        }
        for (j in 4:5) {
          child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
        }
        if (5 %in% debug) {
          display_sol(child[i,1:7], voxel, showscene = F, col='green')
          Sys.sleep(debug_time)
        }
      }
    }
    # random dilation/elongation
    for (j in 6:7) {
      if (j %in% allowed_operators) {
        if (runif(1) < mprob) {
          selected_operators[j]=T
          if (j %in% debug) {
            cat('\nRandom ', if(j==6){'elongation'}else{'dilation'}, '\n')
            display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
            backupchild = child[i,1:7]
          }
          child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
          if (j %in% debug) {
            display_sol(child[i,1:7], voxel, showscene = F, col='green')
            Sys.sleep(debug_time)
          }
        }
      }
    }
  }
  # enforce the boundaries (is it still really useful? the mutation operator should be good enough now)
  # child[,1:7] = pmax(child[,1:7], matrix(lowerBounds, nrow=nrow(child), ncol=7, byrow = T))
  # child[,1:7] = pmin(child[,1:7], matrix(upperBounds, nrow=nrow(child), ncol=7, byrow = T))
  if (verbose_return) {
    return(list(child=child, selected_operators=selected_operators))
  } else {
    return(child)
  }
}


mutation3Dprobs = function (parent_chromosome, varNo, objDim, hidDim,
                            lowerBounds, upperBounds, mprob, mum, probs_operators = c(rep(1/8,8),0))
{
  debug = NULL
  debug_plotnew = F
  debug_time = 5
  popSize = nrow(parent_chromosome)
  child <- parent_chromosome
  mutate_choice = cumsum(probs_operators)
  for (i in 1:popSize) {
    mutate_prob = which(probs_operators >= runif(1))[1]
    # move the center toward the best contact area
    if ((hidDim == 0) | (is.na(child[i,7+objDim+4]))) stop('only possible if using hidDim')
    if (mutate_prob == 1) {
      if (1 %in% debug) {
        cat('\nRecenter Translation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
        spheres3d(child[i, 1:3],radius=0.01, col='red')
      }
      u = c(
        child[i,6]/2*sin(child[i,4])*cos(-pi/2+child[i,5]),
        child[i,6]/2*cos(child[i,4])*cos(-pi/2+child[i,5]),
        sin(-pi/2+child[i,5])*child[i,6]/2
      )
      u = u / sqrt(sum(u^2))
      for (j in 1:3) {
        child[i, j] = child[i, j] + u[j]*child[i,7+objDim+4]
        # ensure that we still are within the right boundaries
        if (child[i,j] < lowerBounds[j]) {
          child[i,j] = lowerBounds[j]
        } else if (child[i,j] > upperBounds[j]) {
          child[i,j] = upperBounds[j]
        }
      }
      # for (j in 4:5) {
      #   child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
      # }
      if (1 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        spheres3d(child[i, 1:3],radius=0.01, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # symmetric flip along orient (only if using hidDim)
    if ((hidDim == 0) | (is.na(child[i,7+objDim+1]))) stop('only possible if using hidDim')
    if (mutate_prob == 2) {
      if (2 %in% debug) {
        cat('\nSymmetric Flip\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
      }
      r = child[i, 7]
      for (j in 1:3) {
        child[i, j] = child[i, j] + 2*r*child[i, 7+objDim+j]
        # ensure that we still are within the right boundaries
        if (child[i,j] < lowerBounds[j]) {
          child[i,j] = lowerBounds[j]
        } else if (child[i,j] > upperBounds[j]) {
          child[i,j] = upperBounds[j]
        }
      }   
      if (2 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # orient-based dilation (only if using hidDim)
    if ((hidDim == 0) | (is.na(child[i,7+objDim+1]))) stop('only possible if using hidDim')
    if (mutate_prob == 3) {
      if (3 %in% debug) {
        cat('\nOrient-based dilation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
      }
      r = child[i, 7]
      r_new = mutate(r, lowerBounds[7], upperBounds[7], 1) # 1 instead of mum, to really explore far away
      child[i, 7] = r_new
      for (j in 1:3) {
        child[i, j] = child[i, j] + (r-r_new)*child[i, 7+objDim+j]
        # ensure that we still are within the right boundaries
        if (child[i,j] < lowerBounds[j]) {
          child[i,j] = lowerBounds[j]
        } else if (child[i,j] > upperBounds[j]) {
          child[i,j] = upperBounds[j]
        }
      }
      if (3 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # random translation
    if (mutate_prob == 4) {
      if (4 %in% debug) {
        cat('\nRandom translation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
      }
      for (j in 1:3) {
        child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
      }
      if (4 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # random rotation
    if (mutate_prob == 5) {
      if (5 %in% debug) {
        cat('\nRandom rotation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
      }
      for (j in 4:5) {
        child[i, j] = mutate(child[i, j], lowerBounds[j], upperBounds[j], mum)
      }
      if (5 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # random elongation
    if (mutate_prob == 6) {
      if (6 %in% debug) {
        cat('\nRandom elongation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
        backupchild = child[i,1:7]
      }
      child[i, 6] = mutate(child[i, 6], lowerBounds[6], upperBounds[6], mum)
      if (6 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    # random dilation
    if (mutate_prob == 7) {
      if (7 %in% debug) {
        cat('\nRandom dilation\n')
        display_sol(child[i,1:7], voxel, showscene = debug_plotnew, col='red')
        backupchild = child[i,1:7]
      }
      child[i, 7] = mutate(child[i, 7], lowerBounds[7], upperBounds[7], mum)
      if (7 %in% debug) {
        display_sol(child[i,1:7], voxel, showscene = F, col='green')
        Sys.sleep(debug_time)
      }
    }
    
    ## TODO: symmetric rotation along the orient vector
    #     child = indiv #garbage
    #     u = apply(centers,2,diff)
    #     u = u / sqrt(sum(u^2))
    #     center_middle = apply(centers,2,mean)
    #     rot_axis1 = center_middle + u * proj_length
    #     rot_axis2 = center_middle + u * proj_length + orient * r
    #     rot_axis = rot_axis2 - rot_axis1
    #     
    #     xy_over_z = sqrt(sum((r*rot_axis[1:2])^2)) / abs(r*rot_axis[3])
    #     angle = pi/2
    #     
    #     rotindiv = c(rot_axis1+rotate3d(-indiv[1:3]+rot_axis1, angle, rot_axis[1], rot_axis[2], rot_axis[3]),
    #                  (indiv[4]+angle*(1-xy_over_z)) %% (2*pi),
    #                  (indiv[5]+angle*xy_over_z) %% pi,
    #                  indiv[6:7])
    #     display_sol(indiv,showscene=T,voxel = voxel, col='blue')
    #     display_sol(rotindiv,showscene=F,voxel = voxel, col='black')
  }
  # enforce the boundaries (is it still really useful? the mutation operator should be good enough now)
  # child[,1:7] = pmax(child[,1:7], matrix(lowerBounds, nrow=nrow(child), ncol=7, byrow = T))
  # child[,1:7] = pmin(child[,1:7], matrix(upperBounds, nrow=nrow(child), ncol=7, byrow = T))
  return(child)
}

Xover3D = function (parent_chromosome, varNo, objDim, hidDim,
                    lowerBounds, upperBounds, cprob, mu)
{
  popSize = nrow(parent_chromosome)
  child <- parent_chromosome
  p <- 1
  debug = F
  for (i in 1:(popSize/2)) {
    possible_alterations = list(xyz=c(1,2,3), orient=c(4,5), l=6, r=7)
    selected_alterations = lapply(which(runif(length(possible_alterations)) < cprob), function(il)possible_alterations[[il]])
    for (target in seq_along(selected_alterations)) {
      rnd2 = runif(1)
      xoverdir = sample.int(2)[1] # which solution should be the base (for x-over)?
      if (debug) {
        display_sol(child[p,1:7], voxel, showscene = T, col=grey(0.25))
        display_sol(child[p+1,1:7], voxel, showscene = F, col=grey(0.75))
        cat('X-over on ',names(possible_alterations)[target],' with parameter ',rnd2,'\n')
      }
      for (j in selected_alterations[[target]]) {
        if (xoverdir == 1) {
          parent2 <- child[p + 1, j]
          parent1 <- child[p, j]
        } else {
          parent2 <- child[p, j]
          parent1 <- child[p + 1, j]
        }
        yl <- lowerBounds[j]
        yu <- upperBounds[j]
        if (abs(parent1 - parent2) > 1e-06) {
          if (parent2 > parent1) {
            y2 <- parent2
            y1 <- parent1
          } else {
            y2 <- parent1
            y1 <- parent2
          }
          if ((y1 - yl) > (yu - y2)) {
            beta = 1 + (2 * (yu - y2)/(y2 - y1))
          } else {
            beta = 1 + (2 * (y1 - yl)/(y2 - y1))
          }
          alpha = 2 - (beta^(-(1 + mu)))
          if (rnd2 <= 1/alpha) {
            alpha = alpha * rnd2
            betaq = alpha^(1/(1 + mu))
          } else {
            alpha = alpha * rnd2
            alpha = 1/(2 - alpha)
            betaq = alpha^(1/(1 + mu))
          }
          child1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
          child2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
        } else {
          betaq = 1
          y1 = parent1
          y2 = parent2
          child1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
          child2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
        }
        if (child1 > yu) {
          if (debug) cat('Out of bound X-over\n')
          child1 = yu
        } else if (child1 < yl) {
          if (debug) cat('Out of bound X-over\n')
          child1 = yl
        }
        if (child2 > yu) {
          if (debug) cat('Out of bound X-over\n')
          child2 = yu
        } else if (child2 < yl) {
          if (debug) cat('Out of bound X-over\n')
          child2 = yl
        }
        child[p, j] = child1
        child[p + 1, j] = child2
      }
      if (debug) {
        display_sol(child[p,1:7], voxel, showscene = F, col='brown')
        display_sol(child[p+1,1:7], voxel, showscene = F, col='yellow')
        # Sys.sleep(6)
        browser()
      }
    }
    p <- p + 2
  }
  return(child)
}


boundedSBXover = function (parent_chromosome, varNo, objDim, hidDim,
                           lowerBounds, upperBounds, cprob, mu) 
{
  popSize = nrow(parent_chromosome)
  child <- parent_chromosome
  p <- 1
  for (i in 1:(popSize/2)) {
    if (runif(1) < cprob) {
      for (j in 1:varNo) {
        parent1 <- child[p, j]
        parent2 <- child[p + 1, j]
        yl <- lowerBounds[j]
        yu <- upperBounds[j]
        rnd = runif(1)
        if (rnd <= 0.5) {
          if (abs(parent1 - parent2) > 1e-06) {
            if (parent2 > parent1) {
              y2 <- parent2
              y1 <- parent1
            }
            else {
              y2 <- parent1
              y1 <- parent2
            }
            if ((y1 - yl) > (yu - y2)) {
              beta = 1 + (2 * (yu - y2)/(y2 - y1))
            }
            else {
              beta = 1 + (2 * (y1 - yl)/(y2 - y1))
            }
            alpha = 2 - (beta^(-(1 + mu)))
            rnd = runif(1)
            if (rnd <= 1/alpha) {
              alpha = alpha * rnd
              betaq = alpha^(1/(1 + mu))
            }
            else {
              alpha = alpha * rnd
              alpha = 1/(2 - alpha)
              betaq = alpha^(1/(1 + mu))
            }
            child1 = 0.5 * ((y1 + y2) - betaq * (y2 - 
                                                   y1))
            child2 = 0.5 * ((y1 + y2) + betaq * (y2 - 
                                                   y1))
          }
          else {
            betaq = 1
            y1 = parent1
            y2 = parent2
            child1 = 0.5 * ((y1 + y2) - betaq * (y2 - 
                                                   y1))
            child2 = 0.5 * ((y1 + y2) + betaq * (y2 - 
                                                   y1))
          }
          if (child1 > yu) {
            child1 = yu
          }
          else if (child1 < yl) {
            child1 = yl
          }
          if (child2 > yu) {
            child2 = yu
          }
          else if (child2 < yl) {
            child2 = yl
          }
        }
        else {
          child1 = parent1
          child2 = parent2
        }
        child[p, j] <- child1
        child[p + 1, j] <- child2
      }
    }
    p <- p + 2
  }
  return(child)
}

## This function implements the termination criterion. Returns True to continue, False to terminate.
continue_search = function(fitnesses, ite, min_ite, max_ite, constant_ite_termination, constant_tolerance=1e-6)
{
  # do at least min_ite and do at least constant_ite_termination
  if (ite <= max(c(min_ite, constant_ite_termination))) {
    return(T)
  }
  # stop after max_ite
  if (ite > max_ite) {
    message(paste0('Stopping search at iteration ',ite,': <max_ite> = ',max_ite,' has been reached.\n'))
    return(F)
  }
  for (i in (ite-constant_ite_termination):(ite-1)) {
    if (abs(fitnesses[i] - fitnesses[ite]) > constant_tolerance) {
      return(T)
    }
  }
  message(paste0('Stopping search at iteration ',ite,': score has not improved in <constant_ite_termination> = ',constant_ite_termination,'.\n'))
  return(F)
}



crowdingDist4frnt = function (pop, rnk, rng, varNo, objDim) 
{
  popSize <- nrow(pop)
  # objDim <- length(rng)
  # varNo <- ncol(pop) - 1 - length(rng) # now passed as an argument (to handle the hidden variable)
  cd <- matrix(Inf, nrow = popSize, ncol = objDim)
  for (i in 1:length(rnk)) {
    selectRow <- pop[, ncol(pop)] == i
    len <- length(rnk[[i]])
    if (len > 2) {
      for (j in 1:objDim) {
        originalIdx <- rnk[[i]][order(pop[selectRow, 
                                          varNo + j])]
        cd[originalIdx[2:(len - 1)], j] = abs(pop[originalIdx[3:len], 
                                                  varNo + j] - pop[originalIdx[1:(len - 2)], 
                                                                   varNo + j])/rng[j]
      }
    }
  }
  return(cd)
}


fastNonDominatedSorting = function (inputData) 
{
  popSize = nrow(inputData)
  idxDominators = vector("list", popSize)
  idxDominatees = vector("list", popSize)
  for (i in 1:(popSize - 1)) {
    for (j in i:popSize) {
      if (i != j) {
        xi = inputData[i, ]
        xj = inputData[j, ]
        if (all(xi <= xj) && any(xi < xj)) {
          idxDominators[[j]] = c(idxDominators[[j]], 
                                 i)
          idxDominatees[[i]] = c(idxDominatees[[i]], 
                                 j)
        }
        else if (all(xj <= xi) && any(xj < xi)) {
          idxDominators[[i]] = c(idxDominators[[i]], 
                                 j)
          idxDominatees[[j]] = c(idxDominatees[[j]], 
                                 i)
        }
      }
    }
  }
  noDominators <- lapply(idxDominators, length)
  rnkList <- list()
  rnkList <- c(rnkList, list(which(noDominators == 0)))
  solAssigned <- c()
  solAssigned <- c(solAssigned, length(which(noDominators == 
                                               0)))
  while (sum(solAssigned) < popSize) {
    Q <- c()
    noSolInCurrFrnt <- solAssigned[length(solAssigned)]
    for (i in 1:noSolInCurrFrnt) {
      solIdx <- rnkList[[length(rnkList)]][i]
      hisDominatees <- idxDominatees[[solIdx]]
      for (i in hisDominatees) {
        noDominators[[i]] <- noDominators[[i]] - 1
        if (noDominators[[i]] == 0) {
          Q <- c(Q, i)
        }
      }
    }
    rnkList <- c(rnkList, list(sort(Q)))
    solAssigned <- c(solAssigned, length(Q))
  }
  return(rnkList)
}

tournamentSelection = function (pop, pool_size, tour_size) 
{
  popSize = nrow(pop)
  Dim = ncol(pop)
  f <- NULL
  counter <- 1
  while (counter <= pool_size) {
    candidate = sample(popSize, tour_size)
    tmp <- pop[candidate, ]
    f <- rbind(f, tmp[order(tmp[, Dim - 1], -tmp[, Dim])[1], 
                      ])
    counter <- counter + 1
  }
  return(f)
}

nsga2R = function (fn, 
                   varNo, objDim,
                   sort_after_eval=0, decreasing_sort=F, div_fn=NULL,
                   hidDim = 0, 
                   crossover=boundedSBXover, mutation=boundedPolyMutation,
                   showBest = 0,
                   lowerBounds = rep(-Inf, varNo), upperBounds = rep(Inf, varNo), popSize = 100, tourSize = 2, 
                   generations = 20, cprob = 0.7, XoverDistIdx = 5, mprob = 0.2, 
                   MuDistIdx = 10, initialPop = NULL, verbose=T,
                   dosleeps=F, # this triggers a very small pause every generation, useful for visualization (to process 3D calls for example)
                   ...) 
{
  if (showBest) {
    n = ceiling(generations/showBest)+1
    cm = rev(rainbow(n*1.5)[(n/2):(1.5*n)])
    icm = 0
    plot3d(points[which(label_voxels %in% voxels_isaboveground[which(segmented_voxels[,4]==voxel)]),], col = 'black', aspect=T)
  }
  if (verbose)
    cat("********** NSGA2R Custom Version Implementing 3D Operators *********\ninitializing the population... ")
  if (is.null(initialPop)) {
    parent <- t(sapply(1:popSize, function(u) array(runif(length(lowerBounds), 
                                                          lowerBounds, upperBounds))))
    # parent <- cbind(parent, t(apply(parent, 1, function(x) fn(x,...)))) # doesn't work with one obj
  } else {
    parent = initialPop
  }
  # enforce boundaries
  parent <- cbind(parent, matrix(apply(parent, 1, function(x) fn(x,...)), nrow=popSize, byrow=T) )
  if (sort_after_eval != 0) {
    parent = parent[order(parent[,sort_after_eval], decreasing = decreasing_sort),]
  }
  if (!is.null(div_fn)) {
    if (verbose)
      cat("evaluating diversity objective... ")
    # the function takes into input (a) the matrix of indivs sorted for the first objective, and (b) the vector of first objective values
    parent[,varNo+objDim] = div_fn(parent[,1:varNo], parent[,varNo+1]) 
  }
  if (verbose)
    cat("ranking the initial population... ")
  ranking <- fastNonDominatedSorting(matrix(parent[, (varNo + 1):(varNo + objDim)], nrow=nrow(parent), byrow = F))
  rnkIndex <- integer(popSize)
  i <- 1
  while (i <= length(ranking)) {
    rnkIndex[ranking[[i]]] <- i
    i <- i + 1
  }
  parent <- cbind(parent, rnkIndex)
  if (verbose)
    cat("crowding distance calculation...\n")
  objRange <- apply(matrix(parent[, (varNo + 1):(varNo + objDim)],nrow=nrow(parent), byrow = F), 2, max) -
    apply(matrix(parent[, (varNo + 1):(varNo + objDim)],nrow=nrow(parent), byrow = F), 2, min)
  cd <- crowdingDist4frnt(parent, ranking, objRange, varNo, objDim)
  parent <- cbind(parent, apply(cd, 1, sum))
  inter.handler = function(interrupt){
    browser()
    cat('ended with the interrupt: ',str(interrupt),'\n')
  }
  withCallingHandlers(
    tryCatch({
      iter = 0
      while (iter < generations) {
        iter = iter + 1
        if (verbose)
          cat("Generation ", iter, "starts: \ntournament selection... ")
        childAfterX = NULL
        while (is.null(childAfterX) || (nrow(childAfterX) < popSize)) { # ensure that we have unique childrens
          if (is.null(childAfterX)) {
            poolMat = tournamentSelection(parent, popSize, tourSize)
          } else {
            poolMat = rbind(tournamentSelection(parent, popSize-nrow(childAfterX), tourSize), childAfterX)
          }
          if (verbose)
            cat("mutation operator... ")
          childAfterM <- mutation(poolMat, varNo, objDim, hidDim, lowerBounds, upperBounds, mprob, MuDistIdx)
          if (verbose)
            cat("crossover operator... ")
          childAfterX <- crossover(childAfterM, varNo, objDim, hidDim, lowerBounds, upperBounds, cprob, XoverDistIdx)
          # remove duplicated elements
          childAfterX = childAfterX[!duplicated(rbind(parent[,1:varNo],childAfterX[,1:varNo]))[(popSize+1):(popSize+nrow(childAfterX))],]
        }
        if (verbose)
          cat("evaluate the objective fns of childAfterM... ")
        # childAfterM <- cbind(childAfterM, t(apply(childAfterM, 1, function(x) fn(x,...) ))) # see above
        childAfterX <- cbind(childAfterX[,1:varNo], matrix(apply(childAfterX[,1:varNo], 1, function(x) fn(x,...) ), nrow=nrow(childAfterX), byrow = T)) 
        if (verbose)
          cat("Rt = Pt + Qt... ")
        parentNext <- rbind(parent[,1:(varNo+objDim+hidDim)], childAfterX)
        if (sort_after_eval != 0) {
          parentNext = parentNext[order(parentNext[,sort_after_eval], decreasing = decreasing_sort),]
        }
        if (!is.null(div_fn)) {
          if (verbose)
            cat("evaluating diversity objective... ")
          # the function takes into input (a) the matrix of indivs sorted for the first objective, and (b) the vector of first objective values
          parentNext[,varNo+objDim] = div_fn(parentNext[,1:varNo], parentNext[,varNo+1]) 
        }
        if (verbose)
          cat("ranking again... ")
        ranking <- fastNonDominatedSorting(matrix(parentNext[, (varNo + 1):(varNo + objDim)], nrow(parentNext), byrow = F))
        i <- 1
        while (i <= length(ranking)) {
          rnkIndex[ranking[[i]]] <- i
          i <- i + 1
        }
        parentNext <- cbind(parentNext, rnkIndex)
        if (verbose)
          cat("crowded comparison again... ")
        objRange <- apply( matrix(parentNext[, (varNo + 1):(varNo + objDim)], nrow=nrow(parentNext), byrow = F), 2, max) - 
          apply( matrix(parentNext[, (varNo + 1):(varNo + objDim)], nrow=nrow(parentNext), byrow = F), 2, min)
        cd <- crowdingDist4frnt(parentNext, ranking, objRange, varNo, objDim)
        parentNext <- cbind(parentNext, apply(cd, 1, sum))
        parentNext.sort <- parentNext[order(parentNext[, varNo + objDim + 1 + hidDim], 
                                            -parentNext[, varNo + objDim + 2 + hidDim]), ]
        if (verbose)
          cat("environmental selection\n")
        parent <- parentNext.sort[1:popSize, ]
        if ((showBest) & (iter %% showBest == 0)) {
          icm = icm + 1
          display_sol(parent[which.min(parent[,varNo+1]),1:varNo], voxel, showscene = F, col=cm[icm])
        }
        if (dosleeps) {
          Sys.sleep(0.001)
        }
      }
    },  interrupt = inter.handler
    ) # end tryCatch
  ) # end withCallingHandlers
  if (hidDim != 0) {
    hiddenVariables = parent[, (varNo + objDim + 1):(varNo + objDim + hidDim)]
  } else {
    hiddenVariables = NA
  }
  result = list(functions = fn, argfunctions = ..., parameterDim = varNo, objectiveDim = objDim, 
                lowerBounds = lowerBounds, upperBounds = upperBounds, 
                popSize = popSize, tournamentSize = tourSize, generations = generations, 
                XoverProb = cprob, XoverDistIndex = XoverDistIdx, mutationProb = mprob, 
                mutationDistIndex = MuDistIdx, 
                parameters = parent[, 1:varNo], 
                objectives = parent[, (varNo + 1):(varNo + objDim)], 
                hiddenVariables = hiddenVariables, 
                paretoFrontRank = parent[, varNo + objDim + hidDim + 1],
                crowdingDistance = parent[, varNo + objDim + hidDim + 2],
                iter = iter
  )
  class(result) = "nsga2R"
  return(result)
}

## LIST OF IMPROVEMENTS
# - handles single objective NSGA
# - handles interrupts
# - handles custom operators (specified in the function call)
# - handles hidden dimensions
# - add "iter" to the returned list


# this functions ranks the solutions and select the optimal ones
get_ranks = function(scores_cpp, acceptance_ratio, minPopSize, maxPopSize, Solutions, maxdist=NULL)
{
  ranks = list()
  
  if (is.null(maxdist)) {
    # keep all "good enough" solutions + minPopSize/2 other solutions
    new_order_cpp = order(scores_cpp, decreasing = T)
    kept_solutions_cpp = sum(scores_cpp > acceptance_ratio, na.rm=T)
    newPopSize_cpp = minPopSize + 2*kept_solutions_cpp # new popSize is twice the number of "good enough" solutions
    kept_solutions_cpp = kept_solutions_cpp + round(minPopSize/2) # keep at least minPopSize/2 "failed" solutions (in addition to the satisfactory solutions)
    kept_solutions_cpp = min(kept_solutions_cpp, maxPopSize)
  } else {
    # keep all "good enough" solutions, and all the other solutions that are farther away than maxdist
    new_order_cpp = order(scores_cpp, decreasing = T)
    kept_solutions_cpp = sum(scores_cpp > acceptance_ratio, na.rm=T)
    if (kept_solutions_cpp == 0) {
      kept_solutions_cpp = 1
    }
    kept_i = 1:kept_solutions_cpp
    dmat = as.matrix(dist(Solutions[new_order_cpp[1:nrow(Solutions)], 1:3]))
    if (sum(scores_cpp != 0) > 1) {
      # bypass adding solutions if their score is null
      for (i in (kept_solutions_cpp+1):sum(scores_cpp != 0)) {
        # keep more solutions if they are far enough, unless they are too crappy
        if ((min(dmat[i, kept_i]) > maxdist) &  (scores_cpp[new_order_cpp[i]] >= acceptance_ratio/10)) {
          # we add that far away solution
          kept_i = c(kept_i, i)
        }
      }
    }
    new_order_cpp = new_order_cpp[c(kept_i, (1:length(new_order_cpp))[-kept_i])]
    kept_solutions_cpp = length(kept_i)
    newPopSize_cpp = minPopSize + 2*kept_solutions_cpp # new popSize is twice the number of "good enough"+"far enough" solutions
  }
  
  # enforce min and max population size constraints:
  newPopSize_cpp = min(c(newPopSize_cpp, maxPopSize))
  newPopSize_cpp = max(c(newPopSize_cpp, minPopSize))
  
  ranks$order = new_order_cpp
  ranks$newPopSize = newPopSize_cpp
  ranks$keptSize = kept_solutions_cpp
  return(ranks)
}


