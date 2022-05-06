require('rgl')

# Initialization of search loop variables
fitnesses = rep(NA, max_ite)
best_fitness = rep(NA, max_ite)
ite = 0

while (continue_search(fitnesses, ite, min_ite, max_ite, constant_ite_termination)) {
  ite = ite + 1;
  Sys.sleep(0.001)
  
  if ('reduce_width' %in% ls() == T) {
    if (reduce_width) {
      lower[7] = 0.02 + (0.1-0.02) * 0.9*(max_ite-ite)/max_ite
    }
  }
  
  Scene_cpp$iterate()
  # cat('iterate done\n')
  
  scores_cpp = Scene_cpp$get_scores()
  Solutions = cbind(Scene_cpp$get_indivs(), Scene_cpp$get_scores(), Scene_cpp$get_hid_indivs())
  # cat('score done\n')
  
  ranks = get_ranks(scores_cpp, acceptance_ratio, minPopSize, maxPopSize, Solutions, maxdist = maxdist)
  
  # cat('cpp done\n')
  new_order_cpp = ranks$order
  newPopSize_cpp = ranks$newPopSize
  kept_solutions_cpp = ranks$keptSize
  
  # 3D plot of the solutions
  if  ((ite-1) %% ite_modulo_show3d == 0) {
    scene_plotted = F
    if (nrow(Points) > points_cutoff) {
      p_subsample = sample(1:nrow(Points), size = points_cutoff)
    } else {
      p_subsample = 1:nrow(Points)
    }
    if (('show_color_model' %in% ls() == F) || show_color_model) {
      if (('model_light' %in% ls() == F) || (!is.null(model_light))) {
        scene_plotted = T
        # plot3d(model_light, type='dots', alpha=0.5, top=F)
        clear3d(type='all') # we will plot the points after the cylinders (this is quicker)
      }
    }
    for (i in 1:min(max_3d_cylinders,popSize_cpp)) {
      sol_score = Scene_cpp$get_scores()[new_order_cpp[i]]
      if (sol_score >= acceptance_ratio) {
        scene_plotted = display_sol(Scene_cpp$get_indivs()[new_order_cpp[i],], Points[p_subsample,], 
                                    # col='orange', 
                                    col='purple', 
                                    showscene = !scene_plotted)
      } else if (sol_score > 0.05) { # arbitrary threshold to avoid displaying bad solutions
        scene_plotted = display_sol(Scene_cpp$get_indivs()[new_order_cpp[i],], Points[p_subsample,], 
                                    # col= 'purple', alpha = sol_score/acceptance_ratio, # transparent purple
                                    col= rainbow(100)[round(sol_score/acceptance_ratio*80)+1], # rainbow (quicker to draw)
                                    showscene = !scene_plotted)
      }
    }
    if (('show_bounding_box' %in% ls() == F) || show_bounding_box) {
      try({
        bounding_box = structure(list(x = c(1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0), y = c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1), z = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)), .Names = c("x", "y", "z"), row.names = c(NA, -24L), class = "data.frame")
        bounding_box = bounding_box * matrix(upper[1:3] - lower[1:3], ncol=3, nrow=24, byrow = T) + matrix(lower[1:3], ncol=3, nrow=24, byrow = T)
        segments3d(bounding_box, line_antialias = TRUE, col = "blue")
        bounding_box_points = structure(list(x = c(1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0), y = c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1), z = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)), .Names = c("x", "y", "z"), row.names = c(NA, -24L), class = "data.frame")
        bounding_box_points = bounding_box_points * matrix(upper[1:3] - lower[1:3] + max_length, ncol=3, nrow=24, byrow = T) + matrix(lower[1:3]-0.5*max_length, ncol=3, nrow=24, byrow = T)
        segments3d(bounding_box_points, line_antialias = TRUE, col = "green")
      })
    }
    if (('show_color_model' %in% ls() == F) || show_color_model) {
      if (('model_light' %in% ls() == F) || (!is.null(model_light))) {
        if ('vb' %in% names(model_light)) {
          # mesh3d object
          plot3d(model_light, type='dots', top=F, add=T)
        } else {
          # simple point list
          plot3d(model_light, top=F, add=T)
        }
      }
    }
  }
  
  # Adjust population size
  # cat('\nSuggesting a new population size of', newPopSize_cpp, 'for iteration',ite)
  if (newPopSize_cpp < popSize_cpp) {
    for (sol_i in (newPopSize_cpp+1):popSize_cpp) {
      Scene_cpp$remove_indiv(new_order_cpp[sol_i]-1) # remove sub-optimal solution
      cat('Removed solution from C++ because model is shrinking\n')
    }
  }
  popSize_cpp = newPopSize_cpp
  
  # stores the fitness
  fitnesses[ite] = sum(scores_cpp[new_order_cpp[1:kept_solutions_cpp]], na.rm=T) # sum scores of the kept solutions
  # fitnesses[ite] = sum(scores_cpp, na.rm=T) # sum all scores
  best_fitness[ite] = max(scores_cpp[new_order_cpp[1:kept_solutions_cpp]], na.rm=T) # sum scores of the kept solutions
  
  
  # select the solutions to change (stochastic process, has a probability to alter non-kept solutions based on their fitness)
  to_change_cpp = runif(min=0, max=acceptance_ratio, n = (popSize_cpp-randomOffsprings) - kept_solutions_cpp) >= 
    scores_cpp[new_order_cpp][(kept_solutions_cpp+1):(popSize_cpp-randomOffsprings)]
  to_change_cpp = which(to_change_cpp)
  if (length(to_change_cpp) == 0) {
    cat('No solution was selected for change, skipping this iteration...\n')
    next
  } else {
    # cat('Generating n=',length(to_change_cpp),' new solutions.\n')
  }
  
  # generates the offsprings
  Solutions = cbind(Scene_cpp$get_indivs(), Scene_cpp$get_scores(), Scene_cpp$get_hid_indivs())
  new_sol = fill_russian_roulette(n_to_fill=length(to_change_cpp),
                                  fitness= sapply(Scene_cpp$get_scores()[new_order_cpp[1:kept_solutions_cpp]], function(f)min(acceptance_ratio, f)),
                                  Solutions[new_order_cpp[1:kept_solutions_cpp],], 7, objDim, hidDim, lower, upper, mprob=0.25, cprob=0.5, mum=5, muc=5, m_c_ratio=0.8)
  # cat('generates the offsprings\n')
  
  # add them to the C++ model
  for (sol_i in 1:nrow(new_sol$childs)) {
    change_i = to_change_cpp[sol_i]
    if (new_sol$childs[sol_i, 6] < 0) {
      browser() # QAD => this is strange and should not happen
    }
    Scene_cpp$change_indiv(new_sol$childs[sol_i, 1:7], new_order_cpp[kept_solutions_cpp+change_i]-1)
  }
  # cat('add them to the C++ model\n')
  
  # generates random offspring (if enabled)
  if (randomOffsprings > 0) {
    if ((kept_solutions_cpp+length(to_change_cpp)+1) <= popSize_cpp) {
      indices_to_create = new_order_cpp[(kept_solutions_cpp+length(to_change_cpp)+1):popSize_cpp]
      for (i in indices_to_create) {
        Scene_cpp$change_indiv(generate_random_solution(lower, upper, objDim, hidDim, heuristic=NULL, Points=NULL)[1:7], i-1)
      }
    } else {
      message('New offsrping contraint not respected!\n')
    }
  }
  
  # barplot the fitness of the population
  if ((ite-1) %% ite_modulo_showbp == 0) {
    # barplot(t(Scene_cpp$get_scores())[new_order_cpp[1:popSize_cpp]], ylim=c(0, 1), las=1, main = paste('Iteration',ite, ' (population ',newPopSize_cpp,')'))
    col_barplot = rep('black', popSize_cpp)
    col_barplot[new_sol$parents] = 'red'
    # col_barplot[to_change_cpp+kept_solutions_cpp] = 'grey'
    # col_barplot[(kept_solutions_cpp+length(to_change_cpp)+1):popSize_cpp] = 'grey'
    barplot(t(Scene_cpp$get_scores())[new_order_cpp[1:popSize_cpp]], ylim=c(0, 1), las=1, main = paste('Iteration',ite, ' (population ',newPopSize_cpp,')'),
            border=col_barplot, col=col_barplot)
    abline(h=acceptance_ratio, lty=2)
  }
}
