##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- parent.frame(2)$ofile
if (!is.null(this.dir)) {
  this.dir <- dirname(this.dir)
  setwd(this.dir)
}
print(this.dir)
library(jsonlite)
require('rgl')
source('toy_examples_fct.R')



closeAllConnections()

runname = readline(prompt="Enter a Run Name: ")
if(!file.exists(runname)){
  print("hey now that doesn't exist")
  print(runname)
  stop()
}
pngs = sort(list.files(path = runname, pattern = ".*png", all.files = FALSE,
                  full.names = TRUE, recursive = FALSE,
                  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE))
for(png in pngs){
  file.remove(png)
}
closeAllConnections()

files = list.files(path = runname, pattern = ".*txt", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
nums = unlist(lapply(files, function(x){
  t = strsplit(x,".", fixed = TRUE)[[1]][2]
  t = strtoi(t)
  return(t)
}))
onums = order(nums)
nums = nums[onums]
files = files[onums]

acceptance_ratio = 0.4
tmp = file.path(runname, "points.json")
MyPoints = fromJSON(readChar(tmp, file.info(tmp)$size))

for(fi in 1:length(files)){
  fname = files[fi]
  fileName = file.path(runname, fname)
  print(fileName)
  con = file(fileName, 'r')
  max_3d_cylinders = strtoi(readLines(con, n = 1))
  popSize_cpp = strtoi(readLines(con, n = 1))
  points_cutoff = strtoi(readLines(con, n = 1))
  if (nrow(MyPoints) > points_cutoff) {
    p_subsample = sample(1:nrow(MyPoints), size = points_cutoff)
  } else {
    p_subsample = 1:nrow(MyPoints)
  }
  #print(json)
  clear3d()
  
  for (i in 1:min(max_3d_cylinders,popSize_cpp)) {
    #[1] "B"
    #[1] 0.2
    #[1] -2.62224686 -0.27731394  1.49343946  4.30757416  2.60753989  0.08210138  0.03024792
    #[1] "B"
    sol_score = as.numeric(readLines(con, n = 1))
    indivs = lapply(strsplit(readLines(con, n = 1), " "), as.numeric)[[1]]
    if (sol_score >= acceptance_ratio) {
      scene_plotted = display_sol(indivs, MyPoints[p_subsample,], 
                                  # col='orange', 
                                  col='purple', 
                                  showscene = F)
      
    } else if (sol_score > 0.05) { # arbitrary threshold to avoid displaying bad solutions
      scene_plotted = display_sol(indivs, MyPoints[p_subsample,], 
                                  # col= 'purple', alpha = sol_score/acceptance_ratio, # transparent purple
                                  col= rainbow(100)[round(sol_score/acceptance_ratio*80)+1], # rainbow (quicker to draw)
                                  showscene = F)
    }
    
    
  }
  close(con)
  print(gettextf("%s/output.%04d.png", runname, nums[fi]))
  # continue = readline(prompt="continue? y/y")
  rgl.snapshot(gettextf("%s/output.%04d.png", runname, nums[fi]))
  
  #stop()
}