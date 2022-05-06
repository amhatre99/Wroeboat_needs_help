##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- parent.frame(2)$ofile
if (!is.null(this.dir)) {
  this.dir <- dirname(this.dir)
  setwd(this.dir)
}
print(this.dir)
library(jsonlite)
#require('rgl')
#source('toy_examples_fct.R')



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
average_solution_score = 1:length(files)
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
  su = 0
  tot = 0
  for (i in 1:min(max_3d_cylinders,popSize_cpp)) {
    #[1] "B"
    #[1] 0.2
    #[1] -2.62224686 -0.27731394  1.49343946  4.30757416  2.60753989  0.08210138  0.03024792
    #[1] "B"
    sol_score = as.numeric(readLines(con, n = 1))
    indivs = lapply(strsplit(readLines(con, n = 1), " "), as.numeric)[[1]]
    
    tot = tot + 1
    su = su + sol_score
  }
  average_solution_score[fi] = su / tot
  close(con)
  #print(gettextf("%s/output.%04d.png", runname, nums[fi]))
  # continue = readline(prompt="continue? y/y")
  #rgl.snapshot(gettextf("%s/output.%04d.png", runname, nums[fi]))
  
  #stop()
}
print(average_solution_score)
iteration = nums
plot(y = average_solution_score, x = iteration, main = runname)