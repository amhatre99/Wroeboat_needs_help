##steph: https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source('fn_pointcloud.R')
model_nocolor = color.read.ply('Tree4_DB_1Hz_30kkey_3ktie_high_cropped_0.005.ply', ShowSpecimen = F)
model_color = color.read.ply('Tree4_DB_1Hz_30kkey_3ktie_high_cropped.ply', ShowSpecimen = F)

require('FNN')
corresp = get.knnx(data=t(model_color$vb[1:3,]), query=t(model_nocolor$vb[1:3,]), k=1)
model_nocolor$material=list()
model_nocolor$material[['color']] = model_color$material$color[corresp$nn.index]
plot3d(model_nocolor, type='dots')
saveRDS(model_nocolor, 'Tree4_DB_1Hz_30kkey_3ktie_high_cropped_0.005.rds')
