# 3D Cylinder Estimation

 Added a uniform sampling module for polygons to repo by steph@snevets.institute.
 Credit for original repo: steph@snevets.institute

## Installation

1. Install [R](https://cran.r-project.org/)
2. Install [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
3. Install [RTools](https://cran.r-project.org/)
4. Run install_required.R from RStudio

I already had Visual Studio 2013 installed on my computer, but I do not think it is necessary to have this to compile as RTools appear to install MinGW things.

After you mentioned that Linux ran fine, I figure it's probably the MinGW compiler mucking up the speed of execution.


## Model Generation

Run Tree1.R from RStudio. This file has been modified so that it will output only the models from a single run of main_loop.R (5000 iterations) and then terminate.

The model generation script will ask for a 'Run Name'. Each run must be given a unique name. The output will be stored in a folder with that name relative to Tree1.R.

Model snapshots will be taken every 50 iterations.

When you're done running Tree1.R go ahead and restart the R session via the R Studio toolbar so the viewer doesn't draw things in black.


## Viewing

Run viewer.R from RStudio. This will take in the output generated in the model generation step and display it on screen.

The model generation script will ask for a 'Run Name'. The viewer will fetch the models stored in the folder with that name relative to viewer.R.

At present, it just converts the models to images, but if you uncomment line 86 of viewer.R it will allow you to interact with the model before the snapshot is taken (you'll have to type y to continue then).

You can view our five runs, or create your own. Try typing in `run1` when it prompts you for a run name, and watch it work.


## Solution Score Graphing

Run graph.R from RStudio. This will take in the output generated in the model generation step and create a graph of average solution score per frame.

The model generation script will ask for a 'Run Name'. The grapher will fetch the models stored in the folder with that name relative to graph.R. Unlike viewer.R, it will not save the images for you automatically because I am lazy and there were not literally 25000 images to save. 5 images were easier to do by hand than figure out how to code that in. So save the images yourself.

You can view our five runs, or create your own. Try typing in `run1` when it prompts you for a run name, and watch it work.


## Gifs


### Gif ToolChain

For Gif making, I assume you have ffmpeg installed and python 3 (I used python 3.6). If not, download ffmpeg [here](https://www.ffmpeg.org/) and python 3 [here](https://www.python.org/downloads/).


### Making the Gif

Assuming your run name is `run1` and you've run viewer.R

From the main folder:

	python gifreduce.py run1
	cd run1-temp
	ffmpeg -f image2 -i %04d.png output.gif
	
Which leaves you with a folder full of renamed pngs and output.gif. Feel free to delete the pngs.
