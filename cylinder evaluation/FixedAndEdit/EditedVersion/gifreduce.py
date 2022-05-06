import sys
import os
import shutil

directory = sys.argv[1]
tempdir = directory + "-temp"
try:
    os.stat(tempdir)
except:
    os.mkdir(tempdir)

for i, f in enumerate(sorted([x for x in os.listdir(directory) if x.endswith(".png")])):
    print(f)
    shutil.copy(os.path.join(directory, f), os.path.join(tempdir,"%04d.png" % i))

    
