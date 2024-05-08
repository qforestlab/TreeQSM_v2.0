import os
import os.path
import glob
import numpy as np
import subprocess
import sys

files = glob.glob('pointcloud/*.txt')

for i in files:

    if not os.path.exists("results"):
        os.makedirs("results")
    
    filename = i.split(os.sep)[-1].split('.txt')[0]
    print('Starting %s' % filename)
        
    tracker = "results/"+filename+".tracker"
    if not os.path.isfile(tracker):
        subprocess.call(['touch', tracker])
        subprocess.call(['python', './runCylinderModel.py', '-i',i, '-m', '/usr/local/MATLAB/R2021b/bin/matlab', '-p', sys.argv[1]])
        subprocess.call(['python', './optimiseCylinderModel.py', '-c', 'pointcloud/', '-m', 'results'])