# August George - April 2020

import os
import subprocess
from time import sleep,time
from shutil import copy, copyfileobj

# configuration
trial = 5
run_length = 1e2
subdirectory = 'runs_%s_%s' %(int(run_length),trial)
config_file = 'config.txt'
cluster_folder = "cluster_analysis"
n = 50

root_dir = os.path.abspath(os.path.curdir)

start_time = time() 

### initialize  - run latin square once to get values
config_settings = "\nINITIALIZE_ONLY 1\nLOAD_MODEL 0\ninitial_model_n 0\ndirectory %s\nnsteps %s\n" % (subdirectory, run_length)
f = open(config_file, 'a')
f.write(config_settings)
f.close()


pipe = subprocess.Popen(["perl", "ModelExplorer.prl"])  # run perl program
print("intializing latin square...running perl...ModelExplorer \n")
pipe.communicate()

# go to 'runs' subdirectory 
os.chdir(os.path.join(os.path.abspath(os.path.curdir),subdirectory))  # now working in subdirectory

try:
    os.mkdir(cluster_folder)
    cluster_analysis_path = os.path.join(os.path.abspath(os.path.curdir),cluster_folder)
    copy("cluster_analyzer.py",cluster_analysis_path)
    f = open(os.path.join(cluster_analysis_path, "cluster_data_agg.dat"), "w")
    f.close()
except:
    print ("Error. Could not create %s and/or aggregate data file" % cluster_folder)



### run simulations from initial models
for i in range(n):
    print(i)
    config_settings = "\nINITIALIZE_ONLY 0\nLOAD_MODEL 1\ninitial_model_n %s\ndirectory run_%s\n" % (i,i)
    
    f = open(config_file, 'a')
    f.write(config_settings)
    f.close()
    
    print(config_settings)

    print(os.getcwd())
    
    pipe = subprocess.Popen(["perl", "ModelExplorer.prl"])  # run perl program
    print("running perl...\n")
    pipe.communicate()
    
    # copy cluster files into one folder

    cluster_file = "cluster_data.dat"
    new_cluster_file = "cluster_data_%s.dat" % i

    src_dir=os.path.join(os.path.join(os.path.abspath(os.path.curdir), "run_%s" %i), cluster_file)
    dst_dir=os.path.join(cluster_analysis_path,new_cluster_file)
    copy(src_dir,dst_dir)
    
    # append cluster data to one aggregate file
    with open(os.path.join(cluster_analysis_path, "cluster_data_agg.dat"), 'a') as outfile:
            with open(src_dir) as infile:
                for line in infile:
                    outfile.write(line)

runtime = time()-start_time  # testing for runtime data
print("%s s" % runtime)

### clustering analysis
# change working directory to analysis 
os.chdir(cluster_analysis_path)

pipe = subprocess.Popen(["python",  "cluster_analyzer.py"])  # run perl program
print("running python clustering script...\n")
pipe.communicate()

  



