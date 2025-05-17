# submit to run over all files of a generator
import sys
import os 
import glob


input_dir = sys.argv[1].rstrip('/')

files = glob.glob(input_dir + "/miniAOD*.root")

files = ["file:" + f for f in files]
files = ",".join(files)

i=0
for j in range(10):
  if not os.path.exists('output_{}.root'.format(j)):
    break
  i += 1

cmd = "cmsRun $CMSSW_BASE/src/pOAnalysis/Analyzer/python/ConfFile_cfg.py inputFiles={} applyFilt=False".format(files)
print(cmd)
os.system(cmd)

if os.path.exists("OUTPUT/output_{}.root".format(input_dir.split("/")[-1].split("_")[-2] + "_"+input_dir.split("/")[-1].split("_")[-1])):
    os.system('rm OUTPUT/output_{}.root'.format(input_dir.split("/")[-1].split("_")[-2] + "_"+input_dir.split("/")[-1].split("_")[-1]))
os.system('mv output_{}.root OUTPUT/output_tracks_{}.root'.format(i,input_dir.split("/")[-1].split("_")[-2] + "_"+input_dir.split("/")[-1].split("_")[-1]))
