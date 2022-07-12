import os
import subprocess as subprocess
import commands
import re
from glob import glob
import sys
from shutil import move
test = False

INDIR = "/afs/cern.ch/user/c/ccosby/lowMassZPrime/biasStudyLimit/CMSSW_10_3_2/src/DarkPhotonLimits/DimuonAnalysis2017/limits/"
OUTDIR = INDIR+"combine_output/bothYears/"
nsubjobs = 4

shfile = open("condor_combine_task.sh", "w") 
shfile.write("#!/bin/sh \n")
shfile.write("ulimit -s unlimited \n")
shfile.write("set -e \n")
shfile.write("cd "+INDIR+" \n")
shfile.write("export SCRAM_ARCH=slc7_amd64_gcc700 \n")
shfile.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n")
shfile.write("eval `scramv1 runtime -sh` \n")
shfile.write("cd - \n")
shfile.write("\n")


year = "2017"
files = glob("output/dpCard_"+year+"IterV3_*.txt")

job=0
subjob=0
for fname in files:
	m = re.search(r"_m(\d+\.?\d+)_", fname).group(1)
	d = re.search(r"_(\d+).txt", fname).group(1)
	print "ID {0},  m={1}".format(d,m)

        paramLimits = ""

        subjob+=1

	if test and (not int(d)%10==0): continue
	if os.path.isfile(fname):
                subprocess.call("combineCards.py output/dpCard_2017IterV3_m"+m+"_"+d+".txt output/dpCard_2018IterV3_m"+m+"_"+d+".txt > output/dpCard_IterV3_m"+m+"_"+d+".txt",shell=True)

                if (subjob==1):
                        shfile.write("if [ $1 -eq "+str(job)+" ]; then \n")
                        shfile.write("  mkdir output \n")
                shfile.write("  cp "+INDIR+"/output/dp*m"+m+"* output/ \n")
                shfile.write("  cp "+INDIR+"/output/dp*_"+d+".root output/ \n")
                shfile.write("  combine -M AsymptoticLimits output/dpCard_IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+" --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.01 --X-rtd MINIMIZER_freezeDisassociatedParams --rRelAcc 0.001 --rAbsAcc=0.001  --robustFit=1 \n")
                shfile.write("  cp higgsCombineasympMassIndex_*.root "+OUTDIR+" \n") 
                if (subjob==nsubjobs): shfile.write("fi \n")

                
        if (subjob==nsubjobs):
                subjob=0
                job+=1
                
#if (subjob<nsubjobs): shfile.write("fi \n")
shfile.close()

subfile = open("condor_combine_task.sub", "w") 
subfile.write("executable = condor_combine_task.sh \n")
subfile.write("arguments = $(ProcId) \n")
subfile.write("output                = combine_task.$(ClusterId).$(ProcId).out \n")
subfile.write("error                 = combine_task.$(ClusterId).$(ProcId).err \n")
subfile.write("log                   = combine_task.$(ClusterId).log \n")
subfile.write(" \n")
subfile.write("+JobFlavour = \"espresso\" \n")
subfile.write(" \n")
subfile.write("# Send the job to Held state on failure. \n")
subfile.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0) \n")
subfile.write(" \n")
subfile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries. \n")
subfile.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600) \n")
subfile.write(" \n")
subfile.write(" \n")
subfile.write("queue "+str(job)+" \n")

subfile.close()

subprocess.call("condor_submit condor_combine_task.sub",shell=True)
