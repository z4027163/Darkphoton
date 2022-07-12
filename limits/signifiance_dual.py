import os
import subprocess
import commands
import re
from glob import glob
import sys
from shutil import move
test = False

#year = "2017"
year = ""
#mass = "1.584_207"
mass = "1.716_215"
#mass = "1.750_217"
if not os.path.isdir("combine_output/bothYears"):
        os.makedirs("combine_output/bothYears")

files = glob("output_dual/dpCard_"+year+"IterV3_m"+mass+".txt")
for fname in files:
	m = re.search(r"_m(\d+\.?\d+)_", fname).group(1)
	d = re.search(r"_(\d+).txt", fname).group(1)
	print "ID {0},  m={1}".format(d,m)

        paramLimits = ""

	if test and (not int(d)%10==0): continue
	if (os.path.isfile(fname)):
                #os.system("combineCards.py output_dual/dpCard_2017IterV3_m"+m+"_"+d+".txt output_dual/dpCard_2018IterV3_m"+m+"_"+d+".txt > output_dual/dpCard_IterV3_m"+m+"_"+d+".txt")
                os.system("combine -M Significance --pvalue output_dual/dpCard_"+year+"IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+ " --cminDefaultMinimizerStrategy 1  --setParameters pdf_index_2017=0,pdf_index_2018=0  --rMin -0.1 --rMax 1") 
outfiles = glob("higgsCombineasympMassIndex_*.root")
#for of in outfiles:	
#	move(of, "combine_output/bothYears/")
