import os
import subprocess
import commands
import re
from glob import glob
import sys
from shutil import move
test = False

year = "2018"

#filelist = ["1.584_207"]
filelist = ["5.126_325","5.177_326","5.229_327","5.282_328","5.334_329","5.388_330","5.442_331","5.496_332"]
for mass in filelist:    
    
    files = glob("output_highpvd/dpCard_"+year+"IterV3_m"+mass+".txt")
    for fname in files:
    	m = re.search(r"_m(\d+\.?\d+)_", fname).group(1)
    	d = re.search(r"_(\d+).txt", fname).group(1)
    	print "ID {0},  m={1}".format(d,m)
    
    	if test and (not int(d)%10==0): continue
    	if (os.path.isfile(fname)):
            os.system("combineCards.py output_highpvd/dpCard_2017IterV3_m"+m+"_"+d+".txt output_dual/dpCard_2018IterV3_m"+m+"_"+d+".txt > output_dual/dpCard_IterV3_m"+m+"_"+d+".txt")
#            os.system("combine -M Significance --pvalue output_dual/dpCard_"+year+"IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+ " --cminDefaultMinimizerStrategy 1  --setParameters pdf_index_2017=0,pdf_index_2018=0  --rMin -0.1 --rMax 1")
#            os.system("combine -M AsymptoticLimits dual/dkpi_voig/output_dual/dpCard_IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+" --setParameters pdf_index_2017=0,pdf_index_2018=0 --rAbsAcc=0.0001 --rRelAcc=0.0001 --cminDefaultMinimizerStrategy 0 --rMin -1 --rMax 1")
            os.system("combine -M AsymptoticLimits output_highpvd/dpCard_IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+" --setParameters pdf_index_2017=0,pdf_index_2018=0 --rAbsAcc=0.0001 --rRelAcc=0.0001 --cminDefaultMinimizerStrategy 0 --rMin -1 --rMax 1")
 
outfiles = glob("higgsCombineasympMassIndex_*.root")
#    for of in outfiles:	
#    	move(of, "combine_output/bothYears/")
