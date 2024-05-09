import os
import subprocess
import commands
import re
from glob import glob
import sys
from shutil import move
test = False

year = "2018"

folder_path = 'output_dual'
file_list = os.listdir(folder_path)

filelist = []

for filename in file_list:
    if filename.endswith('.txt'):
        match = re.search(r'dpCard_2018IterV3_m(.+?)\.txt', filename)
        if match:
            content = match.group(1)
            filelist.append(content)

print(filelist)


#filelist = ["1.584_207"]
#filelist = ["1.508_202","1.569_206","1.616_209","1.682_213","1.523_203","1.584_207","1.632_210","1.538_204","1.649_211","1.553_205","1.600_208","1.665_212","1.682_213","1.699_214","1.716_215","1.733_216","1.750_217","1.768_218","1.785_219","1.803_220","1.821_221","1.839_222"]
for mass in filelist:    
    
    files = glob("output_dual/dpCard_"+year+"IterV3_m"+mass+".txt")
    for fname in files:
    	m = re.search(r"_m(\d+\.?\d+)_", fname).group(1)
    	d = re.search(r"_(\d+).txt", fname).group(1)
    	print "ID {0},  m={1}".format(d,m)
    
    	if test and (not int(d)%10==0): continue
    	if (os.path.isfile(fname)):
            os.system("combineCards.py output_dual/dpCard_2017IterV3_m"+m+"_"+d+".txt output_dual/dpCard_2018IterV3_m"+m+"_"+d+".txt > output_dual/dpCard_IterV3_m"+m+"_"+d+".txt")
#            os.system("combine -M Significance --pvalue output_dual/dpCard_"+year+"IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+ " --cminDefaultMinimizerStrategy 1  --setParameters pdf_index_2017=0,pdf_index_2018=0  --rMin -0.1 --rMax 1")
#            os.system("combine -M AsymptoticLimits dual/dkpi_voig/output_dual/dpCard_IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+" --setParameters pdf_index_2017=0,pdf_index_2018=0 --rAbsAcc=0.0001 --rRelAcc=0.0001 --cminDefaultMinimizerStrategy 0 --rMin -1 --rMax 1")
            os.system("combine -M AsymptoticLimits output_dual/dpCard_IterV3_m"+m+"_"+d+".txt -m "+m+" -n asympMassIndex_"+d+" --setParameters pdf_index_2017=0,pdf_index_2018=0 --rAbsAcc=0.0001 --rRelAcc=0.0001 --cminDefaultMinimizerStrategy 0 --rMin -1 --rMax 1")
 
outfiles = glob("higgsCombineasympMassIndex_*.root")

for of in outfiles:
        move(of, "combine_output/dual/")
#    for of in outfiles:	
#    	move(of, "combine_output/bothYears/")
