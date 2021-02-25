#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/CMSSW_10_3_1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/

root -l -b generateTree.C\(\"2018/MCDY2018LowMass/DarkPhoton/crab_20200706_185414/200706_165513/0000/merged.root\"\,\"out.root\"\,\"true\"\)
