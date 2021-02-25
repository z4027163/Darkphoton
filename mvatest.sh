#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/CMSSW_7_3_4/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

source /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/TMVA-v4.2.0/test/setup.sh
cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/
in=/$1
dir=$2
outdir=$2
out2="${in/tree_/tree_}"
#out="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/tree/"$3"/"$out2
out="/eos/user/w/wangz/darkphoton/mva/test/"$2"/"$out2
echo $out

root -l -b -q generate_test.C\(\"$2$in\"\,\"$out\"\)
