#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/CMSSW_10_3_1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/
in=/$1
dir=$2
outdir=$3
out2="${in/scout_/tree_}"
#out="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/tree/"$3"/"$out2
out="/eos/user/w/wangz/darkphoton/"$3"/"$out2
echo $out

root -l -b -q generateTree_2017full.C\(\"$2$in\"\,\"$out\"\,false\)

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/CMSSW_7_3_4/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

source /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/TMVA-v4.2.0/test/setup.sh
cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/
in2=/$out2

root -l -b -q add_mva.C\(\"$3$in2\"\)

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/CMSSW_10_3_1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/

out3="${out2/tree_/hists_}"
out4="/eos/user/w/wangz/darkphoton/limit/"$3"/"$out3

root -l -b -q trimscoutFromfulltree.C\(\"$3$in2\"\,\"$out4\"\)
