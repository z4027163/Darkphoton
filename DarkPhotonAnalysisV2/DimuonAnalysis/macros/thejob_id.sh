#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/CMSSW_10_3_1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/

in=/$1
dir=$2
outdir=$2
out2="${in/tree_/hists_}"
#out="/afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/tree/"$3"/"$out2
out="/eos/user/w/wangz/darkphoton/limit_id/"$2"/"$out2
echo $out

root -l -b -q trimscoutFromfulltree_ID.C\(\"$2$in\"\,\"$out\"\)
