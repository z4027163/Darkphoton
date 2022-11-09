#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

#cp -r /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/MVA/darkphoton/diffvar /tmp/wangz/

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/CMSSW_7_3_4/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

source /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/tnp/TMVA-v4.2.0/test/setup.sh
cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/
in=/$1
dir=$2
echo $out

root -l -b -q add_mva.C\(\"$2$in\"\,$3\)
