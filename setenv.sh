#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd /afs/cern.ch/user/w/wangz/YOURWORKINGAREA/CMSSW_10_3_1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd /afs/cern.ch/work/w/wangz/DarkPhotonAnalysisV2/DimuonAnalysis/macros/

