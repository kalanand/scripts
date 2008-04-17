#!/bin/tcsh
source /uscmst1/prod/sw/cms/cshrc prod
source /uscmst1/prod/grid/gLite_SL5.csh
source /uscmst1/prod/grid/CRAB/crab.csh
eval `scramv1 runtime -csh`
project CMSSW
cmscvsroot CMSSW
alias cvsd cvs -d $CVSROOT
alias grid voms-proxy-init -voms cms
