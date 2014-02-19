

 source /uscmst1/prod/sw/cms/cshrc cmslpc
 scram project CMSSW $CMSSW
 cd $CMSSW/src/
 cmscvsroot CMSSW
 scram b
 eval `scramv1 runtime -csh`
