sourceFoundFile GeneratorsQA/EvtGen.tcl

mod talk GfiEvtGen
  DECAY set ../BToDKTo3piK/cleo.dec
exit

module talk PepBuildEnv  
    pepBoostCalFile  set  PepCond/pepBoostCal.raw  
    pepBeamSpotCalFile set PepCond/pepBeamSpotCal.raw 
    pepEnergiesFile set PepCond/pepEnergies.raw 
    pepEnergiesCorrFile set PepCond/pepEnergiesCorr.raw
    pepPackedBunchesFile set PepCond/pepPackedBunches.raw 
    pepFillPatternFile set PepCond/pepFillPattern.raw
exit

mod talk BtaLoadMcCandidates
  numEventsToPrint set 0
exit

mod talk RandomControl
  maxEventsPerRun set 20000
exit

ev begin -nev 20000
exit

