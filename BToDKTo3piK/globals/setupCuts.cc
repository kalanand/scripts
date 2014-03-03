// sets up the global cuts

#include "../BToDKTo3piK/globals/cuts.hh"

// Setup composite cuts
// If you change any of the basic cuts in cuts.hh rerun setupCuts() !
void setupCuts() {

  cutNN = cutNNq+cutNNd;
  cutNN.SetName("cutNN");
  
  cutBasic = cutNN+cutDtoKpi+cutKsVeto;
  cutBasic.SetName("basic cut"); 
  
  cutSigReg = cutDeltaE && cutmES && cutMD && cutBasic;
  cutSigReg.SetName("signal region"); 


  cutSBUpperDE = TCut("0.06<Deltae&&Deltae<0.140") && cutmES && cutMD && cutBasic;
  cutSBUpperDE.SetName("upper #DeltaE SB");

  cutSBLowerDE = TCut("-0.140<Deltae&&Deltae<-0.07") && cutmES && cutMD && cutBasic;
  cutSBLowerDE.SetName("lower #DeltaE SB");

  cutSBmES = TCut("5.2<mes&&mes<5.272") && cutDeltaE && cutMD && cutBasic;
  cutSBmES.SetName("m_{ES} SB");

  cutSBUpperMD = TCut("1.9<d0mass&&d0mass") && cutDeltaE && cutmES && cutBasic;
  cutSBUpperMD.SetName("upper m_{D} SB");

  cutSBLowerMD = TCut("d0mass&&d0mass<1.82") && cutDeltaE && cutmES && cutBasic;
  cutSBLowerMD.SetName("lower m_{D} SB");

  cutSBMD = cutSBUpperMD || cutSBLowerMD;
  cutSBMD.SetName("m_{D} SB");



  // The cut actually used in reading files. Responsibility of reader to 
  // change appropriately before reading.
  readCut = cutSigReg;
}
