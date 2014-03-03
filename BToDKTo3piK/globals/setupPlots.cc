// sets up the global plots

#include "../BToDKTo3piK/globals/plots.hh"

void setupPlots() {
  cout << "--- START setupPlots() ---" << endl;

  mesFrame = mes->frame();
  DeltaeFrame = Deltae->frame();
  d0massFrame = d0mass->frame();
  nnFrame     = nnout->frame();
  bknnFrame   = bknnout->frame();
  R2Frame     = R2->frame();
  qprimeFrame = qprime->frame(-5,5,50);
  dprimeFrame = dprime->frame(-5,5,50);
  m12Frame    = m12->frame();
  m13Frame    = m13->frame();

  cout << "--- END setupPlots() ---" << endl;
}
