// #inclued'ing this file sets up the global objects (see also utils/utils.cc)
// This is a reduced version of setup.cc only containing scripts that compile

#include "../BToDKTo3piK/globals/flags.hh"
#include "../BToDKTo3piK/globals/cuts.hh"
#include "../BToDKTo3piK/globals/setupCuts.cc"
#include "../BToDKTo3piK/globals/setupVars.cc"
#include "../BToDKTo3piK/globals/setupPdfHolder.cc"
#include "../BToDKTo3piK/globals/setupPdfs.cc"
#include "../BToDKTo3piK/globals/setupPlots.cc"
#include "../BToDKTo3piK/globals/setupChains.cc"

#include "../BToDKTo3piK/utils/TypeBits.hh"
#include "../BToDKTo3piK/utils/printEff.cc"

//#include "../BToDKTo3piK/utils/utils.cc"
//#include "../BToDKTo3piK/analysis/analysis.cc"
//#include "../BToDKTo3piK/toy/toy.cc"


void setupred() {
  cout << "--- START setupred() ---" << endl;

  setupCuts();
  setupVars();
  setupPdfHolder();
  setupPdfs();
  setupChains();
  //  setupPlots();

  cout << "--- END setupred() ---" << endl;
}
