// #including this file sets up the global objects (see also utils/utils.cc)

using namespace RooFit;

#include "../BToDKTo3piK/globals/flags.hh"
#include "../BToDKTo3piK/globals/cuts.hh"

#include "../BToDKTo3piK/globals/setupCuts.cc"
#include "../BToDKTo3piK/globals/setupVars.cc"
#include "../BToDKTo3piK/globals/setupPdfHolder.cc"
#include "../BToDKTo3piK/globals/setupPdfs.cc"
#include "../BToDKTo3piK/globals/setupPlots.cc"
#include "../BToDKTo3piK/globals/setupChains.cc"
#include "../BToDKTo3piK/utils/utils.cc"
#include "../BToDKTo3piK/analysis/analysis.cc"
#include "../BToDKTo3piK/toy/toy.cc"


void setup(BdkPdfDKDalitz::COORD coord = BdkPdfDKDalitz::POLAR) {
  cout << "--- START setup() ---" << endl;

  setupCuts();
  setupVars();
  setupPdfHolder(coord, m12, m13);
  setupPdfs();
  setupChains();
  setupPlots();

  cout << "--- END setup() ---" << endl;
}
