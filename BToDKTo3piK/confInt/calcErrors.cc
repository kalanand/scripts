// Functions in this file:
//
// The two calcErrors() functions calculate the best value and the 
// errors for a 1D histogram. One of the functions requires less 
// arguments and prints the results. The other returns the results'
// values as arguments.
//
// calcAllErrors(h3readDir + projFileName) calculates and prints the
// errors on rB, gamma, delta for a particular file.
//
// calcQSubAllErrors() calculates the all errors for all the relevant
// files and does the appropriate quadrature subtractions to give the
// statistical, systematic, and model errors separately.


#include <iomanip>
#include "../BToDKTo3piK/confInt/makeProjs.cc"  // for readProjHists()

void calcErrors(const TH1F * hist, double cl1sigma = 0.199) {
  double val, errLo, errHi;
  int  valBin, errLoBin, errHiBin;
  double minCL;

  calcErrors(hist, val, errLo, errHi, valBin, errLoBin, errHiBin,
	     minCL, cl1sigma);

  cout << "Lowest CL is at " << val << " at bin " << valBin
       << ", CL=" << minCL << endl
       << "Upper error is " << errHi << " at bin " << errHiBin 
       << ", at value " << errHi + val << endl
       << "Lower error is " << errLo << " at bin " << errLoBin 
       << ", at value " << val - errLo << endl;
}



void calcErrors(const TH1F * hist, 
		double & val, double & errLo, double & errHi,
		int & valBin, int & errLoBin, int & errHiBin, 
		double & minCL, double cl1sigma = 0.199) {

  val = -999;
  errLo = -999;
  errHi = -999;
  valBin = -999;
  errLoBin = -999;
  errHiBin = -999;
  minCL = 100;

  int b;

  // First, find the minimum point:
  for (b = 1; b <= hist->GetNbinsX(); ++b) {
    if (minCL >= hist->GetBinContent(b)) {
      minCL = hist->GetBinContent(b);
      valBin = b;
      val = hist->GetBinCenter(b);
    }
  }

  // Now start at valBin and move from there to find the upper error
  for (b = valBin; b <= hist->GetNbinsX(); ++b) {
    if (hist->GetBinContent(b) >= cl1sigma - minCL) {
      errHiBin = b;
      errHi = hist->GetBinCenter(b) - val;
      break;  // stop at the first occurrence 
    }
  }

  // Now start at valBin and move from there to find the lower error
  for (b = valBin; b >= 1; --b) {
    if (hist->GetBinContent(b) >= cl1sigma - minCL) {
      errLoBin = b;
      errLo = val - hist->GetBinCenter(b);
      break;  // stop at the first occurrence 
    }
  }
}




void calcAllErrors(const char * fileName, // set to  h3readDir + projFileName,
		   const char * name = h3HistName,
		   double cl1sigma = 0.199) {

  readProjHists(fileName, name);
  
  cout << endl << "--rB--" << endl;
  calcErrors(hR, cl1sigma);

  cout << endl << "--gamma--" << endl;
  calcErrors(hG, cl1sigma);

  cout << endl << "--delta--" << endl;
  calcErrors(hD, cl1sigma);
}



calcQSubAllErrors(double cl1sigma = 0.199) {
  int  valBin, errLoBin, errHiBin;
  double minCL;

  const int NFILES = 3;  // for stat, stat+syst, all
  const int S(0), S2(1), S2M(2); // stat, stat+syst, stat+syst+model (=all)

  double valR[NFILES], valG[NFILES], valD[NFILES];
  double errLoR[NFILES], errLoG[NFILES], errLoD[NFILES];
  double errHiR[NFILES], errHiG[NFILES], errHiD[NFILES];

  // all errors:
  readProjHists(h3readDir + projFileName, h3HistName);

  calcErrors(hR, valR[S2M], errLoR[S2M], errHiR[S2M],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hG, valG[S2M], errLoG[S2M], errHiG[S2M],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hD, valD[S2M], errLoD[S2M], errHiD[S2M],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);


  // stat+syst errors:
  readProjHists(h3readDir + projFileNameStatSyst, h3HistNameStatSyst);

  calcErrors(hR, valR[S2], errLoR[S2], errHiR[S2],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hG, valG[S2], errLoG[S2], errHiG[S2],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hD, valD[S2], errLoD[S2], errHiD[S2],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);


  // stat errors:
  readProjHists(h3readDir + projFileNameStat, h3HistNameStat);

  calcErrors(hR, valR[S], errLoR[S], errHiR[S],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hG, valG[S], errLoG[S], errHiG[S],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);

  calcErrors(hD, valD[S], errLoD[S], errHiD[S],
	     valBin, errLoBin, errHiBin, minCL, cl1sigma);


  // Now calculate from this just the syst and the model errors:
  double errLoRSyst = sqrt(pow(errLoR[S2], 2) - pow(errLoR[S], 2));
  double errLoGSyst = sqrt(pow(errLoG[S2], 2) - pow(errLoG[S], 2));
  double errLoDSyst = sqrt(pow(errLoD[S2], 2) - pow(errLoD[S], 2));

  double errHiRSyst = sqrt(pow(errHiR[S2], 2) - pow(errHiR[S], 2));
  double errHiGSyst = sqrt(pow(errHiG[S2], 2) - pow(errHiG[S], 2));
  double errHiDSyst = sqrt(pow(errHiD[S2], 2) - pow(errHiD[S], 2));


  double errLoRModel = sqrt(pow(errLoR[S2M], 2) - pow(errLoR[S2], 2));
  double errLoGModel = sqrt(pow(errLoG[S2M], 2) - pow(errLoG[S2], 2));
  double errLoDModel = sqrt(pow(errLoD[S2M], 2) - pow(errLoD[S2], 2));

  double errHiRModel = sqrt(pow(errHiR[S2M], 2) - pow(errHiR[S2], 2));
  double errHiGModel = sqrt(pow(errHiG[S2M], 2) - pow(errHiG[S2], 2));
  double errHiDModel = sqrt(pow(errHiD[S2M], 2) - pow(errHiD[S2], 2));

  // Print:
  cout << setw(7) 
       << "Errors: Stat,   Syst,   Model:" << endl
       << "------------------------------" << endl
       << "rB+:    " << errHiR[S] << ", " << errHiRSyst << ", " << errHiRModel
       << endl 
       << "rB- :   " << errLoR[S] << ", " << errLoRSyst << ", " << errLoRModel
       << endl << endl
       << "gamma+: " << errHiG[S] << ", " << errHiGSyst << ", " << errHiGModel
       << endl 
       << "gamma-: " << errLoG[S] << ", " << errLoGSyst << ", " << errLoGModel
       << endl << endl
       << "delta+: " << errHiD[S] << ", " << errHiDSyst << ", " << errHiDModel
       << endl
       << "delta-: " << errLoD[S] << ", " << errLoDSyst << ", " << errLoDModel
       << endl;
}



