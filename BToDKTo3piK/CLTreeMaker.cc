
#include <math.h>
#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TTree.h"
#include "TH3F.h"
#include "TRandom.h"
#include "TString.h"

#include "BToDKTo3piK/CLTreeMaker.hh"

ClassImp(CLTreeMaker)

CLTreeMaker::CLTreeMaker(bool bigTree,
			 int resolution,
			 double maxNDof,
			 double dataRhoN,
			 double dataThetaN,
			 double dataRhoP,
			 double dataThetaP,
			 double dataRhoErrN,
			 double dataRhoErrP,
			 double dataThetaErrN,
			 double dataThetaErrP) :
  _bigTree(bigTree),

  // parms for chi^2 integral calculation resolution:
  _resolution(resolution),
  _maxNDof(maxNDof),
  _chi2Step_4(maxNDof * 4 / resolution),

  // The data results:
  _dataRhoN(dataRhoN),
  _dataThetaN(dataThetaN),
  _dataRhoP(dataRhoP),
  _dataThetaP(dataThetaP),

  // The data errors:
  _dataRhoErrN(dataRhoErrN),
  _dataRhoErrP(dataRhoErrP),
  _dataThetaErrN(dataThetaErrN),
  _dataThetaErrP(dataThetaErrP),

  // Variation of stat errs with true rho:
  _rhoErr(0.2345, -0.2818, 0.1899, -0.04834),
  _thetaErr(80.59, -141.4, 117.8, -35.52),

  // Variation of stat err spreads with true rho:
  _rhoErrSpread(0.08218, -0.1961, 0.1706, -0.05117),
  _thetaErrSpread(40.9, -98.06, 94.26, -31.66),

  // Variation of stat err pull sigmas with true rho:
  _rhoPullSigma(1.003, -0.256, 0.3339, -0.1355),
  _thetaPullSigma(1.188, -0.06302, 0.08119, -0.04653),

  // Variation of stat err pull biases with true rho:
  _biasPullRho(-0.02635, -0.0845, 0.05691, 0.009782),
  _biasPullTheta(0.006253, -0.04283, 0.07425, -0.05055),

  // Syst err calculators:
  _errMatSyst(kFALSE),          // without model
  _errMatModel(kFALSE, kTRUE),  // model only

  // Other:
  _verbose(kFALSE)
{
  // Stat err calculator:
  setupErrMatStatGauss();

  // Prepare chi^2 lookup table:
  makeTable();
}



CLTreeMaker::~CLTreeMaker() {}



void CLTreeMaker::setupErrMatStatGauss(Bool_t plotPoly3) {
  if (plotPoly3){
    // check that we did the Poly3 objects:right:
    TCanvas * can = new TCanvas("poly3", "poly3 objects", 800, 800);
    can->Divide(2,2);
    
    can->cd(1);
    _rhoErr.draw();
    
    can->cd(2);
    _thetaErr.draw();
    
    can->cd(3);
    _rhoErrSpread.draw();
    
    can->cd(4);
    _thetaErrSpread.draw();
  }

  // The number of spreads the errors are from the average in the data fit.
  // The rho results of the data fit are taken as true rho for the purpose
  // of calculating the spreads:
  double rhoErrNSpreads[2] = {
    (_dataRhoErrN - _rhoErr.getVal(_dataRhoN)) / _rhoErrSpread.getVal(_dataRhoN),
    (_dataRhoErrP - _rhoErr.getVal(_dataRhoP)) / _rhoErrSpread.getVal(_dataRhoP)
  };
  
  double thetaErrNSpreads[2] = {
    (_dataThetaErrN - _thetaErr.getVal(_dataRhoN)) / _thetaErrSpread.getVal(_dataRhoN),
    (_dataThetaErrP - _thetaErr.getVal(_dataRhoP)) / _thetaErrSpread.getVal(_dataRhoP)
  };

  // Print out how the theta NSpreads are calculated:
  //  cout << dataThetaErrN << " - " << thetaErr.getVal(dataRhoN)
  //       << " / " << thetaErrSpread.getVal(dataRhoN) 
  //       << " = " << thetaErrNSpreads[0] << endl;
  //  cout << dataThetaErrP << " - " << thetaErr.getVal(dataRhoP)
  //       << " / " << thetaErrSpread.getVal(dataRhoP) 
  //       << " = " << thetaErrNSpreads[1] << endl;
  

  if (plotPoly3){  
    cout << "Spreads are "  << endl
	 << "  rho-: " << rhoErrNSpreads[0] << endl
	 << "  theta-: " << thetaErrNSpreads[0] << endl
	 << "  rho+: " << rhoErrNSpreads[1] << endl
	 << "  theta+: " << thetaErrNSpreads[1] << endl;
  }

  _errMatStat.init(_rhoErr, _thetaErr, 
		   _rhoErrSpread, _thetaErrSpread,
		   _rhoPullSigma, _thetaPullSigma,
		   rhoErrNSpreads, thetaErrNSpreads);
}





// Make a confidence-level and chi^2 tree from the raw MC experiments:
void CLTreeMaker::makeCLTree(int nBinsR,
			     int nBinsGD,
			     const char * writeFilename,
			     double rBMin,
			     double rBMax,
			     double gammaMin,
			     double gammaMax,
			     double deltaMin,
			     double deltaMax) {

  TString fileName = writeFilename;
  fileName += ".root";
  TFile * file = 0;
  if (0 != writeFilename) {
    // Open the file to write the _tree:
    file = new TFile (fileName, "RECREATE");
  }

  // Create the tree. The branch vars are float to save space, but 
  // the intermediate vars are double for high precision. So copy
  // from double to float before filling ntuple:
  _tree = new TTree("makeCLTree", "makeCLTree");

  double rB;
  float rBF;
  _tree->Branch("rB", &rBF, "rB/F");

  double gammaVal;  // needs to be != "gamma", which appears in math.h
  float gammaValF;
  _tree->Branch("gamma", &gammaValF, "gamma/F");

  double delta;
  float deltaF;
  _tree->Branch("delta", &deltaF, "delta/F");

  double genRhoN;
  float genRhoNF;
  if (_bigTree) _tree->Branch("genRhoN", &genRhoNF, "genRhoN/F");

  double genThetaN;
  float genThetaNF;
  if (_bigTree) _tree->Branch("genThetaN", &genThetaNF, "genThetaN/F");

  double genRhoP;
  float genRhoPF;
  if (_bigTree) _tree->Branch("genRhoP", &genRhoPF, "genRhoP/F");

  double genThetaP;
  float genThetaPF;
  if (_bigTree) _tree->Branch("genThetaP", &genThetaPF, "genThetaP/F");

  double chi2Stat[2];
  float chi2StatF[2];
  if (_bigTree) {
    _tree->Branch("chi2Stat", &chi2StatF[0], "chi2Stat/F");
    _tree->Branch("chi2StatBias", &chi2StatF[1], "chi2StatBias/F");
  }

  double clStat[2];
  float clStatF[2];
  _tree->Branch("clStat", &clStatF[0], "clStat/F");
  _tree->Branch("clStatBias", &clStatF[1], "clStatBias/F");

  double chi2StatSyst[2];
  float chi2StatSystF[2];
  if (_bigTree) {
    _tree->Branch("chi2StatSyst", &chi2StatSystF[0], "chi2StatSyst/F");
    _tree->Branch("chi2StatSystBias", &chi2StatSystF[1], "chi2StatSystBias/F");
  }

  double clStatSyst[2];
  float clStatSystF[2];
  _tree->Branch("clStatSyst", &clStatSystF[0], "clStatSyst/F");
  _tree->Branch("clStatSystBias", &clStatSystF[1], "clStatSystBias/F");

  double chi2All[2];
  float chi2AllF[2];
  if (_bigTree) {  
    _tree->Branch("chi2All", &chi2AllF[0], "chi2All/F");
    _tree->Branch("chi2AllBias", &chi2AllF[1], "chi2AllBias/F");
  }

  double clAll[2];
  float clAllF[2];
  _tree->Branch("clAll", &clAllF[0], "clAll/F");
  _tree->Branch("clAllBias", &clAllF[1], "clAllBias/F");

  // the rB step is the same as the gamma and delta step, calculated
  // with the same nBinsGD. This isbecause each job does a different 
  // section of the rB range, determined with rBMin:
  const double rBStep = (rBMax - rBMin) / nBinsGD;
  const double gammaStep = (gammaMax - gammaMin) / nBinsGD;
  const double deltaStep = (deltaMax - deltaMin) / nBinsGD;

  // loop over points:
  for (int iR = 0; iR < nBinsR; ++iR){
    rB = rBMin + (iR + 0.5) * rBStep;
    for (int iG = 0; iG < nBinsGD; ++iG){
      gammaVal = gammaMin + (iG + 0.5) * gammaStep;
      for (int iD = 0; iD < nBinsGD; ++iD){
	delta = deltaMin + (iD + 0.5) * deltaStep;

	// Calculate the experimental coordinates:
	coordTrans(rB, gammaVal, delta, &genRhoN, &genThetaN, &genRhoP, &genThetaP);
	
	
	
	if (_verbose) {
	  cout << "rB=" << rB << " gamma=" << gammaVal << " delta=" << delta 
	       << endl
	       << "rM=" << genRhoN << " tM=" << genThetaN 
	       << " rP=" << genRhoP << " tP=" << genThetaP << endl;
	}
	
	// The stat error matrix:
	TMatrixD errMatStat = 
	  _errMatStat.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// The syst error matrix:
	TMatrixD errMatSyst = 
	  _errMatSyst.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// The model error matrix:
	TMatrixD errMatModel = 
	  _errMatModel.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// From the true experimental coordinates, calculate the difference
	// vectors to the point we found in the data fit:
	const int NEXPCOORD = 4;
	
	TMatrixD diff(NEXPCOORD, 1);
	diff(0,0) = _dataRhoN - genRhoN; 
	diff(1,0) = _dataThetaN - genThetaN;
	diff(2,0) = _dataRhoP - genRhoP;
	diff(3,0) = _dataThetaP - genThetaP;
	
	// Add the bias, which is a product of the measured error for 
	// the appropriate variable times the pull bias on that variable,
	// given the true value of rho:
	for (int b = 0; b < 2; ++b) {
	  diff(0,0) -= b 
	    * _errMatStat.measuredErrRho(genRhoN, ErrMatStatGauss::MINUS)
	    * _biasPullRho.getVal(genRhoN);
	  
	  diff(1,0) -= b 
	    * _errMatStat.measuredErrTheta(genRhoN, ErrMatStatGauss::MINUS)
	    * _biasPullTheta.getVal(genRhoN);
	  
	  diff(2,0) -= b 
	    * _errMatStat.measuredErrRho(genRhoP, ErrMatStatGauss::PLUS)
	    * _biasPullRho.getVal(genRhoP);
	  
	  diff(3,0) -= b 
	    * _errMatStat.measuredErrTheta(genRhoP, ErrMatStatGauss::PLUS)
	    * _biasPullTheta.getVal(genRhoP);
	  
	  // Take the transpose of the diff:
	  TMatrixD diffT(1, NEXPCOORD);
	  diffT.Transpose(diff);
	  
	  // Calculate the statistical chi^2 and 4D confidence level:
	  const int NDOF = 4;
	  calcChi2Cl(diff, diffT, errMatStat, NDOF, chi2Stat[b], clStat[b]);
	  
	  // The stat+syst chi2 and conf level:
	  TMatrixD errMatStatSyst = errMatStat + errMatSyst;
	  calcChi2Cl(diff, diffT, errMatStatSyst, NDOF, chi2StatSyst[b], clStatSyst[b]);
	  
	  // The stat+syst+model chi2 and conf level:
	  TMatrixD errMatAll = errMatStatSyst + errMatModel;
	  calcChi2Cl(diff, diffT, errMatAll, NDOF, chi2All[b], clAll[b]);
	  
	  if (_verbose) {
	    cout << "Bias factor = " << b << endl;
	    diff.Print();
	    diffT.Print();
	    errMatStat.Print();
	    errMatAll.Print();
	    
	    cout << "chi2All=" << chi2All[b] << ", clAll=" << clAll[b];
	    cout << ", chi2Stat=" << chi2Stat[b] << ", clStat=" << clStat[b];
	    cout << ", chi2StatSyst=" << chi2StatSyst[b] << ", clStatSyst=" 
		 << clStatSyst[b];
	    cout << endl << endl;
	  }
	}
	
	// Fill the _tree and histograms:
	
	rBF = rB;
	gammaValF = gammaVal;   
	deltaF = delta;
	genRhoNF = genRhoN;
	genThetaNF = genThetaN;
	genRhoPF = genRhoP;
	genThetaPF = genThetaP;
	for (int i = 0; i < 2; ++i) {
	  chi2StatF[i] = chi2Stat[i];
	  clStatF[i] = clStat[i];
	  chi2StatSystF[i] = chi2StatSyst[i];
	  clStatSystF[i] = clStatSyst[i];
	  chi2AllF[i] = chi2All[i];
	  clAllF[i] = clAll[i];
	}
	
	_tree->Fill();

      }
    }
  } 

  if (0 != writeFilename) {
    // Write the _tree:
    _tree->Write();
    file->Close();
  }
}


// Make a confidence-level and chi^2 TH3F from the raw MC experiments:
void CLTreeMaker::makeCLHist3(int nBinsR,
			      int nBinsGD,
			      const char * writeFilename,
			      double rBMin,
			      double rBMax,
			      double rBFullRange,
			      double gammaMin,
			      double gammaMax,
			      double deltaMin,
			      double deltaMax) {
  
  TString fileName = writeFilename;
  fileName += ".root";
  TFile * file = 0;
  if (0 != writeFilename) {
    // Open the file to write the _tree:
    file = new TFile (fileName, "RECREATE");
  }

  // Create the TH3F's: 
  _hCLStat = new TH3F("hCLStat", "hCLStat",
		    nBinsR, rBMin, rBMax,
		    nBinsGD, gammaMin, gammaMax,
		    nBinsGD, deltaMin, deltaMax);

  _hCLStatSyst = new TH3F("hCLStatSyst", "hCLStatSyst",
		    nBinsR, rBMin, rBMax,
		    nBinsGD, gammaMin, gammaMax,
		    nBinsGD, deltaMin, deltaMax);

  _hCLAll = new TH3F("hCLAll", "hCLAll",
		    nBinsR, rBMin, rBMax,
		    nBinsGD, gammaMin, gammaMax,
		    nBinsGD, deltaMin, deltaMax);

  _hChi2All = new TH3F("hChi2All", "hChi2All",
		    nBinsR, rBMin, rBMax,
		    nBinsGD, gammaMin, gammaMax,
		    nBinsGD, deltaMin, deltaMax);

  _hCLAllBias = new TH3F("hCLAllBias", "hCLAllBias",
		    nBinsR, rBMin, rBMax,
		    nBinsGD, gammaMin, gammaMax,
		    nBinsGD, deltaMin, deltaMax);


  // the rB step is the same as the gamma and delta step, calculated
  // with the same nBinsGD. This is because each job does a different 
  // section of the rB range, determined with rBMin:
  const double rBStep = rBFullRange / nBinsGD;
  const double gammaStep = (gammaMax - gammaMin) / nBinsGD;
  const double deltaStep = (deltaMax - deltaMin) / nBinsGD;

  // loop over points:
  for (int iR = 1; iR <= nBinsR; ++iR){
    double rB = rBMin + (iR - 0.5) * rBStep;
    for (int iG = 1; iG <= nBinsGD; ++iG){
      double gammaVal = gammaMin + (iG - 0.5) * gammaStep;
      for (int iD = 1; iD <= nBinsGD; ++iD){
	double delta = deltaMin + (iD - 0.5) * deltaStep;

	// Calculate the experimental coordinates:
	double genRhoN;
	double genThetaN;
	double genRhoP;
	double genThetaP;

	coordTrans(rB, gammaVal, delta, 
		   &genRhoN, &genThetaN, &genRhoP, &genThetaP);
	
	if (_verbose) {
	  cout << "rB=" << rB << " gamma=" << gammaVal << " delta=" << delta 
	       << endl
	       << "rM=" << genRhoN << " tM=" << genThetaN 
	       << " rP=" << genRhoP << " tP=" << genThetaP << endl;
	}
	
	// The stat error matrix:
	TMatrixD errMatStat = 
	  _errMatStat.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// The syst error matrix:
	TMatrixD errMatSyst = 
	  _errMatSyst.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// The model error matrix:
	TMatrixD errMatModel = 
	  _errMatModel.matrix(genRhoN, genThetaN, genRhoP, genThetaP);
	
	// From the true experimental coordinates, calculate the difference
	// vectors to the point we found in the data fit:
	const int NEXPCOORD = 4;
	
	TMatrixD diff(NEXPCOORD, 1);
	diff(0,0) = _dataRhoN - genRhoN; 
	diff(1,0) = _dataThetaN - genThetaN;
	diff(2,0) = _dataRhoP - genRhoP;
	diff(3,0) = _dataThetaP - genThetaP;
	
	double chi2Stat[2];
	double clStat[2];
	double chi2StatSyst[2];
	double clStatSyst[2];
	double chi2All[2];
	double clAll[2];

	// Add the bias, which is a product of the measured error for 
	// the appropriate variable times the pull bias on that variable,
	// given the true value of rho:
	for (int b = 0; b < 2; ++b) {
	  diff(0,0) -= b 
	    * _errMatStat.measuredErrRho(genRhoN, ErrMatStatGauss::MINUS)
	    * _biasPullRho.getVal(genRhoN);
	  
	  diff(1,0) -= b 
	    * _errMatStat.measuredErrTheta(genRhoN, ErrMatStatGauss::MINUS)
	    * _biasPullTheta.getVal(genRhoN);
	  
	  diff(2,0) -= b 
	    * _errMatStat.measuredErrRho(genRhoP, ErrMatStatGauss::PLUS)
	    * _biasPullRho.getVal(genRhoP);
	  
	  diff(3,0) -= b 
	    * _errMatStat.measuredErrTheta(genRhoP, ErrMatStatGauss::PLUS)
	    * _biasPullTheta.getVal(genRhoP);
	  
	  // Take the transpose of the diff:
	  TMatrixD diffT(1, NEXPCOORD);
	  diffT.Transpose(diff);
	  
	  // Calculate the statistical chi^2 and 4D confidence level:
	  const int NDOF = 4;
	  calcChi2Cl(diff, diffT, errMatStat, NDOF, chi2Stat[b], clStat[b]);
	  
	  //if (1==iG && 1==iD && 0==b) {
	  //  cout << "------------------------" << endl;
	  //  cout << "iR="<< iR << " rB=" << rB << endl;
	  //  cout << "iG="<< iG << " gamma=" << gammaVal << endl;
	  //  cout << "iD="<< iD << " delta=" << delta << endl;
	  //  cout << "diff= ";
	  //  diff.Print();
	  //  cout << "clStat=" << clStat[0] << endl;
	  //}

	  // The stat+syst chi2 and conf level:
	  TMatrixD errMatStatSyst = errMatStat + errMatSyst;
	  calcChi2Cl(diff, diffT, errMatStatSyst, NDOF, chi2StatSyst[b], clStatSyst[b]);
	  
	  // The stat+syst+model chi2 and conf level:
	  TMatrixD errMatAll = errMatStatSyst + errMatModel;
	  calcChi2Cl(diff, diffT, errMatAll, NDOF, chi2All[b], clAll[b]);
	  
	  if (_verbose) {
	    cout << "Bias factor = " << b << endl;
	    diff.Print();
	    diffT.Print();
	    errMatStat.Print();
	    errMatAll.Print();
	    
	    cout << "chi2All=" << chi2All[b] << ", clAll=" << clAll[b];
	    cout << ", chi2Stat=" << chi2Stat[b] << ", clStat=" << clStat[b];
	    cout << ", chi2StatSyst=" << chi2StatSyst[b] << ", clStatSyst=" 
		 << clStatSyst[b];
	    cout << endl << endl;
	  }
	}
	
	// Fill the histograms:
	_hCLStat->SetBinContent(iR, iG, iD, clStat[0]);
	_hCLStatSyst->SetBinContent(iR, iG, iD, clStatSyst[0]);
	_hCLAll->SetBinContent(iR, iG, iD, clAll[0]);
	_hChi2All->SetBinContent(iR, iG, iD, chi2All[0]);
	_hCLAllBias->SetBinContent(iR, iG, iD, clAll[1]);
      }
    }
  } 

  if (0 != writeFilename) {
    // Write the _tree:
    _hCLStat->Write();
    _hCLStatSyst->Write();
    _hCLAll->Write();
    _hChi2All->Write();
    _hCLAllBias->Write();

    file->Close();
  }
}


void 
CLTreeMaker::coordTrans(double rB, double gamma, double delta, 
			double * rM, double * tM, double * rP, double * tP,
			double * xM, double * yM, double * xP, double * yP) 
{
  double g = gamma * pi() / 180;
  double d = delta * pi() / 180;
  
  double xMinus = rB * cos (d - g);
  double yMinus = rB * sin (d - g);
  double xPlus  = rB * cos (d + g);
  double yPlus =  rB * sin (d + g);

  if (0 != xM) {*xM = xMinus;}
  if (0 != yM) {*yM = yMinus;}
  if (0 != xP) {*xP = xPlus;}
  if (0 != yP) {*yP = yPlus;}
  
  const double x0 = 0.850;
  double rMinus = sqrt((xMinus - x0) * (xMinus - x0) + yMinus * yMinus);
  double rPlus  = sqrt((xPlus - x0)  * (xPlus - x0)  + yPlus  * yPlus);
  double tMinus = atan2(yMinus, (xMinus - x0)) * 180./pi();
  double tPlus  = atan2(yPlus,  (xPlus - x0))  * 180./pi();

  if (tPlus < 0)  tPlus  += 360;
  if (tMinus < 0) tMinus += 360;

  if (0 != rM) {*rM = rMinus;}
  if (0 != tM) {*tM = tMinus;}
  if (0 != rP) {*rP = rPlus;}
  if (0 != tP) {*tP = tPlus;}
}




// calculates Gamma(ndof/2), needed to normalize the chi2 distribution.
// To get Gamma(x), call gammaHalf(x*2):
double CLTreeMaker::gammaHalf(int ndof) {
  double result = 1;  
  int nhalf = ndof / 2;
  if (nhalf * 2 == ndof) {
    // ndof is even, so make use of the fact that the Gamma function for
    // an integer n is (n-1)!.
    for (int n = 2; n <= nhalf - 1; ++n) {
      result *= (double)n;
    }
  }
  else {
    // Otherwise, use the facts that Gamma(x+1) = x * Gamma(x) and
    // Gamma(1/2) = sqrt(pi):
    static const double SQRTPI = sqrt(pi());
    result = SQRTPI;

    for (double x = 0.5; x < ((double)ndof)/2; ++x) {
      result *= x;
    }
  }
  return result;
}


// The chi^2 distribution value for chi2 and ndof:
double CLTreeMaker::chi2dist(double chi2, int ndof) {
  double ndofHalf = ((double)ndof)/2;
  double denom = pow((double)2.0, ndofHalf) * gammaHalf(ndof);
  double num = pow(chi2, ndofHalf - 1) * exp(-chi2/2);
  return num / denom;
}


// The integral of the chi2dist for ndof DOF, from chi2 to
// chi2 + ndof * maxNdof (taken to estimate infinity), calculated in
// steps of size ndof * stepNfod:
double CLTreeMaker::pValue(double chi2,
			   int ndof,
			   double stepNdof,
			   double maxNdof) {
  double result = 0;
  for (double x = chi2; x < chi2 + maxNdof * ndof; x += stepNdof * ndof) {
    result += chi2dist(x, ndof);
  }
  return result * stepNdof * ndof;
}


// cl() returns 1-pValue(), recalculated for better accuracy:
double CLTreeMaker::cl(double chi2,
		       int ndof,
		       double stepNdof,
		       double maxNdof) {
  // For 1 DOF, avoid the divergence at x=0:
  switch(ndof) {
  case 1: 
    return 1.0 - pValue(chi2, ndof, stepNdof, maxNdof);
    break;

  default:    
    double result = 0;
    for (double x = 0; x < chi2; x += stepNdof * ndof) {
      result += chi2dist(x, ndof);
    }
    return result * stepNdof * ndof;
  }
}


// Calculate the chi^2 and CL given 2 vectors and an error matrix: 
void CLTreeMaker::calcChi2Cl(const TMatrixD & diff, 
			     const TMatrixD & diffT, 
			     const TMatrixD & errMat,
			     int nDOF,
			     double & chi2, 
			     double & clevel) const {
  TMatrixD invErrMat = errMat;
  invErrMat.Invert();
  TMatrixD chi2Mat = diffT * invErrMat * diff;
  chi2 = chi2Mat(0,0);

  if (4 == nDOF) { // then we have a lookup table:
    clevel = lookupCL_4(chi2);
  }
  else {
    // For efficiency, only calculate if the chi^2/nDOF is different from
    // unity by less than _maxNDof widths of the chi^2 distribution, 
    // a width being sqrt(2*nDOF). Otherwise, return 1.1:    
    if ((chi2/nDOF - 1) / sqrt(2.0 * nDOF) < _maxNDof) {
      clevel = cl(chi2, nDOF); 
    }
    else {
      clevel = 1.1;
    }
  }
}


void CLTreeMaker::makeTable() {
  // Make the chi^2 lookup table:
  _chi2Int_4 = new double[_resolution];

  _chi2Int_4[0] = 0;  // first entry is 0
  for (int i = 1; i < _resolution; ++i){  // integrate from next entry
    double chi2_4 = i * _chi2Step_4;          // chi^2 for this entry
    double dist = chi2dist(chi2_4, 4);          // the chi^2 distribution:
    _chi2Int_4[i] = dist * _chi2Step_4 + _chi2Int_4[i-1];   // the integral 
  }
}


double CLTreeMaker::lookupCL_4(double chi2) const {
  // Find the entry in the nDOF=4 table corresponding to this chi^2:
  int i = (int)(chi2 / _chi2Step_4);
  if (i < _resolution) {
    return _chi2Int_4[i];
  }
  return 1.0;  // chi2 was so big that the integral is 1
}

  
