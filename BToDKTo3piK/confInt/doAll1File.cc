// Combines one file, projects it, plots it, and calculates its limits


#include "../BToDKTo3piK/utils/strDouble.cc"
#include "../BToDKTo3piK/confInt/calcErrors.cc"


void doAll1File(double trueRhoN, double trueThetaN,
		double trueRhoP, double trueThetaP,
		const TString & dir = h3writeDir, 
		double minRB = 0, double maxRB = 0.8) {
  
  TString fileName = strDouble(trueRhoN);
  fileName += "_";
  fileName += strDouble(trueThetaN);
  fileName += "_";
  fileName += strDouble(trueRhoP);
  fileName += "_";
  fileName += strDouble(trueThetaP);

  TString rawFile      = dir + fileName;
  TString rawFileRoot  = dir + fileName + ".root";
  TString sumFileRoot  = dir + fileName + "-sum.root";
  TString projFileRoot = dir + fileName + "-projections.root";

  
  cout << "Working on R-=" 
       << trueRhoN << " T-=" << trueThetaN << " R+="
       << trueRhoP << " T+=" << trueThetaP << endl
       << "Raw file name = " << rawFile << endl;

  CLTreeMaker maker(false, 1000000, 200, 
		    trueRhoN, trueThetaN, trueRhoP, trueThetaP);

  maker.makeCLHist3(100, 100, rawFile, minRB, maxRB, maxRB);


  const char * fileStr[1] = {(const char *)rawFileRoot};

  combineFiles(sumFileRoot, "hCLAll", fileStr, 1);

  makeProjs(sumFileRoot, projFileRoot);

  readProjHists(projFileRoot);

  plotProjs(0, rawFile);


  const double cl1sigma = 0.199;

  cout << endl;
  cout << "************ Limits for file ************" << fileName << endl;
  cout << "--rB--" << endl;
  calcErrors(hR, cl1sigma);

  cout << endl << "--gamma--" << endl;
  calcErrors(hG, cl1sigma);

  cout << endl << "--delta--" << endl;
  calcErrors(hD, cl1sigma);
  cout << "************ END Limits for file ************" << fileName << endl
       << endl;

}

