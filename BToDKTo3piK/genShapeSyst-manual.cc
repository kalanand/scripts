// Generates the csh files needed to submit jobs to propagate errors
// associated with shapes.
// To compile, need to add genShapeSyst to EXTRABINS in GNUmakefile,
// and add CLHEP to link_BToDKTo3piK

#include <iostream>
#include <string>
using std::string;
using std::cout;
using std::endl;
using std::cerr;


#include <iomanip.h>
#include <fstream.h>
#include <strstream.h>
#include <assert.h>
#include <math.h>
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"


void check(const char slash, const char trueSlash) {
  if (trueSlash != slash){
    cerr << "read \'" << slash << "\' which is not \'" << trueSlash
	 << "\'. Did you give the wrong # of lines?" << endl;
    assert(0);
  }
}


void readMatrix(HepSymMatrix & matrix, int size, istream & stream = cin){
  for (int i = 0; i < size; ++i){
    // remove the first '//':
    char slash;
    stream.get(slash);
    check(slash, '/');
    stream.get(slash);
    check(slash, '/');
    
    // read rest of line:
    const int SIZE = 1000;
    char line[SIZE];
    stream.getline(line, SIZE);
    istrstream lineStream(line);

    // then read the #'s:
    for (int j = 0; j < size; ++j){
      double d;
      lineStream >> d;
      matrix(i+1, j+1) = d;
    }
  }
}
      
// -----------------------------------------------------------------------
int genShapeSyst() {
  cout << "enter base name for output file " << ends;
  string baseName;
  cin >> baseName;

  cout << "enter # of parameters " << ends;
  int npar;
  cin >> npar;

  cout << "enter " << npar << " parameter lines" << endl;
  string * names = new string[npar];
  HepMatrix tempParams(npar, 1);
  HepMatrix tempErrors(npar, 1);
  const int SIZE = 1000;
  char line[SIZE];
  int nFreePar = 0;

  for (int p = 0; p < npar; ++p) {
    // read full line to a temp char* that is then discarded:
    cin.getline(line, SIZE, '=');

    // read name and value from line:
    istrstream lineStream(line);
    char name[SIZE];
    lineStream.getline(name, SIZE, '=');

    double param;
    cin >> param;

    char plusMinus[20];
    cin >> plusMinus;

    double error;
    cin >> error;

    // figure out if the param is fixed:
    cin.getline(line, SIZE, '(');
    if (0 == strchr(line, 'C')) {
      // not fixed, take it and increment index of free pars::
      tempParams(nFreePar+1, 1) = param;
      tempErrors(nFreePar+1, 1) = error;
      names[nFreePar] = name;    
      ++nFreePar;
    }
    // read to the end of the cin line:
    cin.getline(line, SIZE);    
  }

  // Copy tempParams into a final vector of the right size:
  HepMatrix params(nFreePar, 1);  
  HepMatrix errors(nFreePar, 1);  
  int fp;
  for (fp = 0; fp < nFreePar; ++fp) {
    params(fp+1, 1) = tempParams(fp+1, 1);
    errors(fp+1, 1) = tempErrors(fp+1, 1);
  }

  cout << "Found " << nFreePar << " floating parameters." << endl;
  cout << "enter " << nFreePar << "x" << nFreePar 
       << " error matrix with leading //" 
       << endl;
  
  HepSymMatrix mat(nFreePar);
  readMatrix(mat, nFreePar);
  HepSymMatrix diag = mat;
  HepMatrix trans = diagonalize(&diag);

  cout << "I got " << endl;
  for (fp = 0; fp < nFreePar; ++fp) {
    cout << names[fp] << " = " << params(fp+1, 1) << " +/- " 
	 << errors(fp+1, 1) << endl;
  }
  cout << "matrix= " << mat << endl
       << "diag  = " << diag << endl
       << "trans  = " << trans << endl;

  // Check errors and matrix:
  for (fp = 0; fp < nFreePar; ++fp) {
    const double TOL = 0.05;
    double diff = 1 - errors(fp+1, 1) / sqrt(mat(fp+1, fp+1));
    if (fabs(diff) > TOL) {
      cerr << "found relative diff of " << diff 
	   << " between error (" << errors(fp+1, 1) << ") and sqrt of "
	   << "  matrix element (" << sqrt(mat(fp+1, fp+1)) 
	   << " )for parameter " << fp << ": " 
	   << names[fp] << endl;
      assert(0);
    }
  }

  // for each +/- 1 sigma change in diag, generate a new par vector:
  for (fp = 0; fp < nFreePar; ++fp) {
    // get fp into a string (they should have had an operator+(int) for this!):
    ostrstream fpStream;
    fpStream << fp << ends;
    string fpString = fpStream.str();

    // shift param fp by + 1 sigma & write new par file:
    string parFilePlusName = baseName + "+";
    parFilePlusName += fpString;
    parFilePlusName += ".par";

    ofstream parFilePlus(parFilePlusName.c_str());

    parFilePlus << "setVar(" << "paramsOnResDK, "
		<< names[fp] << ", " 
		<< params(fp+1, 1) + sqrt(diag(fp, fp)) << endl;
    
    // shift by -1 sigma:
    string parFileMinusName = baseName + "-";
    parFileMinusName += fpString;
    parFileMinusName += ".par";

    ofstream parFileMinus(parFileMinusName.c_str());

    parFileMinus << "setVar(" << "paramsOnResDK, "
		 << names[fp] << ", " 
		 << params(fp+1, 1) - sqrt(diag(fp, fp)) << endl;    

    cout << "Wrote files " << parFilePlusName << " and " << parFileMinusName 
	 << endl;
  }
  
  return 0;
}


int main() {
  genShapeSyst();
}
