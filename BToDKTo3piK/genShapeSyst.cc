
// Generates the csh files needed to submit jobs to propagate errors
// associated with shapes.
// To compile, need to add genShapeSyst to EXTRABINS in GNUmakefile,
// and add CLHEP to link_BToDKTo3piK

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;


#include <stdlib.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream.h>
#include <assert.h>
#include <string>
#include <math.h>
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"


//----------------------------------------------------
void check(const char slash, const char trueSlash, int occurrence, int line) {
  if (trueSlash != slash){
    cerr << " on occurrence " << occurrence << " line " << line
	 << " read \'" << slash << "\' which is not \'" << trueSlash
	 << "\'. Did you give the wrong # of lines?" << endl;
    assert(0);
  }
}


//----------------------------------------------------
void readMatrix(HepSymMatrix & matrix, int size, istream & stream){
  for (int i = 0; i < size; ++i){
    // remove the first '//':
    char slash;
    stream.get(slash);
    check(slash, '/', 1, i);
    stream.get(slash);
    check(slash, '/', 2, i);
    
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
int genShapeSyst(istream & stream, 
		 string baseName, 
		 int npar, double scale = 1) {

  cout << "Reading " << npar << " parameter lines" << endl;
  string * names = new string[npar];
  HepMatrix tempParams(npar, 1);
  HepMatrix tempErrors(npar, 1);
  const int SIZE = 1000;
  char line[SIZE];
  int nFreePar = 0;

  for (int p = 0; p < npar; ++p) {
    // read full line to a temp char* that is then discarded:
    stream.getline(line, SIZE, '=');

    // read name and value from line:
    istrstream lineStream(line);
    char name[SIZE];
    lineStream.getline(name, SIZE, '=');

    double param;
    double error;
    stream >> param;

    // figure out if the param is fixed:
    stream.getline(line, SIZE, '(');
    if (0 == strchr(line, 'C')) {
      // not fixed, take it and increment index of free pars::
      tempParams(nFreePar+1, 1) = param;
      names[nFreePar] = name;    
      ++nFreePar;

      // read the error:
      istrstream lineStream(line);
      char plusMinus[20];
      lineStream >> plusMinus;      
//      lineStream >>  tempErrors(nFreePar+1, 1);
      lineStream >>  tempErrors(nFreePar+0, 1);
    }

    // read to the end of the stream line:
    stream.getline(line, SIZE);    
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
  cout << "Reading " << nFreePar << "x" << nFreePar 
       << " error matrix with leading //" 
       << endl;
  
  HepSymMatrix mat(nFreePar);
  readMatrix(mat, nFreePar, stream);
  HepSymMatrix diag = mat;
  HepMatrix trans = diagonalize(&diag);
   
  // the collelation matrix:
  HepSymMatrix corr(nFreePar);
  for (int row = 1; row <= nFreePar; ++row){
    for (int col = 1; col <= nFreePar; ++col) {
      corr(row, col) = mat(row, col) / sqrt(mat(row, row) * mat(col, col));
    }
  }

  cout << "I got " << endl;
  for (fp = 0; fp < nFreePar; ++fp) {
    cout << names[fp] << " = " << params(fp+1, 1) << " +/- " 
	 << errors(fp+1, 1) << endl;
  }
  cout << "matrix= " << mat << endl
       << "corr  = " << corr << endl
       << "diag  = " << diag << endl
       << "trans  = " << trans << endl;

  HepMatrix invTrans = trans;
  int fail;
  invTrans.invert(fail);
  if (0 != fail) {
    cerr << "*********************************************" << endl
	 << "Error inverting transformation matrix. Doing nothing." << endl
	 << "*********************************************" << endl;
  }

  // the parameters in the diagonal basis:
  HepMatrix paramsDiag = invTrans * params;

  cout << "recalculating diag: " << invTrans * mat * trans << endl;

  cout << "invTrans= " << invTrans << endl
       << "paramsDiag= " << paramsDiag << endl
       << "trans*paramsDiag = "<< trans*paramsDiag <<endl;
  
  // Check errors and matrix:
  for (fp = 0; fp < nFreePar; ++fp) {
    const double TOL = 0.05;
    double diff = 1 - errors(fp+1, 1) / sqrt(mat(fp+1, fp+1));
    if (fabs(diff) > TOL) {
      cerr << "WARNING: found relative diff of " << diff 
	   << " between error (" << errors(fp+1, 1) << ") and sqrt of "
	   << "  matrix element (" << sqrt(mat(fp+1, fp+1)) 
	   << " )for parameter " << fp << ": " 
	   << names[fp] << endl;
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


    cout<<" diag matrix elements "<<fp+1<<" - "<<fp+1<<" value "
        <<diag(fp+1,fp+1)<<endl;

    ofstream parFilePlus(parFilePlusName.c_str());

    HepMatrix paramsDiagShiftedP = paramsDiag;
    paramsDiagShiftedP(fp + 1, 1) += sqrt(diag(fp+1, fp+1)) * scale;
    HepMatrix paramsShiftedP = trans * paramsDiagShiftedP;

    cout<<" invTrans*paramsDiagP "<< paramsShiftedP <<endl;

    // fix this:
    for (int fpp = 0; fpp < nFreePar; ++fpp) {
      parFilePlus << names[fpp] << "  =  "
		  << paramsShiftedP(fpp+1, 1) << endl;
    }

    // shift by -1 sigma:
    string parFileMinusName = baseName + "-";
    parFileMinusName += fpString;
    parFileMinusName += ".par";

    ofstream parFileMinus(parFileMinusName.c_str());

    HepMatrix paramsDiagShiftedM = paramsDiag;
    paramsDiagShiftedM(fp + 1, 1) -= sqrt(diag(fp+1, fp+1)) * scale;
    HepMatrix paramsShiftedM = trans * paramsDiagShiftedM;

    cout<<" invTrans*paramsDiagM "<< paramsShiftedM <<endl;


    // fix this:
    for (int fpm = 0; fpm < nFreePar; ++fpm) {
      parFileMinus << names[fpm] << "  =  " 
		   << paramsShiftedM(fpm+1, 1)  << endl;
    }

    cout << "Wrote files " << parFilePlusName << " and " << parFileMinusName 
	 << endl;
  }
  
  return 0;
}

//----------------------------------------------------
void cleanStream(istream & inStream, string & outString) {
  // Clean up inStream and put the results in outString for one PDF fit:
  int nDollars = 0;
  while(nDollars < 2) {
    // read 1 char at a time:
    char c1 = 'Z';
    inStream.get(c1);

    if (inStream.eof()) {
      cout << "eof reached" << endl;
      break;
    }

    if (inStream.bad()) {
      cout << "inStream BAD" << endl;
      break;
    }

    assert('$' != c1);
    
    if ('/' != c1) {
      // non comment, copy it:
      outString += c1;
    }
    else {
      // see if it's a "//":
      char c2;
      inStream.get(c2);
      
      assert('$' != c2);

      if ('/' != c2) {
	// not a "//", copy it:
	outString += c1;
	outString += c2;
      }
      else {
	// a "//". See if it's a "//$":
	char c3; 
	inStream.get(c3);
	
	if ('$' != c3) {
	  // not a "//$", copy it:
	  outString += c1;
	  outString += c2;
	  outString += c3;
	}
	else {
	  // it's a "//$". Count $'s and don't put in outString:
	  ++nDollars;

	  // remove the rest of the line from the inStream:
	  char c4[1000];
	  inStream.getline(c4, 1000);
	  cout << "on cleaning line read \"" << c4 << "\"" << endl;
	} // end of c3 check
      } // end of c2 check
    } // end of c1 check
  } // end main loop
}

//----------------------------------------------------
int count(const char * key, const string & inStr, string & pdfName) {
  // count occurrences of key in inStr and figure out pdfName:
  int result = 0;
  int pos = 0;
  string lastPdfName;

  while (1) {
    pos = inStr.find(key, pos);
    if (string::npos == pos) {
      break;
    }
    ++result; // count
    ++pos;     // start next search on next char
    
    // get the next word, which is the pdfName:
    int startPos = inStr.find('.');
    int endPos = inStr.find('.', startPos+1);

    pdfName = inStr.substr(startPos, endPos-startPos);

    // check that his didn't chang:
    if (0 != lastPdfName.length()) {
      assert(lastPdfName == pdfName);
    }

    // store for comparison later next loop:
    lastPdfName = pdfName;
  }
  return result;
}

//----------------------------------------------------
int main() {
  cout << "enter path to input files (    params/DK/    ): " << ends;
  string inPath;
  cin >> inPath;

  cout << "enter path for output files (    syst/    ): " << ends;
  string outPath;
  cin >> outPath;

  cout << "enter input file: " << ends;
  string inFileName;
  cin >> inFileName;

  cout << "Enter # of sigmas of change you want to make: " << ends;
  double scale;
  cin >> scale;

  string fullInFileName = inPath + inFileName;
  cout << "working on " << fullInFileName << endl;
  system(string("ls " + fullInFileName).c_str());

  ifstream inFile(fullInFileName.c_str());

  while(1) { // until end of file:
    // extract only good info from inFile into cleanOstream for one PDF:
    string cleanString;
    cleanStream(inFile, cleanString);

    cout << "clean string is \"" << cleanString << "\"" << endl;

    // Count params and ger pdf name:
    string pdfName;
    int npar = count("pdfOnResDK.", cleanString, pdfName);

    if (0 == npar) {
      break;
    }

    // pass this info on to make the new input files for this PDF:
    istrstream cleanIstream(cleanString.c_str());    
    string fileBaseName = outPath + inFileName + pdfName;
    genShapeSyst(cleanIstream, fileBaseName, npar, scale);
  }
}

