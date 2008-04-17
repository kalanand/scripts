// $Id: printFitResult.cc,v 1.7 2007/03/02 09:48:27 fwinkl Exp $
// prints fit result with error and cov matrices

// #include <iomanip.h>
#include "TMatrix.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooFitResult.hh"

void printCovMatrix(RooFitResult * res, ostream & stream);
void getCovCorMatrix(RooFitResult* res, TMatrix& cov, TMatrix& cor);

void printCorMatrix(const TMatrix& cor, ostream& stream = cout,
                    Bool_t printLatex = kFALSE, const char* fmt = "%.3f");

void printCorMatrix(RooFitResult* res, ostream& stream = cout,
                    Bool_t printLatex = kFALSE, const char* fmt = "%.3f");

//--------------------------------------------------------------------
void printChanges(RooFitResult * res, ostream & os = cout) 
{
  for (int i = 0 ; i < res->floatParsFinal().getSize(); ++i) {
  
    RooRealVar *final = (RooRealVar*)res->floatParsFinal().at(i);
    
    os << "Change in " << setw(20) << final->GetName();
    const double init  = ((RooRealVar*)res->floatParsInit().at(i))->getVal();
    double diff = final->getVal() - init;
    
    double error = -1;
    if (final->hasError())
      error = ((RooRealVar*)res->floatParsFinal().at(i))->getError();
    else if (final->hasAsymError()) {
      if (diff<0) error = fabs(final->getAsymErrorLo());
      else error = fabs(final->getAsymErrorHi());
    }
    if (error > 0) diff /= error;
    else diff = -999;

    os << " = " << diff << endl;
  }
}

//--------------------------------------------------------------------
void printFitResult(RooFitResult * res) {
  if (0 == res) {
    return;
  }

  TMatrix cov;
  TMatrix cor;

  res->Print();
  cout << endl << "Fit status = "<<res->status()<<endl<<endl;

  getCovCorMatrix(res, cov, cor);
  cout << "covariance matrix" << endl;
  cov.Print();
  cout << "correlation matrix" << endl;
  cor.Print();

  printChanges(res);
}


//--------------------------------------------------------------------
void printFitResult(RooFitResult * res, ostream & stream) {
  if (0 == res) {
    return;
  }

  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);

  printCovMatrix(res, stream);
  printCorMatrix(cor, stream);

  res->printToStream(stream, RooPrintable::Verbose);
  stream << endl << "Fit status = "<<res->status()<<endl<<endl;
  printChanges(res, stream);

}

//--------------------------------------------------------------------
// Print the correlation matrix
void printCorMatrix(RooFitResult* res, ostream& stream, 
                    Bool_t printLatex, const char* fmt)
{
  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);
  printCorMatrix(cor, stream, printLatex, fmt);
}


void printCorMatrix(const TMatrix& cor, ostream& stream, 
                    Bool_t printLatex, const char* fmt)
{
  for (int i = 0; i < cor.GetNrows(); ++i) {
    for (int j = 0; j < cor.GetNcols(); ++j) {
      TString s;
      if (cor(i,j)==1) s = "1";
      else s.Form(fmt,(double)cor(i,j));

      if (printLatex) {
        if (j>0) stream << " & ";
        stream << s;      
      }
      else stream << setw(13) << s;
    }
    if (printLatex) stream << "\\\\";
    stream << endl;
  }
}

//--------------------------------------------------------------------
// Print matrix
void printMatrix(TMatrix& cov, ostream & stream,
       	   int precision = 8, int width = 15) {

  stream << "//$" << endl;
  for (int i = 0; i < cov.GetNrows(); ++i) {
    stream << "// ";
    for (int j = 0; j < cov.GetNcols(); ++j) {
      stream << setprecision(precision) << setw(width) << cov(i,j);
    }
    stream << endl;
  }
  stream << "//$" << endl << endl;
}

//--------------------------------------------------------------------
// Print the covariance matrix
void printCovMatrix(RooFitResult * res, ostream & stream) {

  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);

  printMatrix(cov, stream);
}


//--------------------------------------------------------------------
// Fill the correlation and covariance matrix from the fit result
void getCovCorMatrix(RooFitResult* res, TMatrix& cov, TMatrix& cor)
{
  if (!res) return;
  const int size = res->floatParsFinal().getSize();
  cov.ResizeTo(size, size);
  cor.ResizeTo(size, size);
  for (int i = 0; i < size; ++i) {
    RooRealVar * vari = (RooRealVar *)(res->floatParsFinal().at(i));
    for (int j = 0; j < size; ++j) {
      RooRealVar * varj = (RooRealVar *)(res->floatParsFinal().at(j));
      cor(i,j) = res->correlation(*vari, *varj);
      cov(i,j) = cor(i,j) * vari->getError() * varj->getError();
    }
  }
}
