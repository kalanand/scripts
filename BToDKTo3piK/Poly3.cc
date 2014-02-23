

#include "TF1.h"
#include "BToDKTo3piK/Poly3.hh"

ClassImp(Poly3)

TF1 * Poly3::tf1() const {
  TString formula;
  formula += _p0;
  formula += "+";
  formula += _p1;
  formula += "*x +";
  formula += _p2;
  formula += "*x*x +";
  formula += _p3;
  formula += "*x*x*x";
  
  TF1* result = new TF1(formula, formula);
  result->SetRange(0.45, 1.25);
  
  return result;
}

void Poly3::draw() const {
  TF1 * func = tf1();
  func->Draw();
}

