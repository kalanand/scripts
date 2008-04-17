#include <iostream>
using namespace std;

#include <math.h>
#include "TLine.h"
#include "TEllipse.h"

#include "BToDKTo3piK/BdkDP.hh"
#include "BToDKTo3piK/BdkP2.hh"


double BdkP2::s23() const {
  return BdkDP::sumMasses2() -_s12 -_s13;
}

void BdkP2::setS12S23(double s12, double s23) {
  _s12 = s12;
  _s13 = BdkDP::sumMasses2() - s12 - s23;
}

void BdkP2::setS13S23(double s13, double s23) {
  _s13 = s13;
  _s12 = BdkDP::sumMasses2() - s13 - s23;
}


void BdkP2::Print(std::ostream & str) const {
  str << *this << endl;
}

ostream & operator<<(ostream & str, const BdkP2 & p2) {
  str << "BdkP2(" << p2.s12() << ", " << p2.s13() << ", " << p2.s23() << ")";
  return str;
}

Bool_t BdkP2::inDalitz() const {
  return BdkDP::inDalitz(*this);
}


void BdkP2::Draw(int fillColor, double radius, int fillStyle) const {
  TEllipse * el = new TEllipse(s12(), s13(), radius, radius);
  el->SetFillColor(fillColor);
  el->SetFillStyle(fillStyle);

  // Prevent the line from covering a small-radius fill:
  el->SetLineWidth(0);
  el->SetLineColor(fillColor);

  el->Draw();
}


bool BdkP2::operator==(const BdkP2&p) const {
  return (fabs(_s12 - p._s12) < BdkDP::tolerance() && 
	  fabs(_s13 - p._s13) < BdkDP::tolerance());
}
