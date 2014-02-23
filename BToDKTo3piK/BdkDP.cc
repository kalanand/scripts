#include "BToDKTo3piK/BdkDP.hh"

#include "TLine.h"
#include "BToDKTo3piK/BdkP2.hh"


void BdkDP::drawBoundary(int npoints, 
			 int color, 
			 int width, 
			 int style) {

  double space = edge() / npoints;
  for (int i = 0; i < npoints; ++i){
    double x1 = space * i;
    double max1 = s13Max(x1);
    double min1 = s13Min(x1);

    double x2 = space * (i+1);
    double max2 = s13Max(x2);
    double min2 = s13Min(x2);

    if (min1 != 0 && max1 != 0 && min2 != 0 && max2 != 0) {
      line(x1, max1, x2, max2, color, width, style);
      line(x1, min1, x2, min2, color, width, style);
    }
    else { // for the extreme vertical parts of the boundary:
      line(x1, max1, x1, min1, color, width, style);
      line(x2, max2, x2, min2, color, width, style);
    }
  }
}



TLine * BdkDP::line(const BdkP2 & p1, const BdkP2 & p2, 
		    int color, int width, int style) {

  return line(p1.s12(), p1.s13(), p2.s12(), p2.s13(), 
	      color, width, style);
}



TLine * BdkDP::line(double x1, double y1, double x2, double y2, 
		    int color, int width, int style) {
  TLine * result = new TLine(x1, y1, x2, y2);
  result->SetLineColor(color);
  result->SetLineWidth(width);
  result->SetLineStyle(style);
  result->Draw();
  return result;
}


TLine * BdkDP::lineEdge(double x1, double y1, double x2, double y2, 
			int color, int width, 
			int style) {
  // Draw a line covering the two points and extending all the way to
  // the axes and to edge():
  if (x1==x2) {
    // vertical line:
    if (x1 < edge() && x1 >0){
      return line(x1, 0, x1, edge(), color, width, style);
    }
  }
  else if (y1==y2) {
    // horizontal line:
    if (y1 < edge() && y1 >0){
      return line(0, y1, edge(), y1, color, width, style);
    }
  }
  else {
    // diagonal line:
    double slope = (y2-y1)/(x2-x1);
    return lineEdgeSlope(x1, y1, slope, color, width, style);
  }
  return 0;
}  



TLine * BdkDP::lineEdgeSlope(double x1, double y1, double slope, 
			     int color, int width, 
			     int style) {
  // There are four possible intercepts with the boundaries:
  BdkP2 intLeft  (0,                         y1 - slope * x1);
  BdkP2 intRight (edge(),                    y1 + slope * (edge() - x1));
  BdkP2 intBottom(x1 - y1 / slope,           0);
  BdkP2 intTop   (x1 + (edge() - y1) / slope,  edge());

  // The correct intercepts are the ones whose coordinates are between
  // 0 and edge(). 
  BdkP2 * ints[2] = {0,0};
  int index = -1;

  if (inBox(intLeft))                {ints[++index] = &intLeft;}
  if (inBox(intTop))                 {ints[++index] = &intTop;}
  if (inBox(intRight)  && index < 1) {ints[++index] = &intRight;}
  if (inBox(intBottom) && index < 1) {ints[++index] = &intBottom;}
  // Above, need to check index<1, since if the intercept is at the 
  // corner then two intercepts would be inBox()

  if (1 == index && (0 != ints[0] && 0 != ints[1])) {
    return line(*(ints[0]), *(ints[1]), color, width, style);
  }

  return 0;
}



Bool_t BdkDP::inBox(const BdkP2 & point){ 
  if (point.s12() >= 0      && point.s13() >= 0 &&
      point.s12() <= edge() && point.s13() <= edge()) {
    return kTRUE;
  }
  return kFALSE;
}



Bool_t BdkDP::inDalitz(const BdkP2 & point){ 
  const double s12 = point.s12();
  const double s13 = point.s13();

  // kinematic limits
  Double_t m13High = s13Max(s12);
  if (m13High<=0) return kFALSE;

  Double_t m13Low = s13Min(s12);
  if (m13Low<=0) return kFALSE; 

  // decide
  if ( (s13 < m13Low) || (s13 > m13High)) return kFALSE;

  return kTRUE;
}


// A locally used squaring function:
static double sqr(double d) {return d*d;}


Double_t BdkDP::s13Max(Double_t s12){
  if ((s12 < sqr(m1()+m2())) || (s12>sqr(mD()-m3()))) return 0;
  
  // Energies of particles 1 and 3 in s12 restframe
  Double_t e1star = (s12 - m2sq() + m1sq()) / (2.0*sqrt(s12));
  Double_t e3star = (mDsq() - s12 - m3sq()) / (2.0*sqrt(s12));
  
  return sqr(e1star+e3star)-sqr(sqrt(e1star*e1star-m1sq())
				-sqrt(e3star*e3star - m3sq()));
}

Double_t BdkDP::s13Min(Double_t s12){
  if ((s12<sqr(m1()+m2())) || (s12>sqr(mD()-m3()))) return 0;
  
  // Energies of particles 1 and 3 in s12 restframe
  Double_t e1star = (s12 - m2sq() + m1sq()) / (2.0*sqrt(s12));
  Double_t e3star = (mDsq() - s12 - m3sq()  ) / (2.0*sqrt(s12));
  
  return sqr(e1star+e3star)-sqr(sqrt(e1star*e1star-m1sq())
				+sqrt(e3star*e3star - m3sq()));
}


TLine * BdkDP::drawE12Line() {
  return lineEdgeSlope(symValue(), symValue(), 
		       -0.5, kGreen, 2, 2);
}


TLine * BdkDP::drawE13Line() {
  return lineEdgeSlope(symValue(), symValue(), 
		       -2, kGreen, 2, 2);
}


TLine * BdkDP::drawE23Line() {
  return lineEdgeSlope(symValue(), symValue(), 
		       1, kGreen, 2, 2);
}



