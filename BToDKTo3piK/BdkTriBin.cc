#include <iostream>
using namespace std;

#include "TLine.h"
#include "TPaveText.h"
#include "TRandom3.h"

#include "BToDKTo3piK/BdkDP.hh"
#include "BToDKTo3piK/BdkTriBin.hh"

ClassImp(BdkTriBin)

// STATIC MEMBERS:
int BdkTriBin::calcAreaNSteps() {
  return setCalcAreaNSteps(-1, kFALSE);
}


int BdkTriBin::setCalcAreaNSteps(int n, Bool_t change) {
  static int result = 1000;
  if (change) {
    result = n;
  }
  return result;
}


double BdkTriBin::effCoeffS(int index) {
  return setEffCoeffS(index);
}

double BdkTriBin::setEffCoeffS(int index, double val, Bool_t change) {
  assert(index >=0 && index <= 5);

  // The default values come from Kalanand's entries in the BAD note:
  static double coeffS[6] = {
    1.0,
    35.6,
    -11.8,
    0.79,
    3.6,
    -23.5
  };

  if (change) {
    coeffS[index] = val;
  }

  return coeffS[index];
}
    
  
double BdkTriBin::effCoeffSErr(int index) {
  assert(index >=0 && index <= 5);

  // The values come from Kalanand's entries in the BAD note:
  static const double coeffSErr[6] = {
    0.0,  // c0 is fixed
    5.8,
    2.0,
    0.47,
    1.3,
    5.8
  };

  return coeffSErr[index];
}

    


// INSTANCE MEMBERS:
BdkTriBin::BdkTriBin() :
  _waveFunc(0),
  _waveFuncErr(0),
  _nEventsSR(0),
  _nEventsSB(0),
  _isMirror(kFALSE),
  _area(0),
  _inDalitzFrac(0)
{
  initWgtLists();
}


// This constructor sets 
// p1 at the bottom left (bottom right) of the triangle for 
// mirror == kFALSE (kTRUE),
// p2 above p1,
// p3 to the right (top left) of p1 for mirror == kFALSE (kTRUE):
BdkTriBin::BdkTriBin(const BdkP2 & p1, double space, Bool_t mirror, 
		     const char * name) :
  _name(name),
  _waveFunc(0),
  _waveFuncErr(0),
  _nEventsSR(0),
  _nEventsSB(0),
  _isMirror(mirror)
{
  initWgtLists();

  BdkP2 p2(p1.s12(),           p1.s13() + space);
  BdkP2 p3(p1.s12() + space,   p1.s13());
  
  if (kTRUE == mirror) { // Put p3 to the left of and above p1:
    p3 = BdkP2(p1.s12() - space,   p1.s13() + space);
  }
  
  setPoints(p1, p2, p3);
}
  
  


BdkTriBin::BdkTriBin(const BdkP2 & p1, const BdkP2 & p2, const BdkP2 & p3, 
		     const char * name) :
  _name(name),
  _waveFunc(0),
  _waveFuncErr(0),
  _nEventsSR(0),
  _nEventsSB(0)
{
  initWgtLists();
  setPoints(p1, p2, p3);
}



BdkTriBin::~BdkTriBin() {}


void BdkTriBin::initWgtLists() {
  for (int op = 0; op < BdkIsoOp::NOPS; ++op) {
    _wgtList[op] = BdkBinWgtList(op);
  }
}


void BdkTriBin::setPoints(const BdkP2 &p1, const BdkP2 &p2, const BdkP2 &p3) {
  _p1 = p1;
  _p2 = p2;
  _p3 = p3;
  calcMinMax();
  calcArea();
}



void BdkTriBin::calcWaveFunc(double mDratioSRtoSB) {
  // the factorSB is the product of mDratioSRtoSB and the ratio of 
  // the bin's area() to the area of the corresponding sideband bin.
  // At this point we take this ratio to be 1, so:
  double factorSB = mDratioSRtoSB; 
  _waveFunc = 
    sqrt((_nEventsSR - _nEventsSB * factorSB) 
	 / efficiency() / area()
	 );
  
  _waveFuncErr = 1.0 / (2.0 * _waveFunc * area() * efficiency())
    * sqrt(_nEventsSR + factorSB * factorSB * _nEventsSB);
}
  


const BdkP2 BdkTriBin::center() const {
  double x = (_p1.s12() + _p2.s12() + _p3.s12())/3;
  double y = (_p1.s13() + _p2.s13() + _p3.s13())/3;
  return BdkP2(x, y);
}


BdkTriBin::InDalitzCode  BdkTriBin::inDalitzStat() const {
  Bool_t p1In = BdkDP::inDalitz(_p1);
  Bool_t p2In = BdkDP::inDalitz(_p2);
  Bool_t p3In = BdkDP::inDalitz(_p3);

  if (kFALSE == p1In &&
      kFALSE == p2In &&
      kFALSE == p3In) {
    return NONE;
  }
  
  if (kFALSE == p1In ||
      kFALSE == p2In ||
      kFALSE == p3In) {
    return SOME;
  }
  
  return ALL;
}



void BdkTriBin::calcArea(Bool_t verbose, Bool_t drawPoints) {
  if (kTRUE == verbose) {
    cout << "calcArea() called for bin " << *this << endl;
  }
  
  // First, calculate the "naive" area:
  _area = (_max.s12() - _min.s12()) * (_max.s13() - _min.s13()) / 2;

  if (0 == _area) {
    return;
  }

  // check if the bin is in/out of the DP:
  InDalitzCode in = inDalitzStat();
  switch(in) {
  case NONE:
    if (kTRUE == verbose) {
      cout << "bin is outside the DP, area=" << 0 << endl;
    }
    _inDalitzFrac = 0;
    _area = 0;    // the entire bin is outside the DP
    return;      
  case ALL:
    if (kTRUE == verbose) {
      cout << "bin is inside the DP, area=" << _area << endl;
    }
    _inDalitzFrac = 1.0;
    return;    // the entire bin is inside the DP
  default:
    // calculate the fraction of the bin that's within the DP. We
    // sample points inside the bin and check what fraction of points
    // is inside the DP.

    int nPointsIn = 0;
    /*
      // calculate using random points. This is much slower:
      for (int i = 0; i < calcAreaNPoints(); ++i) {
      BdkP2 point = generatePoint();
      if (BdkDP::inDalitz(point)){
      ++nPointsIn; // count # of points inside the DP
      }      
    }
    */
    int nPoints = 0;
    const int nSteps = calcAreaNSteps();
    const double step = (_max.s12() - _min.s12()) / nSteps;

    // Do the loop on the x points:
    for (int i = 1; i <= nSteps; ++i) {
      double x = _min.s12() + i * step; // first assume not mirror
      if (_isMirror) {  
	x = _max.s12() - i * step;      // from right to left for a mirror bin
      }
      // The y loop. The max on j is easy to understand if one thinks
      // of a non-mirror bin with x starting from the left:
      for (int j = 1; j <= nSteps - i - 1; ++j) { 
	++nPoints;   // count # of points checked
	double y = _min.s13() + j * step;  // for a non-mirror bin
	if (_isMirror) {
	  y = _max.s13() - j * step; // from top to bottom for a mirror bin
	}

	BdkP2 point(x,y);
	if (BdkDP::inDalitz(point)) {
	  ++nPointsIn;  // count this point as inside the DP
	}       	    	  

	if (drawPoints) {
	  // Draw the point:
	  if (BdkDP::inDalitz(point)) {
	    point.Draw(kBlack);
	  }
	  else {
	    point.Draw(kRed);
	  }
	} // end if to draw the point
      } // end loop on j (=y)
    } // end loop on i (=x)

    _inDalitzFrac = (float)nPointsIn / nPoints;

    if (verbose) {
      cout << "A fraction " << _inDalitzFrac << " of bin is inside the DP."
	   << endl;
    }

    _area *= _inDalitzFrac;
  }
}



void BdkTriBin::writeName(double relSize) const {
  if (0 == _name.Length()) {
    return;
  }

  const double size = (_max.s12() - _min.s12()) * relSize;

  BdkP2 cen = center();
  double x = cen.s12();
  double y = cen.s13();

  TPaveText * text = new TPaveText(x - size, y - size, x + size, y + size);

  text->AddText(_name);
  text->SetLineWidth(0);
  text->SetFillColor(10);
  text->SetMargin(0);
  text->Draw();
}



void BdkTriBin::Draw(int color, int width, Bool_t tagLines,
		     Bool_t doWriteName) const {
  int solid = 1;
  int dashed = 2;
  int dotted = 3;

  if (kFALSE == tagLines){
    dashed = solid;
    dotted = solid;
  }

  BdkDP::line(_p1, _p2, color, width, solid);
  BdkDP::line(_p1, _p3, color, width, dashed);
  BdkDP::line(_p2, _p3, color, width, dotted);

  if (doWriteName) {
    writeName();
  }
}


BdkTriBin BdkTriBin::R(Bool_t createName) const {
  BdkTriBin result(_p1.R(), _p2.R(), _p3.R());
  if (kTRUE == createName) {
    TString name = _name + ".R";
    result.setName(name);
  }
  
  return result;
}


BdkTriBin BdkTriBin::E(Bool_t createName) const {
  BdkTriBin result(_p1.E(), _p2.E(), _p3.E());
  if (kTRUE == createName) {
    TString name = _name + ".E";
    result.setName(name);
  }
  
  return result;
}


BdkTriBin BdkTriBin::RR(Bool_t createName) const { 
  return R(createName).R(createName);
}


BdkTriBin BdkTriBin::ER(Bool_t createName) const { 
  return E(createName).R(createName);
}


BdkTriBin BdkTriBin::ERR(Bool_t createName) const { 
  return E(createName).RR(createName);
}


bool BdkTriBin::operator==(const BdkTriBin & bin) const {
  // Bins are defined as identical even if the order of their points
  // is not the same. So it's enough to check their center point and size:
  return (center() == bin.center()  &&
	  fabs(area() - bin.area()) < BdkDP::tolerance());
}



std::vector<BdkTriBin> BdkTriBin::divide(bool createName) const {
  // Divide into four bins:
  // Start by finding the middle points along each line:
  BdkP2 mid12((_p1.s12() + _p2.s12())/2, (_p1.s13() + _p2.s13())/2);
  BdkP2 mid13((_p1.s12() + _p3.s12())/2, (_p1.s13() + _p3.s13())/2);
  BdkP2 mid23((_p2.s12() + _p3.s12())/2, (_p2.s13() + _p3.s13())/2);

  TString name1, name2, name3, name4;
  if (kTRUE == createName) {
    name1 = _name + ".div1";
    name2 = _name + ".div2";
    name3 = _name + ".div3";
    name4 = _name + ".div4";
  }

  // Then use them to create the bins:
  BdkTriBin bin1(_p1, mid12, mid13, name1);
  BdkTriBin bin2(_p2, mid12, mid23, name2);
  BdkTriBin bin3(_p3, mid13, mid23, name3);
  BdkTriBin bin4(mid12, mid13, mid23, name4);

  std::vector<BdkTriBin> result;
  result.push_back(bin1);
  result.push_back(bin2);
  result.push_back(bin3);
  result.push_back(bin4);

  return result;
}




void BdkTriBin::calcMinMax() {
  // Find the square boundaries of the triangle:
  const BdkP2 * points[3] = {&_p1, &_p2, &_p3};

  // start by setting _min and _max to _p1
  _min = _p1;
  _max = _p1;
  
  int p;
  for (p = 1; p < 3; ++p){ // then loop over _p2 and _p3 only
    if (_max.s12() < points[p]->s12()) {
      _max.setS12(points[p]->s12());
    }
    if (_max.s13() < points[p]->s13()) {
      _max.setS13(points[p]->s13());
    }
    if (_min.s12() > points[p]->s12()) {
      _min.setS12(points[p]->s12());
    }
    if (_min.s13() > points[p]->s13()) {
      _min.setS13(points[p]->s13());
    }
  }

  // See if it's a mirrored triangle, meaning that only one point has the 
  // minimum x, and two points have the maximum y:
  int nMinX = 0;
  for (p = 0; p < 3; ++p){ // now loop over all 3 points
    if (fabs(points[p]->s12() - _min.s12()) < BdkDP::tolerance()) {
      ++nMinX;
    }
  }

  switch(nMinX) {
  case 1: _isMirror = kTRUE; break;
  case 2: _isMirror = kFALSE; break;
  default:
    std::cerr 
      << "**** ERROR: BdkTriBin::calcMinMax(): Can't determine if mirrored." 
      << std::endl
      << *this << std::endl;
    assert(0);
  }

  // Check that the lengths of the 2 sides are equal.
  // Since min and max may be the result of a floating-point
  // calculation (if the bin is the result of a symmetry operation
  // on another bin), need to allow for finite resolution:
  if (fabs((_max.s12() - _min.s12()) -
	   (_max.s13() - _min.s13())) > BdkDP::tolerance()) {
    std::cerr 
      << "**** ERROR: BdkTriBin::calcMinMax(): Two sides are not equal:"
      << _max.s12() - _min.s12() << " != " << _max.s13() - _min.s13()
      << std::endl
      << *this << std::endl;
    assert(0);
  }
}



Bool_t BdkTriBin::contains(const BdkP2 & point) const {
  // First see if the point is in the enclosing square:
  if (point.s12() > _max.s12() || point.s12() < _min.s12() ||
      point.s13() > _max.s13() || point.s13() < _min.s13()) {
    return kFALSE;
  }

  // If we're here, it's in the square. See on which side of the 
  // diagonal line it is:
  if (_max.s13() - point.s13() < point.s12() - _min.s12()) {
    // Then point is above the diagonal, and is contained if _isMirror:
    if (_isMirror) {
      return kTRUE;
    }
    else {
      return kFALSE;
    }
  }
  else {
    // Then point is below the diagonal, and is contained if !_isMirror:
    if (_isMirror) {
      return kFALSE;
    }
    else {
      return kTRUE;
    }
  }
}


void BdkTriBin::Print(std::ostream & ostr) const {
  ostr << *this << endl;
}


ostream & operator<<(ostream & ostr, const BdkTriBin & bin) {
  ostr << "BdiTriBin \"" << bin.name() << "\":" << endl
       << "  p1="
       << bin.p1() << ", p2=" << bin.p2() << ", p3=" << bin.p3() << endl
       << "  min=" << bin.min() << ", max=" << bin.max() << ", isMirror=";

  if (bin.isMirror()) {
    ostr << "true";
  }
  else {
    ostr << "false";
  }
  
  ostr << endl << "  waveFunc=" << bin.waveFunc() 
       << ", waveFuncErr=" << bin.waveFuncErr()
       << ", nEventsSR=" << bin.nEventsSR() 
       << ", nEventsSB=" << bin.nEventsSB() 
       << endl
       << "  inDalitzFrac=" << bin.inDalitzFrac()
       << ", area=" << bin.area() << endl;

  return ostr;
}




BdkP2 BdkTriBin::generatePoint() const {
  static TRandom3 rand;  // static so it continues the seed sequence
  // Generate a point inside the covering box and return it if it is 
  // contained in the triangle. This has an efficiency of 50%, but it
  // is necessary for preventing a bias in the point density:
  while(1) {
    double x = rand.Uniform(_min.s12(), _max.s12());
    double y = rand.Uniform(_min.s13(), _max.s13());
    BdkP2 point(x,y);
    if (contains(point)) {
      return point;
    }
  }
}


double BdkTriBin::efficiency() const {
  BdkP2 cen = center();
  const double x = cen.s12(); // define for easier reading
  const double y = cen.s13();
  
  double value 
    = effCoeffS(0) 
    + effCoeffS(1) * (x + y)
    + effCoeffS(2) * (x*x + y*y)
    + effCoeffS(3) * (x*x*x + y*y*y)
    + effCoeffS(4) * (x*x * y + y*y * x)
    + effCoeffS(5) * x * y;
  
  return value;
}
  
