#include "BToDKTo3piK/BdkBinWgt.hh"
#include "BToDKTo3piK/BdkTriBin.hh"

ClassImp(BdkBinWgt)


BdkBinWgt::BdkBinWgt() :
    _bin(0), 
    _weight(0) {
}


BdkBinWgt::BdkBinWgt(const BdkTriBin * bin, double weight) :
    _bin(bin), 
    _weight(weight) {
}


// The copy constructor and equality operator copy the pointer and weight:
BdkBinWgt::BdkBinWgt(const BdkBinWgt & source) {
  *this = source;
}


BdkBinWgt & BdkBinWgt::operator=(const BdkBinWgt & source) {
  // there is no bin to point to after a copy operation.  The BdkBinning's 
  // buildWgtLists() has to be called again
  _bin = 0; 
  _weight = 0;
  return *this;
}
