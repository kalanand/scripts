#include "BToDKTo3piK/BdkIsoOp.hh"

#include <assert.h>

ClassImp(BdkIsoOp)

const char * BdkIsoOp::name(int i) {
  switch(i) {
  case P32: return "P32";
  case P22: return "P22";
  case P12: return "P12";
  case P21: return "P21";
  case P11: return "P11";
  case P01: return "P01";
  case P10: return "P10";
  default: return 0;
  }
}


BdkIsoOp::Op BdkIsoOp::op(int i) {
  switch(i) {
  case P32: return P32;
  case P22: return P22;
  case P12: return P12;
  case P21: return P21;
  case P11: return P11;
  case P01: return P01;
  case P10: return P10;
  default: assert(0);
  }
}
