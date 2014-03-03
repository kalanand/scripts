#include <assert.h>
#include <math.h>
#include "BToDKTo3piK/BdkSymOp.hh"
#include "BToDKTo3piK/BdkBinWgtList.hh"

ClassImp(BdkBinWgtList)

BdkBinWgtList::BdkBinWgtList(int op) :
  _op (BdkIsoOp::op(op)){
}


void BdkBinWgtList::setI(const BdkTriBin & bin) {
  _i.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _i.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P22: _i.setWeight(0.5); break;
  case BdkIsoOp::P12: _i.setWeight(3./sqrt(44.)); break;
  case BdkIsoOp::P21: _i.setWeight(1./sqrt(12.)); break;
  case BdkIsoOp::P11: _i.setWeight(0.5); break;
  case BdkIsoOp::P01: _i.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _i.setWeight(1./sqrt(2.)); break;
  }
}


void BdkBinWgtList::setR(const BdkTriBin & bin) {
  _r.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _r.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P22: _r.setWeight(-0.5); break;
  case BdkIsoOp::P12: _r.setWeight(3./sqrt(44.)); break;
  case BdkIsoOp::P21: _r.setWeight(1./sqrt(12.)); break;
  case BdkIsoOp::P11: _r.setWeight(-0.5); break;
  case BdkIsoOp::P01: _r.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _r.setWeight(0.); break;
  }
}


void BdkBinWgtList::setRR(const BdkTriBin & bin) {
  _rr.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _rr.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P22: _rr.setWeight(0.); break;
  case BdkIsoOp::P12: _rr.setWeight(-2./sqrt(44.)); break;
  case BdkIsoOp::P21: _rr.setWeight(-2./sqrt(12.)); break;
  case BdkIsoOp::P11: _rr.setWeight(0.); break;
  case BdkIsoOp::P01: _rr.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _rr.setWeight(0.); break;
  }
}


void BdkBinWgtList::setE(const BdkTriBin & bin) {
  _e.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _e.setWeight(-1./sqrt(6.)); break;
  case BdkIsoOp::P22: _e.setWeight(-0.5); break;
  case BdkIsoOp::P12: _e.setWeight(-3./sqrt(44.)); break;
  case BdkIsoOp::P21: _e.setWeight(1./sqrt(12.)); break;
  case BdkIsoOp::P11: _e.setWeight(0.5); break;
  case BdkIsoOp::P01: _e.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _e.setWeight(-1./sqrt(2.)); break;
  }
}


void BdkBinWgtList::setER(const BdkTriBin & bin) {
  _er.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _er.setWeight(-1./sqrt(6.)); break;
  case BdkIsoOp::P22: _er.setWeight(0.5); break;
  case BdkIsoOp::P12: _er.setWeight(-3./sqrt(44.)); break;
  case BdkIsoOp::P21: _er.setWeight(1./sqrt(12.)); break;
  case BdkIsoOp::P11: _er.setWeight(-0.5); break;
  case BdkIsoOp::P01: _er.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _er.setWeight(0.); break;
  }
}

 
void BdkBinWgtList::setERR(const BdkTriBin & bin) {
  _err.setBin(&bin);

  switch(_op) {
  case BdkIsoOp::P32: _err.setWeight(-1./sqrt(6.)); break;
  case BdkIsoOp::P22: _err.setWeight(0.); break;
  case BdkIsoOp::P12: _err.setWeight(2./sqrt(44.)); break;
  case BdkIsoOp::P21: _err.setWeight(-2./sqrt(12.)); break;
  case BdkIsoOp::P11: _err.setWeight(0.); break;
  case BdkIsoOp::P01: _err.setWeight(1./sqrt(6.)); break;
  case BdkIsoOp::P10: _err.setWeight(0.); break;
  }
}


const BdkBinWgt & BdkBinWgtList::binWgt(int op) const {
  switch(op) {
  case BdkSymOp::I: return _i;
  case BdkSymOp::R: return _r;
  case BdkSymOp::RR: return _rr;
  case BdkSymOp::E: return _e;
  case BdkSymOp::ER: return _er;
  case BdkSymOp::ERR: return _err;
  default: assert(0);
  }
}

const BdkBinWgt * BdkBinWgtList::binWgt(const BdkTriBin & bin) const {
  // Find the binWgt that points to bin:
  for (int op = 0; op < BdkSymOp::NOPS; ++op) {
    const BdkBinWgt & bw = binWgt(op);
    if (bw.bin() == &bin) {
      return &bw;
    }
  }
  return 0;
}
