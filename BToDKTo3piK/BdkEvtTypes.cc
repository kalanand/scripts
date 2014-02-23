#include "BToDKTo3piK/BdkEvtTypes.hh"

ClassImp(BdkEvtTypes)

const char * BdkEvtTypes::name(int t) {
  const char * result = "UnknownType";
  switch(t) {
  case SIG_BAD_D: result = "SigBadD"; break;
  case SIG_GOOD_D: result = "SigGoodD"; break;
  case DPi_BAD_D: result = "DpiBadD"; break;
  case DPi_GOOD_D: result = "DpiGoodD"; break;
  case DPiX: result = "DpiX"; break;
  case DKX: result = "DKX"; break;
  case BB_BAD_D: result = "BbBadD"; break;
  case BB_GOOD_D: result = "BbGoodD"; break;
  case QQ_BAD_D: result = "QqBadD"; break;
  case QQ_GOOD_D: result = "QqGoodD"; break;
  }

  return result;
}
