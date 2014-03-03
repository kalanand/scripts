/*
* Bdk Pdf holder
*/

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;


#include "TString.h"
#include "BToDKTo3piK/BdkEvtTypes.hh"
#include "BToDKTo3piK/BdkPdfAbsBase.hh"
#include "BToDKTo3piK/BdkPdfHolder.hh"

ClassImp(BdkPdfHolder)

//------------------------------------------------
BdkPdfHolder::BdkPdfHolder() {} //Constructor

BdkPdfHolder::~BdkPdfHolder() {} //Deconstructor

//------------------------------------------------

BdkPdfAbsBase * BdkPdfHolder::pdf(int typeIndex) {

  switch(typeIndex) {
  case BdkEvtTypes::SIG_BAD_D   : return sigBadD0()   ; break;
  case BdkEvtTypes::SIG_GOOD_D  : return sigGoodD0()  ; break;
  case BdkEvtTypes::DPi_BAD_D   : return DpiBadD0()   ; break;
  case BdkEvtTypes::DPi_GOOD_D  : return DpiGoodD0()  ; break;
  case BdkEvtTypes::DPiX        : return DPiX()       ; break;
  case BdkEvtTypes::DKX         : return DKX()  ; break;
  case BdkEvtTypes::BB_BAD_D    : return BBBadD0() ; break;
  case BdkEvtTypes::BB_GOOD_D   : return BBGoodD0(); break;
  case BdkEvtTypes::QQ_BAD_D    : return qqBadD0()  ; break;
  case BdkEvtTypes::QQ_GOOD_D   : return qqGoodD0() ; break;

  }
  if( typeIndex<0 && typeIndex>=BdkEvtTypes::NTYPES ) {
   cerr<< " Index " << typeIndex << " is out of allwoed range(0-9)."
       << " Unrecognized argument \"" << typeIndex << "\". Returning 0."
       << endl;
 
  }


  return 0;
}

BdkPdfAbsBase * BdkPdfHolder::pdf(const char * typeNameChar) {
  TString typeName(typeNameChar);
  if(0 == typeName.CompareTo("SIG_BAD_D"))    return sigBadD0()   ;
  if(0 == typeName.CompareTo("SIG_GOOD_D"))   return sigGoodD0()  ;
  if(0 == typeName.CompareTo("DPi_BAD_D") )   return DpiBadD0()   ;
  if(0 == typeName.CompareTo("DPi_GOOD_D"))   return DpiGoodD0()  ;
  if(0 == typeName.CompareTo("DPiX"))         return DPiX()      ;
  if(0 == typeName.CompareTo("DKX"))          return DKX()  ;
  if(0 == typeName.CompareTo("BB_BAD_D") )    return BBBadD0() ;
  if(0 == typeName.CompareTo("BB_GOOD_D"))    return BBGoodD0();
  if(0 == typeName.CompareTo("QQ_BAD_D"))     return qqBadD0()  ;
  if(0 == typeName.CompareTo("QQ_GOOD_D"))    return qqGoodD0() ;
  
  cerr<< " BdkPdfHolder::pdf(const char * typeNameChar ) :"
      << "    Unrecognized argument \"" << typeName << "\". Returning 0."
      << endl;
  return 0;
}

     
