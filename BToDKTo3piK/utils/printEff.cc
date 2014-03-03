
#ifndef PRINTEFF_CC
#define PRINTEFF_CC

#include "TString.h"
#include "BToDKTo3piK/BdkPdfDalitzBase.hh"
#include "BToDKTo3piK/BdkEvtTypes.hh"
#include "BToDKTo3piK/BdkDalitzEff.hh"
#include "../BToDKTo3piK/globals/globals.hh"

void printEff(Int_t charge = 1)
{
  for (int i=0; i<BdkEvtTypes::NTYPES; i++) {
    BdkPdfDalitzBase* pdf;

    if (charge>0) pdf = (BdkPdfDalitzBase*)pdfOnResDK.prodP(i)->getDalitzPdf();
    else pdf = (BdkPdfDalitzBase*)pdfOnResDK.prodN(i)->getDalitzPdf();

    TString effName = "none";
    if (pdf->efficiencyFunc()) effName = pdf->efficiencyFunc()->GetName();

    cout << setw(30) << pdf->GetName()
         << setw(20) << effName
         << endl;
  }

}

#endif
