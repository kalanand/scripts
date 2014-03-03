// $Id: printFrac.cc,v 1.1 2006/06/25 01:30:05 fwinkl Exp $
// Print the event fractions of BdkPdfOnRes

void printFrac(const BdkPdfOnRes& pdf = pdfOnResDK)
{
  RooArgSet evts(pdf.numEvt());
  Double_t total = pdf.totalNumEvts();

  TIterator *iter = evts.createIterator();
  RooAbsReal* r;
  while (r = (RooAbsReal*)iter->Next()) {
    Double_t f = r->getVal()/total*100;
    TString s;
    s.Form("%.1f",f);
    cout << setw(30) << r->GetName() << " = " << s << endl;
  }
}
