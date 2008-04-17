// $Id: toyErrors.cc,v 1.1 2006/05/29 22:40:34 fwinkl Exp $
// Fit toy MC sample with different variables in the fit to study errors

void toyErrors()
{
  gROOT->cd();
  data = pdfOnResDK.generate(10000);

  pdfOnResDK.useDalitz(true);
  pdfOnResDK.useDE(true);
  pdfOnResDK.useNnCont(true);
  pdfOnResDK.useNnComb(false);

  fitOption = "mer";

  fit(pdfOnResDK,*data);
  RooArgSet set1((RooArgSet&)*pdfOnResDK.parametersFree().snapshot(false));
  
  pdfOnResDK.useDalitz(true);
  pdfOnResDK.useDE(false);
  pdfOnResDK.useNnCont(false);
  pdfOnResDK.useNnComb(false);

  fit(pdfOnResDK,*data);
  RooArgSet set2((RooArgSet&)*pdfOnResDK.parametersFree().snapshot(false));

  pdfOnResDK.useDalitz(false);
  fixxy();
  pdfOnResDK.useDE(true);
  pdfOnResDK.useNnCont(true);
  pdfOnResDK.useNnComb(false);

  fit(pdfOnResDK,*data);
  RooArgSet set3((RooArgSet&)*pdfOnResDK.parametersFree().snapshot(false));

  set1.Print("v");
  set2.Print("v");
  set3.Print("v");

}
