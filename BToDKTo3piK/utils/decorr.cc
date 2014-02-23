// $Id: decorr.cc,v 1.1 2006/04/20 16:56:23 fwinkl Exp $

// Decorrelate var1 and var2 in data.
// Adds new columns called var1_p and var2_p with the uncorrelated values.
void decorr(RooDataSet& data, RooRealVar& var1, RooRealVar& var2)
{
  TH2* h2 = data.createHistogram(var1,var2);

  Double_t rho = h2->GetCorrelationFactor();
  Double_t s1 = sqrt(h2->GetCovariance(1,1));
  Double_t s2 = sqrt(h2->GetCovariance(2,2));

  RooRealVar theta("theta","",0.5*atan(2*rho*s1*s2/(s1*s1-s2*s2)));
  
  RooFormulaVar var1p(TString(var1.GetName())+"_p",
                      TString(var1.GetTitle())+" decorrelated",
                      "cos(@0)*@1+sin(@0)*@2",
                      RooArgList(theta,var1,var2));

  RooFormulaVar var2p(TString(var2.GetName())+"_p",
                      TString(var2.GetTitle())+" decorrelated",
                      "-sin(@0)*@1+cos(@0)*@2",
                      RooArgList(theta,var1,var2));

  data.addColumn(var1p);
  data.addColumn(var2p);

}
