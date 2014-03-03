// $Id: chi2test2d.cc,v 1.2 2006/02/28 02:32:21 fwinkl Exp $
// Script to calculate the chi2 histogram of two 2D histograms

TVector3 chi2test2d(TH2* h1, TH2* h2, TH2* hchi2, TH1* hpull = 0)
{
  if (!h1 || !h2 || !hchi2) return 0;

  TH2* hn1 = h1->Clone("hn1");
  TH2* hn2 = h2->Clone("hn2");
  
  hn1->Sumw2();
  hn2->Sumw2();

  hn1->Scale(1/hn1->Integral());
  hn2->Scale(1/hn2->Integral());
  
  int ndof = hchi2->GetNbinsX() * hchi2->GetNbinsY();
  ndof--;   // histos are normalized to each other

  Double_t chi2Total = 0;    // total chi2
  for (int x=1; x <= hchi2->GetNbinsX(); x++) {
    for (int y=1; y <= hchi2->GetNbinsY(); y++) {

      Double_t bin1 = hn1->GetBinContent(x,y);
      Double_t bin2 = hn2->GetBinContent(x,y); 

      if (bin1==0 && bin2==0) ndof--;    // no data -> one less dof
      else {
        Double_t sqrtChi2 = (bin1-bin2);
        sqrtChi2 /= sqrt(hn1->GetBinError(x,y)**2 + hn2->GetBinError(x,y)**2);        
        hchi2->SetBinContent(x,y,sqrtChi2);
        chi2Total += sqrtChi2**2;
      }
        if (hpull) hpull->Fill(sqrtChi2);
    }
  }
  Double_t chi2Prob = TMath::Prob(0.5*chi2Total,Int_t(0.5*ndof));

  return TVector3(chi2Total,ndof,chi2Prob);
}
