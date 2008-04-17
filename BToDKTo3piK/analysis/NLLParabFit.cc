#include "../BToDKTo3piK/analysis/plotNLL.cc"

void NLLParabFit(RooRealVar * par) {
  
  cout << "plotting " ;
  par->Print();
  double val = par->getVal();
  double range = 0.6;

  plotNLL(30, 
	  val - range,
	  val + range,
	  kTRUE, kTRUE, kTRUE, kFALSE, kTRUE, "fit1d-", kFALSE, 1, par);

  val = par->getVal(); // update to fitted minimum value

  TH1 * histToFit = 0;

  const RooRealVar * xN = pdfOnResDK.xMinus();
  const RooRealVar * yN = pdfOnResDK.yMinus();
  const RooRealVar * xP = pdfOnResDK.xPlus();
  const RooRealVar * yP = pdfOnResDK.yPlus();
  
  if (xN == par) {
    histToFit = LXn;
  }
  else if (yN == par) {
    histToFit = LYn;
  }
  else if (xP == par) {
    histToFit = LXp;
  }
  else if (yP == par) {
    histToFit = LyP;
  }

  TCanvas * can = new TCanvas("NLLParabFit", "NLLParabFit", 1000, 500);
  can->Divide(2,1);
  cout << "********* narrow fit *********" << endl;
  can->cd(1);
  histToFit->Draw();
  histToFit->Fit("pol2", "", "", val - 0.06, val + 0.06); 
  can->Update();

  cout << "********* wide fit *********" << endl;
  can->cd(2);
  histToFit->Draw();
  histToFit->Fit("pol2"); 
  can->Update();

}
