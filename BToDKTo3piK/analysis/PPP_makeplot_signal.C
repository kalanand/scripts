{
  gROOT->Reset();
  gStyle->Reset();
  gROOT->SetStyle("BABAR");

  int numBINS = 100;

  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0}) [GeV^{2}/c^{4}]",0,3); //This is pi+ pi0  
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]",0,3); //This is pi- pi0  
  fit_S23.setBins(numBINS);
  fit_S31.setBins(numBINS);
 

  BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
  dalitzCfg->getParameters(RooArgSet())->readFromFile("../BToDKTo3piK/Dstar/dalitzCfg.par");
  dalitzCfg->getParameters(RooArgSet())->Print("v");


  // define signal pdf
  BdkPdf2DpolyDalitz eff("eff","",fit_S23, fit_S31);
  eff.parameters().readFromFile("PPP_eff.par");
  BdkPdfDDalitz mypdf("mypdf","mypdf",fit_S23,fit_S31,BdkDalitzBase::D0);
  dstar.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());

  RooArgSet allPars = mypdf.parameters();
  mypdf.parameters().readFromFile("PPP_NominalFit.par");


  RooDataSet *fitdata = mypdf.generate(RooArgSet(fit_S23,fit_S31),100000);

  //Make the plots
  //******************************************************
  TH2F* hFit = fitdata->createHistogram(fit_S31,fit_S23);
  hFit->Sumw2();
  hFit->GetXaxis()->SetTitle("m^{2}_{#pi^{-}#pi^{0}}");
  hFit->GetYaxis()->SetTitle("m^{2}_{#pi^{+}#pi^{0}}");
  hFit->SetTitle("Data");


  TCanvas* aCanvas = new TCanvas("c3", "c3", 900, 600);
  gStyle->SetOptStat(0);
  hFit->Draw("colz");
  aCanvas->Draw();
  aCanvas->SaveAs(TString("ppp_chi_nominal.gif"));
  aCanvas->SaveAs(TString("ppp_chi_nominal.eps"));
  aCanvas->SaveAs(TString("ppp_chi_nominal.C"));
}  //end the macro







