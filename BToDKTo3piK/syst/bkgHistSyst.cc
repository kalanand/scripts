// $Id: bkgHistSyst.cc,v 1.5 2006/07/11 21:44:03 fwinkl Exp $
//
// Evaluate impact of background histogram shapes on fit result using toy MC


TH2* hDiff = 0;
TH2* hDefault = 0;
TH2* hMod = 0;

BdkPdfDalitzHist* dalitzHist = 0;


void bkgHistSyst(Int_t samples = 1, const BdkEvtTypes::Type evtType)
{
  const Int_t NCFG = 2;
  RooArgList* fitResults[NCFG];

  loadHistogram(evtType);

  //  increaseDiff(*hDiff,3);

  // Create lists for fit results
  for (int cfg=0; cfg<NCFG; cfg++) fitResults[cfg] = new RooArgList();

  for (;samples>0;samples--) {
    if (data) delete data;
    changePdf(1,evtType);
    cout << "Generating sample "<<samples<<" ..."<<endl;
    pdfOnResDK.setNsigAsymFromXY();
    data = pdfOnResDK.generate();
    if (samples==1) pdfOnResDK.parameters().writeToFile("bkgShapeSyst-init.par");

    for (int cfg=0; cfg<NCFG; cfg++) {
      changePdf(cfg,evtType);
      RooFitResult* fit = fit(pdfOnResDK, *data, true);
      fitResults[cfg]->addClone(fit->floatParsFinal());

      TString file = "bkgShapeSyst-";
      file += cfg;
      file += "-";
      file += samples;
      file += ".fit";

      ofstream of(file);
      printFitResult(fit, of);
      of.close();
    }
  }

  changePdf(0,evtType);

  for (int cfg=0; cfg<NCFG; cfg++) {    
    cout << "Mean fit result for PDF config "<<cfg<<endl;
    RooArgSet* m = mean(*fitResults[cfg]);
    m->Print("v");
    TString file = "bkgShapeSyst-";
    file += cfg;
    file += ".par";
    m->writeToFile(file);
  }
}

void loadHistogram(const BdkEvtTypes::Type evtType)
{
  // Load histogram with MC/data ratio from sidebands
  // Create this histogram with syst/diffDataMC.cc

  if (evtType==BdkEvtTypes::DPiX) 
    dalitzHist = &dalitzHolderN.DPiXType();
  else if (evtType==BdkEvtTypes::SIG_BAD_D)
    dalitzHist = &dalitzHolderN.sigBadD0Type();
  else 
    cout << "Invalid event type."<<endl;

  TString histName = dalitzHist->dataHist().GetName()+TString("_diff");
  
  TFile f(TString("../BToDKTo3piK/params/syst/")+histName+".root");
  hDiff = (TH2*)f.Get(histName);
  hDiff->SetDirectory(gROOT);
  f.Close();
  hDefault = dalitzHist->dataHist();
}


// Fit the data with modified PDF
void bkgHistSystData(const BdkEvtTypes::Type evtType)
{
  loadHistogram(evtType);
  readCut = cutSigReg;
  data = read(dataTree);
 
  changePdf(0, evtType); 
  fit(pdfOnResDK,*data,true);
  RooArgSet* refResult = pdfOnResDK.fitResult();

  changePdf(1, evtType);
  fit(pdfOnResDK,*data,true);
  RooArgSet* modResult = pdfOnResDK.fitResult();

  removeAllRanges(*refResult);
  removeAllRanges(*modResult);

  refResult->Print("v");
  modResult->Print("v");

  sub(*refResult,*modResult);
  refResult->Print("v");
}



// Change PDF to configuration cfg
// cfg = 0 :  default PDF
// type is the event type the difference gets accounted to
void changePdf(int cfg = 0, const BdkEvtTypes::Type evtType)
{
  useYieldFitVars();
  pdfOnResDK.useDalitz();
  readOnResDKPar();
  if (cfg==0) {
    if (hDefault) {
      if (evtType==BdkEvtTypes::DPiX) {
        dalitzHolderN.setDPiXHist(*hDefault);
        dalitzHolderP.setDPiXHist(*hDefault);
      }
      else if (evtType==BdkEvtTypes::SIG_BAD_D) {
        dalitzHolderN.setSigBadD0Hist(*hDefault);
        dalitzHolderP.setSigBadD0Hist(*hDefault);
      }
      else cout << "Invalid event type." << endl;
    }
  }
  else if (cfg==1) {
    if (hMod) delete hMod;
    hMod = (TH2*)hDefault->Clone("hMod");
    hMod->Multiply(hDefault, hDiff);
    if (evtType==BdkEvtTypes::DPiX) {
      dalitzHolderN.setDPiXHist(*hMod);
      dalitzHolderP.setDPiXHist(*hMod);
    }
    else if (evtType==BdkEvtTypes::SIG_BAD_D) {
      dalitzHolderN.setSigBadD0Hist(*hMod);
      dalitzHolderP.setSigBadD0Hist(*hMod);
    }
    else cout << "Invalid event type." << endl;
  }
}



void increaseDiff(TH2& h2, Double_t factor)
{
  for (Int_t xbin=1; xbin<=h2.GetNbinsX(); xbin++) {
    for (Int_t ybin=1; ybin<=h2.GetNbinsY(); ybin++) {
      Double_t bin = h2.GetBinContent(xbin,ybin);
      if (bin!=0) h2.SetBinContent(xbin,ybin,fabs((bin-1)*factor + 1));
    }
  } 
}




void bkgShapeTest()
{
  // Load histogram with MC/data ratio from sidebands
  TString diffFile("hist_dpix_diff.root");
  TFile f(diffFile);
  hDiff = (TH2*)f.Get("hist_dpix_diff");
  hDiff->SetDirectory(gROOT);
  f.Close();
  
  increaseDiff(*hDiff,3);
  TCanvas* can2 = new TCanvas("can2","",400,400);
  hDiff->Draw("colz");

  hDefault = &dalitzHolderN.DPiXType().dataHist();
  

  TCanvas* can = new TCanvas("can","",1200,600);
  can->Divide(4,2);
  can->cd(1);
  changePdf(0);
  dalitzHolderN.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  can->cd(2);
  dalitzHolderP.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  changePdf(1);
  can->cd(3);
  dalitzHolderN.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  can->cd(4);
  dalitzHolderP.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");

  changePdf(0);
  can->cd(5);
  dalitzHolderN.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  can->cd(6);
  dalitzHolderP.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  changePdf(1);
  can->cd(7);
  dalitzHolderN.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");
  can->cd(8);
  dalitzHolderP.DPiXType().generate(20000)->createHistogram(*m12,*m13)->Draw("colz");


  changePdf(0);
}
