// $Id: diffDataMC.cc,v 1.2 2006/06/14 17:04:39 fwinkl Exp $
// Script to get the difference of MC and data histograms
// Remember to use ntuples that contain sidebands with setupChains(false)


void diffDataMC(TCut cutSB = cutSBmES)
{
  TCut cut = cutSB+cutBasic;

  BdkPdfDalitzHist& dalitzHist = dalitzHolderN.DPiXType();
  //BdkPdfDalitzHist& dalitzHist = dalitzHolderN.sigBadD0Type();

  TH2* hDiff = diffDataMCDalitz(dalitzHist.dataHist(),cut,true);

  TString name = dalitzHist.dataHist().GetName()+TString("_diff");
  hDiff->SetName(name);
  TFile f(name+".root","recreate");
  hDiff->Write();
  f.Close();
  cout << "Histogram written to "+name+".root" << endl;

  hDiff->Draw("colz");
  drawBinning(hDiff);
  
}


// Return relative data/MC difference: bin content = #Data/#MC
// hproto is a protoype 2D-histogram that defines the binning
// set norm to true if you want to normalize the MC to the data
TH2* diffDataMCDalitz(const TH2& hproto, TCut cut, Bool_t norm = kFALSE)
{
  TH2* hMC = hproto.Clone("hMC");
  hMC->Reset();     // clear content, keep binning

  TH2* hData = hMC->Clone("hData");

  hMC->Sumw2();
  hData->Sumw2();

  TString dalitzVars = "d0ppmupmass**2:d0pppupmass**2";
  weightedMCHisto(hMC,dalitzVars,cut,false);
  dataTree->Project(hData->GetName(),dalitzVars,cut);

  // optionally normalize MC to data
  if (norm) hMC->Scale((Double_t)hData->Integral()/hMC->Integral());

  TH2* hDiff = hData->Clone("hDiff");
  hDiff->Sumw2();
  hDiff->Divide(hMC);
  delete hMC;
  delete hData;

  return hDiff;
}
