// $Id: dpiFit.cc,v 1.4 2006/07/11 21:42:00 fwinkl Exp $
// Set of functions to do the Dpi fit validation

BdkPdfGauss* dpiPdfGauss;
BdkPdfPolyn* dpiPdfPolyn;
BdkPdfSum* dpiPdfDeltae;

BdkPdfSum* dpiPdfN;
BdkPdfSum* dpiPdfP;

BdkPdfDDalitzInc* dpiPdfBkgIncN;
BdkPdfDDalitzInc* dpiPdfBkgIncP;

BdkPdfDKDalitz* dpiPdfSigN;
BdkPdfDKDalitz* dpiPdfSigP;

BdkPdfDalitzHist* dpiPdfBkgN;
BdkPdfDalitzHist* dpiPdfBkgP;

RooSimultaneous * simDalitzPdf;

TCut dpiCutSigReg;
TCut dpiCutBasic;
TCut dpiCutDeltae;
TCut dpiCutDeltaeSB;


RooRealVar *dpiDeltae;

void dpiSetup()
{
  gStyle->SetPalette(1);
  gROOT->ProcessLine(".x ../BToDKTo3piK/globals/dpiSetupChains.cc");
  dpiSetupCuts();
  dpiSetupPdf();
}

void dpiSetupCuts(Double_t n2cut = 0.5, Double_t n1cut = 0.4)
{
  dpiCutDeltae = TCut("#DeltaE signal region","(0.005<Deltae&&Deltae<0.090)");
  //  dpiCutDeltae = TCut("#DeltaE signal region","(0.045<Deltae&&Deltae<0.130)");
  dpiCutDeltaeSB = TCut("#DeltaE SB","(-0.08<Deltae&&Deltae<-0.02)||(Deltae>0.1)");

  TString s2;  
  s2.Form("bknnout>%f",n2cut);

  TString s1;  
  s1.Form("nnout>%f",n1cut);

  dpiCutBasic = cutmES+cutMD+cutKsVeto+TCut(s1)+TCut(s2);
  dpiCutSigReg = dpiCutDeltae+dpiCutBasic;

}


void dpiSetupPdf()
{

  dpiDeltae = new RooRealVar("Deltae", "#Delta E", -0.08, 0.14, "GeV");
  setBinning(dpiDeltae, 0.005);
  dpiDeltae->setRange("selection",0.005,0.090);

  dpiPdfGauss = new BdkPdfGauss("dpiPdfGauss","",*dpiDeltae);
  dpiPdfPolyn = new BdkPdfPolyn("dpiPdfPolyn","",*dpiDeltae,1,1);
  dpiPdfPolyn->useOffset(kFALSE);

  dpiPdfDeltae = new BdkPdfSum("dpiPdfDeltae","",*dpiPdfGauss,*dpiPdfPolyn,BdkPdfSum::FULL);

  // Dalitz PDFs
  dpiPdfSigN = &(dalitzHolderN.sigGoodD0Type());
  dpiPdfSigP = &(dalitzHolderP.sigGoodD0Type());
  
  TFile f("hist_dpi.root");
  //TFile f("hist_dpi_dataSB.root");
  TH2* hist = (TH2*)f.Get("hist_dpi");
  hist->SetDirectory(gROOT);
  f.Close();

  dpiPdfBkgN = new BdkPdfDalitzHist("dpiPdfBkgN","",*m12,*m13,BdkDalitzBase::D0,*hist);
  dpiPdfBkgP = new BdkPdfDalitzHist("dpiPdfBkgP","",*m13,*m12,BdkDalitzBase::D0,*hist);

  // Inc PDFs for the background (may not be used. See below):
  dpiPdfBkgIncN = new BdkPdfDDalitzInc("dpiPdfBkgIncN","",*m12,*m13,BdkDalitzBase::D0,0);
  dpiPdfBkgIncP = new BdkPdfDDalitzInc("dpiPdfBkgIncP","",*m12,*m13,BdkDalitzBase::D0BAR,0);

  dpiPdfBkgIncN->setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());
  dpiPdfBkgIncP->setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());

  // Build the total PDF using a sum of the signal plus the hist bgd:
  dpiPdfN = new BdkPdfSum("dpiPdfN","",*dpiPdfSigN,*dpiPdfBkgN,BdkPdfSum::FULL);
  dpiPdfP = new BdkPdfSum("dpiPdfP","",*dpiPdfSigP,*dpiPdfBkgP,BdkPdfSum::FULL);
  
  //  dpiPdfN->parameters().readFromFile("../BToDKTo3piK/params/bdpi/dpiN.par");
  //  dpiPdfP->parameters().readFromFile("../BToDKTo3piK/params/bdpi/dpiP.par");


  /*
  // Use CLEO parameters
  for (int i=0; i<dpiPdfSigN->dalitzAmp()->nComps(); i++)
    dpiPdfSigN->dalitzAmp()->ampRes(i)->setVal(0);
  for (int i=0; i<dpiPdfSigP->dalitzAmp()->nComps(); i++)
    dpiPdfSigP->dalitzAmp()->ampRes(i)->setVal(0);

  dpiPdfN->parameters().readFromFile("dpiN-cleo.par");
  dpiPdfP->parameters().readFromFile("dpiP-cleo.par");
  */
}


void dpiReplaceTrees()
{
  dpiTree = dpiSigTree;
  b0Tree = dpiB0Tree;
  bpTree = dpiBpTree;
  dataTree = dpiDataTree;
  ccTree = dpiCcTree;
  udsTree = dpiUdsTree;
  bbTree = dpiBbTree;
  qqTree = dpiQqTree;
}


void dpiCompareDalitz()
{

  dpiReplaceTrees();
  gROOT->LoadMacro("dataMCSBCompare.cc");
  TCanvas* can1 = new TCanvas("can1","can1",900,300);
  can1->Divide(3,1);
  dataMCSBCompareDalitz(dpiCutDeltaeSB+dpiCutBasic,can1);
  can1->SaveAs("dpiCompareDalitzDataMC.eps");

  const int BINS = 15;

  TH2D *hSB = new TH2D("hSB","MC #DeltaE sideband",BINS,0,3,BINS,0,3);
  TH2D *hSig = new TH2D("hSig","MC signal region",BINS,0,3,BINS,0,3);
  TH2D *hChi2 = new TH2D("hChi2","#sqrt{chi2}",BINS,0,3,BINS,0,3);

  dpiTree = 0; // don't want any signal in MC cocktail
  weightedMCHisto(hSB,"d0ppmupmass**2:d0pppupmass**2",dpiCutDeltaeSB+dpiCutBasic);
  weightedMCHisto(hSig,"d0ppmupmass**2:d0pppupmass**2",dpiCutSigReg);
  dpiTree = dpiSigTree;

  hSB->Scale(1/hSB->Integral());
  hSig->Scale(1/hSig->Integral());
  TVector3 chi2test = chi2test2d(hSB,hSig,hChi2);
  
  TCanvas* can = new TCanvas("can","can",900,300);
  can->Divide(3,1);
  can->cd(1);
  hSB->Draw("colz");

  can->cd(2);
  hSig->Draw("colz");
  
  can->cd(3);
  hChi2->Draw("colz");
  
  can->SaveAs("dpiCompareDalitzMC.eps");

  cout << "Chi2          = " << chi2test[0] << endl;
  cout << "ndof          = " << chi2test[1] << endl;
  cout << "Chi2/ndof     = " << chi2test[0]/chi2test[1] << endl;
  cout << "Chi2 prob     = " << chi2test[2] << endl;
}


// Fit DeltaE for different cuts of N2
void dpiFitDeltaE_N2()
{
  TCanvas* can = new TCanvas("can","",1200,600);
  can->Divide(3,2);

  int i = 1;
  for (double n2 = 0.2; n2<=0.7; n2 += 0.1) {
    dpiSetupCuts(n2);
    
    can->cd(i++);
    TString title;
    title.Form("d > %.1f",n2);

    double sig = dpiFitDeltaE(dpiDataTree, dpiCutBasic, title);
    can->Update();
  }
}


// Fit DeltaE and return S/sqrt(S+B) within "selection" range of dpiDeltae var
double dpiFitDeltaE(TTree* tree = dpiDataTree, 
		    TCut cut = dpiCutBasic, 
		    TString plotTitle = "",
		    double * nsig = 0,
		    double * nsigErr = 0,
		    double * nbkg = 0)
{
  RooArgSet vars(*dpiDeltae);
  RooArgSet showPars(*dpiPdfDeltae->fraction(0),*dpiPdfGauss->b(),*dpiPdfGauss->s());
 
  delete data;
  data = read(tree,0,cut,&vars);
  fitOption = "mr";
  fit(*dpiPdfDeltae,*data);
  
  RooPlot* p = dpiDeltae->frame();    
  data->plotOn(p);
  dpiPdfDeltae->getPdf()->plotOn(p);
  //    dpiPdfDeltae->getPdf()->paramOn(p,Parameters(showPars),Layout(0.1,0.5,0.9));
  p->SetTitle(plotTitle);
  p->Draw();    

  // The total # of signal and background:
  RooAbsReal* fracSig  = dpiPdfGauss->getPdf()->createIntegral(*dpiDeltae,*dpiDeltae,"selection");
  RooAbsReal* fracBkg  = dpiPdfPolyn->getPdf()->createIntegral(*dpiDeltae,*dpiDeltae,"selection");
  
  Double_t nSig = fracSig->getVal()*dpiPdfDeltae->fraction(0)->getVal()*data->numEntries();
  Double_t nSigErr = fracSig->getVal()*dpiPdfDeltae->fraction(0)->getError()*data->numEntries();
  Double_t nBkg = fracBkg->getVal()*(1-dpiPdfDeltae->fraction(0)->getVal())*data->numEntries();

  if (0 != nsig) {
    *nsig = nSig;
  }
  if (0 != nsigErr) {
    *nsigErr = nSigErr;
  }
  if (0 != nbkg) {
    *nbkg = nBkg;
  }
  
  double sig = nSig/sqrt(nSig+nBkg);

  TString s;
  s.Form("S/sqrt(S+B) = %.1f",sig);
  TText t;
  t.DrawTextNDC(0.13,0.83,s);    

  cout << nSig << " "<<nBkg<<" "<< sig <<endl;

  return sig;
}


// Define the histogram binning
void getDalitzParamsDPi(TTree *tree = dpiDataTree,
                        TCut cut = dpiCutBasic+dpiCutDeltaeSB)
{
  gROOT->LoadMacro("../BToDKTo3piK/analysis/getDalitzParams.cc");

 // Define binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);

  // We use the same binning for m12 and m13
  Double_t d12[] = {0.25,0.42,0.62,0.8,1.0,1.2,1.4,1.75,
  		    2.15,2.30,2.5,2.65,2.8};


  //  Double_t d13[] = {0.24,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
 
  setBinning(bins12,d12,sizeof(d12)/sizeof(Double_t));
  setBinning(bins13,d12,sizeof(d12)/sizeof(Double_t));

  TH2 *hist = getDalitzHist(tree,cut,bins12,bins13);
  hist->SetName("hist_dpi");

  gStyle->SetPalette(1);
  // Remove the palette since that gives warnings on reading the root file
  TList *l = hist->GetListOfFunctions();
  l->Remove(l->FindObject("palette"));

  const char* file = "hist_dpi.root";
  cout << "Writing DPI histogram PDF to "<<file<<endl;
  TFile f(file,"recreate");
  hist->Write();
  f.Close();
}


// Define the histogram binning
void getDalitzParamsDPiFromMC(TCut cut = dpiCutSigReg)
{
  dpiReplaceTrees();
  gROOT->LoadMacro("getDalitzParams.cc");

  // Define binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);

  // We use the same binning for m12 and m13
  Double_t d12[] = {0.22,0.3,0.4,0.5,0.7,0.85,1.0,1.15,1.4,1.75,
  		    2.15,2.30,2.45,2.60,2.78};

 
  setBinning(bins12,d12,sizeof(d12)/sizeof(Double_t));
  setBinning(bins13,d12,sizeof(d12)/sizeof(Double_t));

  TH2* histN = (TH2*)m12->createHistogram("histN",Binning(bins12),
                                          YVar(*m13,Binning(bins13)));

  TH2* histP = (TH2*)histN->Clone("histP");
  TH2* hist = (TH2*)histN->Clone("hist_dpi");

  dpiTree = 0;   // don't add signal
  weightedMCHisto(histN,"d0ppmupmass**2:d0pppupmass**2",cut+"Hdtrkchge<0");
  weightedMCHisto(histP,"d0pppupmass**2:d0ppmupmass**2",cut+"Hdtrkchge>0");  
  dpiTree = dpiSigTree;

  hist->Add(histN,histP);

  // Weight the bins according to their area inside the Dalitz plot
  ((BdkDalitzBase*)eff.getPdf())->weightBins(hist,true,0.0005);

  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("can","",500,500); 
  hist->Draw("colz");
  drawBinning(hist);
  TGraph* g = ((BdkDalitzBase*)eff.getPdf())->drawBoundary(100);
  g->SetLineWidth(2);
  g->Draw("c same");

  // Remove the palette since that gives warnings on reading the root file
  TList *l = hist->GetListOfFunctions();
  l->Remove(l->FindObject("palette"));

  const char* file = "hist_dpi.root";
  cout << "Writing DPI histogram PDF to "<<file<<endl;
  TFile f(file,"recreate");
  hist->Write();
  f.Close();
}


void dpiFit(Int_t kcharge = -1, TString suffix = "")
{
  BdkPdfAbsBase* pdf;
  RooArgSet vars(*mass12,*mass13,*Hdtrkchge);
  if (kcharge<0) {
    readCut = dpiCutSigReg + cutMinus;
    pdf = dpiPdfN;
    //    pdf->parameters().readFromFile("dpiN.par");
  }
  else {
    readCut = dpiCutSigReg + cutPlus;
    pdf = dpiPdfP;
    //    pdf->parameters().readFromFile("dpiP.par");
  }

  data = read(dpiDataTree,0,readCut,&vars);
  
  fitOption = "mr";
  fit(*pdf, *data, );

  pdf->parameters().writeToFile("dpiFit"+suffix+".par");
}


// Do yield fits followed by a Dalitz + penalty fit.
// Before calling this, call the following: 
/*
  .L ../BToDKTo3piK/analysis/dpiFit.cc
  dpiSetup();
  getDalitzParamsDPi(dpiDataTree, dpiCutBasic+dpiCutDeltaeSB);
  dpiSetupPdf();
    
*/

void dpiFitPenalty(int todo = 3){ // 1=shape, 2=penalty

  // Calculate the nsig, asym, and bgd fractions:
  double nsigPlus, nsigPlusErr, nbkgPlus;
  double nsigMinus, nsigMinusErr, nbkgMinus;

  TCanvas * can = new TCanvas("DEfitCanvas", "DEfitCanvas", 1000, 500);
  can->Divide(2,1);
  
  can->cd(1);
  dpiFitDeltaE(dpiDataTree, dpiCutBasic + cutMinus, "", &nsigMinus, &nsigMinusErr, &nbkgMinus);
  
  can->cd(2);
  dpiFitDeltaE(dpiDataTree, dpiCutBasic + cutPlus, "", &nsigPlus, &nsigPlusErr, &nbkgPlus);
  
  can->Update();
  can->SaveAs("DpiDEFit.eps");
  can->SaveAs("DpiDEFit.root");

  // set the nsig and asym in the pdfOnResdK:
  double nsig = nsigMinus + nsigPlus;
  double nsigErr = sqrt(nsig);

  double asym = (nsigMinus - nsigPlus) / nsig;
  double asymErr = 1.0 / nsigErr / (1 - asym * asym);

  RooRealVar* nsigVar= (RooRealVar *)(pdfOnResDK.numEvt(BdkEvtTypes::SIG_GOOD_D));
  RooRealVar* asymVar= (RooRealVar *)(pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D));

  nsigVar->setVal(nsig);
  nsigVar->setError(nsigErr);

  asymVar->setVal(asym);
  asymVar->setError(asymErr);

  cout << "Using measured nsig and asym: " << endl;
  nsigVar->Print();
  asymVar->Print();

  // and the fraction of bgd in + and - samples:
  dpiPdfN->fraction(0)->setVal(nbkgMinus / (nbkgMinus + nsigMinus));
  dpiPdfP->fraction(0)->setVal(nbkgPlus / (nbkgPlus + nsigPlus));
  
  cout << "Background fractions:" << endl;
  dpiPdfN->fraction(0)->Print();
  dpiPdfP->fraction(0)->Print();
  
  cout << "Parameters of the PDF:" << endl;
  dpiPdfN->parameters().Print("V");
  dpiPdfP->parameters().Print("V");

  // Create the shape PDF for + and - simultaneously:
  RooArgList pdfList(*(dpiPdfP->getPdf()), *(dpiPdfN->getPdf()));

  simDalitzPdf = new RooSimultaneous("simDalitzPdf", "simDalitzPdf", pdfList, *Hdtrkchge);

  // Read data: 
  delete data;
  readCut = dpiCutSigReg;
  data = read(dpiDataTree);

  // Set B->Dpi efficiency and BR before making the yields NLL:
  pdfOnResDK.absoluteEff()->setVal(0.0813); // Jinlong finds 8.8% before corrections. Then apply the same corrections as in the BAD for B->DK
  pdfOnResDK.absoluteEff()->setError(pdfOnResDK.absoluteEff()->getError() 
				     * 0.0335);

  pdfOnResDK.brBtoDK()->setVal(4.92e-03);
  pdfOnResDK.brBtoDK()->setError(0.21e-03);

  // Make the NLLs:
  RooNLLVar nllShape("nllShape","nllShape", *simDalitzPdf, *data, RooArgSet());
  BdkOnResNLLYields nllYields(pdfOnResDK, 0); 
  nllYields.Print("V");
  
  RooFormulaVar nllBoth("nllBoth", "nllBoth", "@0+@1",
			RooArgSet(nllShape, nllYields));

  RooAbsReal * nllFinal = 0;
  if (1 == todo) {
    nllFinal = &nllShape;
  }
  else if (2 == todo) {
    nllFinal = nllYields;
  }
  else if (3 == todo){
    nllFinal = &nllBoth;
  }

  // Fix and float:
  dpiPdfN->fixAll();
  dpiPdfP->fixAll();
  pdfOnResDK.fixAll();

  fixAll(pdfOnResDK.cpParams(),false);

  // Fit:
  if (minuit) delete minuit;
  minuit = new RooMinuit(*nllFinal) ;
  minuit->optimizeConst(1) ;
  minuit->fit(fitOption) ;
}



void dpiFit_DeltaE()
{
  double dEwindow = 0.040;
  double dElow = -0.040;
  double dEhigh = 0.080;

  TCanvas* can = new TCanvas("can","can",1200,800);
  can->Divide(3,2);
  int i = 1;
  for (double dE = dElow; dE<=dEhigh; dE += dEwindow/2) {
    dpiDeltae->setRange("selection",dE,dE+dEwindow);

    TString title;
    title.Form("%.3f_%.3f",dE,dE+dEwindow);
    can->cd(i++);
    double sig = dpiFitDeltaE(dpiDataTree, dpiCutBasic, title);
    
    TString cut;
    cut.Form("%.3f<Deltae&&Deltae<%.3f",dE,dE+dEwindow);
    TCut dEcut("#DeltaE signal region",cut);
    dpiCutSigReg = dEcut+dpiCutBasic;

    dpiFit(1,"_DeltaE_"+title+"P");
    dpiFit(-1,"_DeltaE_"+title+"N");
  }
  can->SaveAs("dpiFit_DeltaE.eps");
}
