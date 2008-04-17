// $Id: fitData.cc,v 1.16 2006/08/04 22:26:34 fwinkl Exp $
// Script to fit the data
// Run in batch with:  
//     bbrroot -q -b setup.cc 'fitData.cc("-Minos")' > & ! fitData-Minos.log &

//#include "../BToDKTo3piK/analysis/mlikecut.cc" 



void fitData(TString postFix = "",
             Bool_t useMinos = kTRUE,           
             Bool_t usePenalty = kTRUE,
             Bool_t doPlot = kFALSE,
             const char* parFile = 0, 
             TTree* tree = dataTree, TCut addCut = "")
{  

  //reload the parameters
  readOnResDKPar();
  if (parFile) pdfOnResDK.parameters().readFromFile(parFile);

  pdfOnResDK.getPdf();

  // read the data
  readCut = TCut(cutSigReg && addCut);

  if (data) delete data;
  data = read(tree);

  // fit
  pdfOnResDK.useMinos(useMinos);
  pdfOnResDK.setMinuitPrintLevel(1);
  //pdfOnResDK.setNLLYieldsSystBit(0);
  fit(pdfOnResDK,*data, usePenalty);

  RooFitResult* yieldFit = pdfOnResDK.yieldFitResult();
  RooFitResult* xyFit = pdfOnResDK.xyFitResult();

  if (yieldFit) printFitResult(yieldFit);
  if (xyFit) printFitResult(xyFit);

  // Save fit result
  TString parResult = TString("fitData")+postFix+".par";
  pdfOnResDK.fitResult()->writeToFile(parResult);
  
  // Append covariance matrix to par file
  ofstream of;
  of.open(parResult, ios_base::out | ios_base::app);
  if (yieldFit) printCovMatrix(yieldFit,of);
  if (xyFit) printCovMatrix(xyFit,of);
  of.close();
  
  // Save fit result in separate file (including LaTeX correlation matrix)
  of.open(TString("fitData")+postFix+".fit");
  if (yieldFit) {
    printCorMatrix(yieldFit,of,true);
    printFitResult(yieldFit,of);
  }
  if (xyFit) {
    printCorMatrix(xyFit,of,true);
    printFitResult(xyFit,of);
  }
  of.close();

  if (doPlot) plotProj(data);

}



// fit different data ranges with different variables
void fitDataRanges()
{
  Bool_t useDA = true;
  Bool_t useQ = true;
  Bool_t useD = true;
  Bool_t useDE = true;

  fitData(useDA,useDE,useQ,!useD,0,"-run15-DA-DE-Q",false,dataTree);
  fitData(!useDA,useDE,useQ,!useD,0,"-run15-DE-Q",false,dataTree);
  fitData(useDA,!useDE,!useQ,!useD,0,"-run15-DA",false,dataTree);
  fitData(!useDA,useQ,useDE,useD,0,"-run15-Q-DE-D",false,dataTree);
  fitData(!useDA,!useQ,!useDE,useD,0,"-run15-D",false,dataTree);

  fitData(useDA,useDE,useQ,!useD,0,"-run14-DA-DE-Q",false,dataTree,cutRun14);
  fitData(!useDA,useDE,useQ,!useD,0,"-run14-DE-Q",false,dataTree,cutRun14);
  fitData(useDA,!useDE,!useQ,!useD,0,"-run14-DA",false,dataTree,cutRun14);
  fitData(!useDA,useQ,useDE,useD,0,"-run14-Q-DE-D",false,dataTree,cutRun14);
  fitData(!useDA,!useQ,!useDE,useD,0,"-run14-D",false,dataTree,cutRun14);

  fitData(useDA,useDE,useQ,!useD,0,"-run5-DA-DE-Q",false,dataTree,cutRun5);
  fitData(!useDA,useDE,useQ,!useD,0,"-run5-DE-Q",false,dataTree,cutRun5);
  fitData(useDA,!useDE,!useQ,!useD,0,"-run5-DA",false,dataTree,cutRun5);
  fitData(!useDA,useQ,useDE,useD,0,"-run5-Q-DE-D",false,dataTree,cutRun5);
  fitData(!useDA,!useQ,!useDE,useD,0,"-run5-D",false,dataTree,cutRun5);
}


void plotNLL() {

  gStyle->SetTitleOffset(1.3,"y");

  RooPlot *p1N, *p2N, *p1P, *p2P;

  if (pdfOnResDK.xMinus()!=0) {
    p1N = pdfOnResDK.xMinus()->frame(-10,10);
    p2N = pdfOnResDK.yMinus()->frame(-10,10);  
    p1P = pdfOnResDK.xPlus()->frame(-10,10);
    p2P = pdfOnResDK.yPlus()->frame(-10,10);  
    p1N->GetXaxis()->SetTitle("x- / x+");
    p2N->GetXaxis()->SetTitle("y- / y+");
  }
  else {
    p1N = pdfOnResDK.rhoMinus()->frame(0.5,2);
    p2N = pdfOnResDK.thetaMinus()->frame(110,220);  
    p1P = pdfOnResDK.rhoPlus()->frame(0.5,2);
    p2P = pdfOnResDK.thetaPlus()->frame(110,220);  
    p1N->GetXaxis()->SetTitle("#rho- / #rho+");
    p2N->GetXaxis()->SetTitle("#theta- / #theta+");
  }

  Bool_t extended = true;
  pdfOnResDK.getPdf()->plotNLLOn(p1N,data,extended);
  pdfOnResDK.getPdf()->plotNLLOn(p2N,data,extended);
  pdfOnResDK.getPdf()->plotNLLOn(p1P,data,extended);
  pdfOnResDK.getPdf()->plotNLLOn(p2P,data,extended);
  p1P->getAttLine()->SetLineStyle(kDashed);
  p2P->getAttLine()->SetLineStyle(kDashed);

  p1P->GetXaxis()->SetTitle("");
  p2P->GetXaxis()->SetTitle("");
  p1P->GetYaxis()->SetTitle("");
  p2P->GetYaxis()->SetTitle("");
 
  p1N->SetTitle("");
  p2N->SetTitle("");
  p1P->SetTitle("");
  p2P->SetTitle("");

  TCanvas* can = new TCanvas("can","dataFit NLL",900,450);
  can->Divide(2,1);
  can->cd(1);
  p1N->Draw();
  p1P->Draw("same");

  can->cd(2);
  p2N->Draw();
  p2P->Draw("same");

  can->SaveAs("fitData-plotNLL.eps");
  can->SaveAs("fitData-plotNLL.root");
}


// Plot PDF and data projections
// Use global dataset "data" by default
void plotProj(RooDataSet* d = 0, Bool_t noTitle = kFALSE)
{
  const Int_t dalitzBINS = 30;

  // Use this times many toy events than data events for projections
  const Int_t scale = 100;

  // Use scaleBins times more bins for PDF than data
  const Int_t scaleBins = 5;

  const Double_t ksLow = 0.489**2;
  const Double_t ksHigh = 0.508**2;

  if (d==0 && data!=0) d = data;
  if (d==0) return;  

  useYieldFitVars();
  pdfOnResDK.useDalitz();


  RooRealVar myqprime(*qprime);
  RooRealVar mydprime(*dprime);
  RooRealVar mym12(*m12);
  RooRealVar mym13(*m13);
  RooRealVar mym23(*m23);
  myqprime.setRange(-5,5);
  mydprime.setRange(-5,5);
  mym12.setBins(dalitzBINS);
  mym13.setBins(dalitzBINS);
  mym23.setBins(dalitzBINS);

  RooArgList vars(*Deltae,myqprime,mydprime,mym12,mym13,mym23);

  Int_t Ntoy = scale*d->numEntries();
  cout << "Generating "<<Ntoy<<" toy events for projections..."<<endl;

  RooDataSet* genData = pdfOnResDK.generate(Ntoy);
  genData->addColumn(*s23);

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleOffset(1.4,"Y");
  gStyle->SetHistMinimumZero();
  TCanvas* can = new TCanvas("can","plotProj",1200,800);
  can->Divide(3,2,0.004,0.004);

  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar& var = (RooRealVar&)vars[i];

    RooBinning dataBins;
    dataBins.setMin(var.getMin());
    dataBins.setMax(var.getMax());
    dataBins.addUniform(var.getBins(),var.getMin(),var.getMax());

    RooBinning pdfBins;
    pdfBins.setMin(var.getMin());
    pdfBins.setMax(var.getMax());
    pdfBins.addUniform(scaleBins*var.getBins(),var.getMin(),var.getMax());

    if (TString(var.GetName())=="m23") {
      // Add bin for m23 veto window
      dataBins.addBoundary(ksLow);
      dataBins.addBoundary(ksHigh);
      pdfBins.addBoundary(ksLow);
      pdfBins.addBoundary(ksHigh);
    }

    TH1* hData = var.createHistogram("hData",Binning(dataBins));
    TH1* hPdf = var.createHistogram("hPdf",Binning(pdfBins));      

    hData->Sumw2();
    hPdf->Sumw2();
    d->fillHistogram(hData,RooArgList(var));
    genData->fillHistogram(hPdf,RooArgList(var));

    hPdf->Scale((double)scaleBins/scale);
  
    if (noTitle) hData->SetTitle("");
    else hData->SetTitle(var.GetTitle()+TString(" projection"));
    hData->SetMarkerStyle(8);
    hData->SetMarkerSize(0.6);

    hPdf->SetMarkerStyle(22);
    hPdf->SetMarkerColor(kBlue);
    hPdf->SetMarkerSize(0.6);
    hPdf->SetLineColor(kBlue);
   
    can->cd(i+1);
    hData->Draw();
    hPdf->Draw("hist same");
    gPad->SaveAs("fitData-"+TString(var.GetName())+".eps");
  }
  
  //  can->SaveAs("fitData-plotProj.eps");
}
