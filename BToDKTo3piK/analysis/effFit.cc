// $Id: effFit.cc,v 1.12 2006/07/11 21:42:01 fwinkl Exp $
//
// Script to fit the signal efficiency
// It will write the fit result to eff.par in the working directory
// If you want, you can read the start configuration from params/eff.flat
//
// Call with fit=false if you only want to plot
// Call with tree = 0 if you want to fit toyMC
// Supply an optional par file that is used to initialize the PDF before the fit

RooDataSet *effData;
RooDataSet *fitData;

void effFitAll(Bool_t fit=true, TTree *tree=sigFlatTree)
{
  effFit(eff,fit,tree,cutDKGoodD+cutSigReg);
  effFit(effOther,fit,tree,cutSigReg);
}


void effFit(BdkPdfAbsBase& effPdf, Bool_t fit=true, 
            TTree *tree=sigFlatTree, TCut cut = cutDKGoodD+cutSigReg,
            const char* effPar = 0) {

  gStyle->SetPalette(1);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(1.0);
  TGaxis::SetMaxDigits(3);
  gROOT->cd();
  
  if (effPar) effPdf.parameters().readFromFile(effPar);
  
  if (tree) effData = readEffData(tree,cut);
  else effData = effPdf.generate(50000);    // Toy MC
  
  // Fitting
  if (fit) {

    fitOption = "mr";
    optOption = "c2";     // for dual CPU fit
    fit(effPdf, *effData);

    TString parFile = TString(effPdf.GetName())+".par";
    cout << "Writing efficiency parameters to "<<parFile<<endl;
    effPdf.parameters().writeToFile(parFile);

    // Append covariance matrix to par file
    ofstream of;
    of.open(parFile, ios_base::out | ios_base::app);
    printCovMatrix(fitResult,of);
    of.close();
  }
  else effPdf.parameters().Print("v");
  
  // Plotting
  TCanvas *canvas = new TCanvas("Dalitz","",0,0,800,800);
  canvas->Divide(2,2);

  canvas->cd(1);
  TH2 *h2 = effData->createHistogram(*m12,*m13);
  h2->SetTitle("Dalitz plot signal MC");
  h2->GetXaxis()->SetTitle(m12->GetTitle());
  h2->GetYaxis()->SetTitle(m13->GetTitle());
  h2->Draw("colz");
  TGraph *g = ((BdkDalitzBase*)effPdf.getPdf())->drawBoundary(100);
  g->SetLineWidth(2);
  g->Draw("c same");

  canvas->cd(2);
  fitData = (RooDataSet*)effPdf.generate(effData->numEntries());
  TH2 *h2mc = fitData->createHistogram(*m12,*m13);
  h2mc->SetTitle("Dalitz plot generated from PDF");
  h2mc->GetXaxis()->SetTitle(m12->GetTitle());
  h2mc->GetYaxis()->SetTitle(m13->GetTitle());
  h2mc->SetMaximum(h2->GetMaximum());
  h2mc->Draw("colz");
  g->Draw("c same");
  
  if (h2 && h2mc) {
    Double_t KS = h2->KolmogorovTest(h2mc);
    cout << "KS prob. of fit result = "<<KS<<endl;
  }
    
  canvas->cd(3);
  RooPlot *plotm12 = m12->frame();
  effData->plotOn(plotm12,DataError(RooAbsData::SumW2));
  effPdf.getPdf()->plotOn(plotm12);
  plotm12->SetTitle("");
  plotm12->Draw();

  canvas->cd(4);
  RooPlot *plotm13 = m13->frame();
  effData->plotOn(plotm13,DataError(RooAbsData::SumW2));
  effPdf.getPdf()->plotOn(plotm13);
  plotm13->SetTitle("");
  plotm13->Draw();

  cout << "Chisquare m12 "<<plotm12->chiSquare()<<endl;
  cout << "Chisquare m13 "<<plotm13->chiSquare()<<endl;

  TString effName(effPdf.GetName());
  canvas->SaveAs(effName+".root");
  canvas->SaveAs(effName+".eps");
}
  
  
void effToyMC(const char* effParFile = "../BToDKTo3piK/params/eff.par")
{
  RooRealVar s12("s12","",0,3);
  RooRealVar s13("s13","",0,3);

  BdkPdf2DpolyDalitz eff("eff","",s12,s13,BdkDalitzBase::D0);
  eff.parameters().readFromFile(effParFile);
  eff.parameters().Print("v");

  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  cfg->method1D().setLabel("RooSegmentedIntegrator1D");
  eff.getPdf()->setIntegratorConfig(*cfg);

  RooMCStudy toyMC(*eff.getPdf(),*eff.getPdf(),RooArgSet(s12,s13),"","2mc");
  toyMC.generateAndFit(200,40000);

  string filename = "../BToDKTo3piK/toy/effToyMC.root";
  TFile f(filename.c_str(),"recreate");
  RooDataSet toyResult = toyMC.fitParDataSet();
  toyResult.Write();
  f.Close();
  cout << "Toy MC result written to "<<filename.c_str()<<endl;

}


// Incrementally fit the efficiency coefficients
void effFitIter(BdkPdf2DpolyDalitz& effPdf,
                TTree *tree=sigFlatTree, TCut cut = cutDKGoodD+cutSigReg,
                Int_t maxCoeff = 9, const char* effPar = 0)
{
  if (maxCoeff>9 || maxCoeff<1) return;

  if (effPar) effPdf.parameters().readFromFile(effPar);
  effPdf.fixAll();
  effPdf.c(0)->setVal(1);
  for (int i=1;i<=9;i++) effPdf.c(i)->setVal(0);

  effData = readEffData(tree,cut);
  
  for (int i=1;i<=maxCoeff;i++) {
    effPdf.c(i)->setConstant(false);   

    fitOption = "mr";
    optOption = "c2";     // for dual CPU fit
    RooFitResult* result = fit(effPdf, *effData);

    fitData = (RooDataSet*)effPdf.generate(effData->numEntries());
    TH2 *h2 = effData->createHistogram(*m12,*m13);
    TH2 *h2mc = fitData->createHistogram(*m12,*m13);
    Double_t KS = h2->KolmogorovTest(h2mc);

    TString file = TString(effPdf.GetName()) + "-";
    file += i;
    ofstream of;
    of.open(file+".fit");
    result->printToStream(of, RooPrintable::Verbose);
    of << "Fit status = " << result->status() << endl;
    of << "KS prob.   = " << KS << endl;
    of.close();

    effPdf.parameters().writeToFile(file+".par");
    // Append covariance matrix to par file
    ofstream of;
    of.open(file+".par", ios_base::out | ios_base::app);
    printCovMatrix(result,of);
    of.close();

    delete fitData;
  }
  delete effData;
}


RooDataSet* readEffData(TTree *tree, TCut cut)
{
  if (tree==0) return 0;

  RooDataSet* data;
  // Apply PID efficiencies if necessary
  if (tree == sigFlatWTree) {
    RooArgSet vars(*allVars);
    vars.add(*pidEffVars);      
    data = read(tree,0,cut,&vars);
    cout << "Applying data PID efficiencies to MC."<<endl;
    data->addColumn(*eventPidEff);
    data->setWeightVar(*eventPidEff);
  }
  else data = read(tree,0,cut,allVars);

  return data;
}
