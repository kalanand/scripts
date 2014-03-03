// $Id: signalToy.cc,v 1.8 2006/05/03 17:36:42 fwinkl Exp $
// Toy Monte Carlo of Dalitz signal PDFs

#include "signalToy.hh"

// setup the RooSimultaneous PDF
void setupSimPdf() {
  pdfBp = (BdkPdfDKDalitz*)dalitzHolderP.sigGoodD0();
  pdfBm = (BdkPdfDKDalitz*)dalitzHolderN.sigGoodD0();

  /*
  ((BdkDKDalitz*)pdfBp->getPdf())->calcNormParams(1e7);
  ((BdkDKDalitz*)pdfBp->getPdf())->printNormParams();
  ((BdkDKDalitz*)pdfBm->getPdf())->calcNormParams(1e7);
  ((BdkDKDalitz*)pdfBm->getPdf())->printNormParams();
  */

  simPdf.addPdf(*pdfBp->getPdf(),"+");
  simPdf.addPdf(*pdfBm->getPdf(),"-");

  // Remove efficiency function and load correct normalization parameters

  //  BdkDalitzBase::setEfficiencyFunc(0);
  //  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/dalitzNorm-noeff.par");

  //BdkDDalitzAmp::normalizeAll();

  
  /*
  pdfBm->x()->setVal(0.077);
  pdfBm->y()->setVal(0.064);
  pdfBp->x()->setVal(-0.129);
  pdfBp->y()->setVal(0.019);
  */
  /*
  pdfBm->x()->setVal(0);
  pdfBm->y()->setVal(0);
  pdfBp->x()->setVal(0);
  pdfBp->y()->setVal(0);
  */

  //  ((Bdk2DpolyDalitz*)pdfBp.getPdf())->setEpsRel(1e-1);
  //  ((Bdk2DpolyDalitz*)pdfBm.getPdf())->setEpsRel(1e-1);
  simPdf.getParameters(RooArgList())->Print("v");

  // This is a global prototype dataset with 100k alternating kcharge entries
  TFile f("awg/toyMC/proto.root");
  proto = new RooDataSet(*(RooDataSet*)f.Get("proto"));
  f.Close();
   
  m12->SetName("*m12");
  m13->SetName("*m13");
}

// create prototype dataset for K charge
RooDataSet* createProto(int Nevents) {
  // Create the prototype dataset with the RooCategory
  RooDataSet* proto = new RooDataSet("proto","",RooArgSet(kcharge));
  
  for (int i=0;i<2*Nevents;i++) {
    if (i%2==0) kcharge.setLabel("+");
    else kcharge.setLabel("-"); 
    proto->add(kcharge);
  }
  //  kcharge.setLabel("-");
  //  for (int i=0;i<Nevents;i++) proto->add(kcharge);

  return proto;
}

// Generate the toy MC samples
// "events" events for each B+ and B-
void generate(int samples, int events, const char* toyFile = "signalToy.root")
{
  setupSimPdf(); 
  BdkBatchMCStudy toyMC(simPdf,simPdf,RooArgSet(*m12,*m13),"","",
                        proto,RooArgSet(kcharge));

  toyMC.generateBatch(toyFile,samples,2*events);
}


// Create the batch scripts to do the toy fits
void batchScripts(int fitsPerJob, const char* toyFile = "signalToy.root",
                  int nSamples = 0, int nEvtPerSample = 0,
                  const char* pathPrefix = "./fitMC",
                  const char* batchSetup = "signalToy-batch.cc")
{ 
  setupSimPdf();
  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(simPdf,simPdf,RooArgSet(*m12,*m13),
                                               "","",proto,RooArgSet(kcharge));

  toyMC->createScripts(toyFile,pathPrefix,batchSetup,fitsPerJob,nSamples,nEvtPerSample);
}


// generate, fit and plot one toy MC sample
void fitOneToy(int Nevents = 1000, Bool_t plot = true)
{
  gStyle->SetPalette(1);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetPadLeftMargin(0.12);

  setupSimPdf();
  RooDataSet* proto = createProto(Nevents);
  RooDataSet *data =  simPdf.generate(RooArgSet(*m12,*m13),*proto);

  RooFitResult *result = simPdf.fitTo(*data,NumCPU(2),Save(true),Minos(true),Hesse(false));
  result->Print("v");

 // Plotting
  if (plot) {
    // Reduce precision for plotting
    RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
    cfg->setEpsAbs(1E-5);
    cfg->setEpsRel(1E-5);
    cfg->method1D().setLabel("RooSegmentedIntegrator1D");
    simPdf.setIntegratorConfig(*cfg);
    
    TCanvas *canvas = new TCanvas("Dalitz","",0,0,1050,350);
    canvas->Divide(3,1);
    
    canvas->cd(1);
    TH2 *h2 = data->createHistogram(*m12,*m13);
    h2->SetTitle("Dalitz plot signal MC");
    h2->GetXaxis()->SetTitle(m12->GetTitle());
    h2->GetYaxis()->SetTitle(m13->GetTitle());
    h2->Draw();
    TGraph *g = ((BdkDalitzBase*)pdfBp->getPdf())->drawBoundary();
    g->SetLineWidth(2);
    g->Draw("c same");
    
    canvas->cd(2);
    RooPlot *plotm12 = m12->frame();
    data->plotOn(plotm12);
    simPdf.plotOn(plotm12,ProjWData(kcharge,*data));
    plotm12->Draw();
    
    canvas->cd(3);
    RooPlot *plotm13 = m13->frame();
    data->plotOn(plotm13);
    simPdf.plotOn(plotm13,ProjWData(kcharge,*data));
    plotm13->Draw();
  }
}


// make plots of toy MC fit results
// "files" contains a list of files that were created with fitToyMC
void plotToyResult2(vector<string> files,
                    RooArgList &params)
{
  RooDataset data("data","",params);

  
  for (int i=0; i<files.size(); i++) {
    TFile f(files[i].c_str());
    TList *list = f.GetListOfKeys();
    TIter nextResult(list);
    TKey *key = 0;
    while (key = (TKey*)nextResult()) {
      RooFitResult *fitResult = (RooFitResult*)key->ReadObj();

      // Loop over all parameters given to function
      TIterator *nextParam = params.createIterator();
      RooRealVar *param = 0;
        RooRealVar *r = (RooRealVar*)fitResult->floatParsFinal()
        if (r) {
          r->Print();
        }
        else {
          cout << "Parameter "<<param.GetName()<<" not found."<<endl;
          return;
        }
      }
  }
  f.Close();
  
}

void plotSignalToy(const char* toyFile,
                   Bool_t fitPull = kTRUE, Bool_t fitErr = kFALSE)
{
  TFile* f = new TFile(toyFile);
  RooDataSet* data = (RooDataSet*)f->Get("fitParData");

  TCanvas* can = new TCanvas("can","Signal toy MC",1000,500);
  can->Divide(4,2);

  RooRealVar Bpx(*dalitzHolderP.sigGoodD0Type().x());
  RooRealVar Bpy(*dalitzHolderP.sigGoodD0Type().y());
  RooRealVar Bmx(*dalitzHolderN.sigGoodD0Type().x());
  RooRealVar Bmy(*dalitzHolderN.sigGoodD0Type().y());
  
  Bpx.SetTitle("B+ x");
  Bpy.SetTitle("B+ y");
  Bmx.SetTitle("B- x");
  Bmy.SetTitle("B- y");
  
  plotToy(data, Bpx, fitPull, fitErr, 0, 0.5, 40, -4, 4, 40, can, 1);
  plotToy(data, Bpy, fitPull, fitErr, 0, 0.5, 40, -4, 4, 40, can, 3);
  plotToy(data, Bmx, fitPull, fitErr, 0, 0.5, 40, -4, 4, 40, can, 5);
  plotToy(data, Bmy, fitPull, fitErr, 0, 0.5, 40, -4, 4, 40, can, 7);

  can->SaveAs("plotSignalToy.eps");
  can->SaveAs("plotSignalToy.root");
}

