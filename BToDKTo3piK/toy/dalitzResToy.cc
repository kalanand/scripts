// $Id: dalitzResToy.cc,v 1.2 2006/03/22 22:48:02 fwinkl Exp $
// Script for toy MC studies of Dalitz resolution

RooAbsPdf* fitPdf;
RooAbsPdf* genPdf;
BdkPdfDalitzSmear* smear;
BdkPdf2DpolyDalitz* dalitzCutOff;

// setup the PDF
void setupPdf(const char* resParFile = 0) {

  dalitzHolderP.sigGoodD0Type().x()->setVal(0);
  dalitzHolderP.sigGoodD0Type().y()->setVal(0);
  
  fitPdf = dalitzHolderP.sigGoodD0Type().getPdf();

  // Dummy resolution models
  /*
  BdkPdfGauss* res12 = new BdkPdfGauss("res12","",*m12);
  BdkPdfGauss* res13 = new BdkPdfGauss("res13","",*m13);
  res12->fixAll();
  res13->fixAll();
  res12.b()->setVal(0);
  res12.s()->setVal(0.1);
  res13.b()->setVal(0);
  res13.s()->setVal(0.1);
  dalitzRes12 = res12;
  dalitzRes13 = res13;
  */
  
  smear = new BdkPdfDalitzSmear("smear","",dalitzHolderP.sigGoodD0Type(),
                                *dalitzRes12, *dalitzRes13);

  if (resParFile) smear.parameters().readFromFile(resParFile);
  smear->setEventBuffer(200);
  
  genPdf = smear->getPdf();
  genPdf->getParameters(RooArgSet()).Print("v");

  // Set correct names for BdkBatchMCStudy
  fitPdf->SetName("*fitPdf");
  genPdf->SetName("*genPdf");
  m12->SetName("*m12");
  m13->SetName("*m13");
}



void batchScripts(int fitsPerJob, const char* toyFile = "dalitzResToy.root",
                  int nSamples = 0, int nEvtPerSample = 0,
                  const char* pathPrefix = "./fitMC",
                  const char* batchSetup = "dalitzResToy-batch.cc")
{ 
  setupPdf();
  BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(*genPdf,*fitPdf,RooArgSet(*m12,*m13),
                                               "","m");

  toyMC->createScripts(toyFile,pathPrefix,batchSetup,fitsPerJob,nSamples,nEvtPerSample);
}
  


void plotResToy(const char* toyFile,
                const char* resToyFile,
                Bool_t fitPull = kFALSE, Bool_t fitErr = kFALSE)
{

  const Double_t pullMin = -4;
  const Double_t pullMax = 4;
  const Int_t pullBins = 40;
  const Double_t errMin = 0.04;
  const Double_t errMax = 0.24;
  const Int_t errBins = 40;
  
  RooRealVar Bpx(*dalitzHolderP.sigGoodD0Type().x());
  RooRealVar Bpy(*dalitzHolderP.sigGoodD0Type().y());
  Bpx.SetTitle("B+ x");
  Bpy.SetTitle("B+ y");
  
  TFile* f1 = new TFile(toyFile);
  TFile* f2 = new TFile(resToyFile);
  RooDataSet *data1 = (RooDataSet*)f1->Get("fitParData");
  RooDataSet *data2 = (RooDataSet*)f2->Get("fitParData");

  // Pulls
  RooPlot* xpull1 = plotToyPull(data1,Bpx.GetName(),Bpx.GetTitle(),
                                fitPull,pullMin,pullMax,pullBins);  
  RooPlot* xpull2 = plotToyPull(data2,Bpx.GetName(),Bpx.GetTitle(),
                                fitPull,pullMin,pullMax,pullBins);
  
  RooPlot* ypull1 = plotToyPull(data1,Bpy.GetName(),Bpy.GetTitle(),
                                fitPull,pullMin,pullMax,pullBins);  
  RooPlot* ypull2 = plotToyPull(data2,Bpy.GetName(),Bpy.GetTitle(),
                                fitPull,pullMin,pullMax,pullBins);  

  // Errors
  RooPlot* xerr1 = plotToyErr(data1,Bpx.GetName(),Bpx.GetTitle(),
                              fitErr,errMin,errMax,errBins);
  RooPlot* xerr2 = plotToyErr(data2,Bpx.GetName(),Bpx.GetTitle(),
                              fitErr,errMin,errMax,errBins);
  RooPlot* yerr1 = plotToyErr(data1,Bpy.GetName(),Bpy.GetTitle(),
                              fitErr,errMin,errMax,errBins);
  RooPlot* yerr2 = plotToyErr(data2,Bpy.GetName(),Bpy.GetTitle(),
                              fitErr,errMin,errMax,errBins);
  
  TCanvas* can = new TCanvas("can","Signal toy MC",800,600);
  gStyle->SetTitleYOffset(1.3);
  can->Divide(2,2);

  can->cd(1);
  xpull1->Draw();
  TAttMarker* marker = 0;
  marker = xpull2->getAttMarker();
  if (marker) {
    marker->SetMarkerColor(kBlue);
    marker->SetMarkerStyle(22);
  }
  xpull2->Draw("same");

  can->cd(2);
  xerr1->Draw();
  marker = xerr2->getAttMarker();
  if (marker) {
    marker->SetMarkerColor(kBlue);
    marker->SetMarkerStyle(22);
  }
  xerr2->Draw("same");
  
  can->cd(3);
  ypull1->Draw();
  marker = ypull2->getAttMarker();
  if (marker) {
    marker->SetMarkerColor(kBlue);
    marker->SetMarkerStyle(22);
  }
  ypull2->Draw("same");

  can->cd(4);
  yerr1->Draw();
  marker = yerr2->getAttMarker();
  if (marker) {
    marker->SetMarkerColor(kBlue);
    marker->SetMarkerStyle(22);
  }

  yerr2->Draw("same");
  
  can->SaveAs("plotResToy.eps");
  can->SaveAs("plotResToy.root");
}

