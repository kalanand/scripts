// $Id: getDalitzParams.cc,v 1.11 2006/06/26 23:09:16 abi Exp $
// Get the Dalitz shape parameters for different event types


// Gets all the Dalitz shapes
void getDalitzParamsAll(Int_t todo = ALL_TYPES, 
                        Bool_t fit = true, Bool_t plot = true)
{
  TCut cut = cutSigReg;

  // Read the PDF parameters, but don't fix any additional variables
  readOnResDKPar(false);

  if (todo & DPI_B_BIT) getDalitzParams(dpiTree, cutDPiBadD+cut, dalitzHolderN.DpiBadD0(), fit, plot);  
  if (todo & DKX_BIT) getDalitzParams(bbTree, cutDKX+cut, dalitzHolderN.DKX(), fit, plot);
  if (todo & BB_B_BIT) getDalitzParams(bbTree, cutBBBadD+cut, dalitzHolderN.BBBadD0(), fit, plot);
  if (todo & QQ_B_BIT) getDalitzParams(qqTree, cutqqBadD+cut, dalitzHolderN.qqBadD0(), fit, plot);

  if (todo & DPIX_BIT) getDalitzParamsDPiX(bbTree, cutDPiX+cut);
  if (todo & SIG_B_BIT) getDalitzParamsDKBadD(sigTree, cutDKBadD+cut);

  /*
  BdkPdf2DpolyDalitz dpixPdf("dpixPdf","dpixPdf",*m12,*m13);
  BdkPdfDDalitzInc dpixPdf2("dpixPdf2","dpixPdf2",*m12,*m13,BdkDalitzBase::D0,0);
  dpixPdf2.parameters().readFromFile("dpixPdf2.par");
  //  BdkDDalitzAmp::normalizeAll();

  if (type==15) getDalitzParams(bbTree, cutDPiX+cut, &dpixPdf2, fit);

  dpixPdf2.parameters().writeToFile("dpixPdf2.par");
  */
}

// Make RooDataSet from tree
// Add B- and flipped B+ Dalitz plots together
RooDataSet* dalitzDataSet(TTree *tree)
{  
  RooArgSet vars(*mass12,*mass13,*Hdtrkchge);
  gROOT->cd();
  RooDataSet *dataBm = new RooDataSet("dataBm","",tree,vars,"Hdtrkchge<0");
  RooDataSet *dataBp = new RooDataSet("dataBp","",tree,vars,"Hdtrkchge>0");
  
  for (int i=0; i<dataBp->numEntries(); i++) {
    
    RooArgSet *args = dataBp->get(i);
    RooRealVar *arg12 = (RooRealVar*)args->find(mass12->GetName());
    RooRealVar *arg13 = (RooRealVar*)args->find(mass13->GetName());
    
    Double_t temp12 = arg12->getVal();
    arg12->setVal(arg13->getVal());
    arg13->setVal(temp12);

    dataBm->add(*args);
  }
  delete dataBp;
  
  dataBm->addColumn(*s12);
  dataBm->addColumn(*s13);
  
  return dataBm;
}

// Fit "pdf" to "tree" after applying "cut"
void getDalitzParams(TTree *tree, TCut cut, BdkPdfAbsBase *pdf = 0,
                     Bool_t fit = true, Bool_t plot = kTRUE, Bool_t plot3D = kFALSE,
                     Int_t toyEvents = -1)
{

  TString baseFile = TString("getDalitzParams_")+cut.GetName();

  TTree *cutTree;
  RooDataSet *data;
  if (toyEvents<0) {
    // Cut the tree
    gROOT->cd();
    cutTree = tree->CopyTree(cut);  
    // Read tree into RooDataSet
    data = dalitzDataSet(cutTree);
  }
  
  //
  // Fitting
  //
  if (pdf) {
    TString parFile = baseFile + ".par";

    if (toyEvents>0) data = pdf->generate(toyEvents);

    if (fit) {
      fitOption = "rm";
      optOption = "c2";
      Int_t outDalitz = dalitzCfg->inDalitzDataSet(*data,*m12,*m13);
      if (outDalitz>0) cout << outDalitz << " events outside the Dalitz plot."
                            << endl;
      result = fit(*pdf, *data);
      pdf->parameters().writeToFile(parFile);

      // Append covariance matrix to par file
      ofstream of;
      of.open(parFile, ios_base::out | ios_base::app);
      printCovMatrix(result,of);
      of.close();
    }
  }

  if (plot) {    
    //
    // Plotting 
    //
    const Int_t BINS = 100;      // Bins for Dalitz plot

    gStyle->SetOptStat(1111);
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.85);
    gStyle->SetTitleYOffset(1.2);

    // Set the binning for 1D projection plots
    m12->setBins(30);
    m13->setBins(30); 

    TCanvas *canvas = new TCanvas("Dalitz","Dalitz",0,0,1200,300);
    canvas->SetTopMargin(0.05);
    canvas->Divide(4,1,0.001,0.001);

    canvas->cd(1);
    TH2 *h1 = data->createHistogram(*m12,*m13,BINS,BINS);
    h1->SetName(cut.GetName());
    h1->SetTitle(TString("a) Monte Carlo events for ")+cut.GetName());
    h1->GetXaxis()->SetTitle(m12->GetTitle());
    h1->GetYaxis()->SetTitle(m13->GetTitle());
    h1->SetMarkerStyle(21);
    h1->SetMarkerColor(kBlue);
    h1->SetMarkerSize(0.2);
    h1->Draw();

    TGraph *g;         // Dalitz boundary
    if (pdf) {
      // We use the eff PDF here just because it is a globally available
      // PDF derived from BdkDalitzBase that knows about the Dalitz boundaries
      g = ((BdkDalitzBase*)eff.getPdf())->drawBoundary(100);
      g->SetLineWidth(2);
      g->Draw("c same");
    }

     // Reduce the accuracy for 1D projections
    RooNumIntConfig oldIntCfg(*pdf->getPdf()->getIntegratorConfig());  
    //    pdf->getPdf()->setIntegratorConfig(plot1dIntCfg);
    
    canvas->cd(2);
    RooPlot *plotm12 = m12->frame();
    plotm12->SetTitle(TString("b) Projection of ")+m12->GetTitle());
    data->plotOn(plotm12,MarkerSize(0.6),XErrorSize(0));
    if (pdf) {
      pdf->getPdf()->plotOn(plotm12);
      cout << "Chisquare m12 "<<plotm12->chiSquare()<<endl;
    }
    plotm12->Draw();

    canvas->cd(3);
    RooPlot *plotm13 = m13->frame();
    plotm13->SetTitle(TString("c) Projection of ")+m13->GetTitle());
    data->plotOn(plotm13,MarkerSize(0.6),XErrorSize(0));
    if (pdf) {
      pdf->getPdf()->plotOn(plotm13);
      cout << "Chisquare m13 "<<plotm13->chiSquare()<<endl;
    }
    plotm13->Draw();

    // Restore integrator setttings
    pdf->getPdf()->setIntegratorConfig(oldIntCfg);
    
    TH2 *h4;
    if (pdf) {
      canvas->cd(4);
      fitData = (RooDataSet*)pdf->generate(data->numEntries());
      TH2 *h4 = fitData->createHistogram(*m12,*m13,BINS,BINS);
      h4->SetName(TString(h1->GetName())+" (gen)");
      h4->SetTitle("d) Generated from pdf");
      h4->GetXaxis()->SetTitle(m12->GetTitle());
      h4->GetYaxis()->SetTitle(m13->GetTitle());
      h4->SetMarkerStyle(21);
      h4->SetMarkerColor(kBlue);
      h4->SetMarkerSize(0.2);
      h4->Draw();
      g->Draw("c same");
    
      // Create another histogram with much bigger statistics
      TH2* hks = pdf->generate(data->numEntries())
	->createHistogram(*m12,*m13,BINS,BINS);

      Double_t ks = h1->KolmogorovTest(hks);
      cout << "Comparing input histogram with toy generated histogram "
	   << "using 100 times the statistics for "<<cut.GetName()<<":"<<endl;
      cout << "KS probability   = "<<ks<<endl;
      cout << "Chi2 probability = "<<h1->Chi2Test(hks,"")<<endl;

      /*
	canvas->cd(5);
	TH2 *hchi2 = h4->Clone("hchi2");
	TH1D *hpull = new TH1D("hpull","pull",24,-3,3);

	chi2test2d(h1,h4,hchi2,hpull);
	hchi2->Draw("colz");

	canvas->cd(6);
	hpull->Draw();
      */
    }
    
    canvas->SaveAs(baseFile+".root");
    canvas->SaveAs(baseFile+".eps");
  } // if (plot)

  if (plot3D) {
    plotDalitz3D(data, pdf->getPdf());
  }

  delete data;
  delete cutTree;
}


// Define the DKBadD histogram binning
void getDalitzParamsDKBadD(TTree *tree, TCut cut)
{
 // Define binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);

  // We use the same binning for m12 and m13
  Double_t d12[] = {0.15,0.25,0.4,0.5,0.6,0.85,1.15,1.4,1.75,
		    2.15,2.30,2.45,2.60,2.75};
  
  //  Double_t d13[] = {0.24,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
 
  setBinning(bins12,d12,sizeof(d12)/sizeof(Double_t));
  setBinning(bins13,d12,sizeof(d12)/sizeof(Double_t));

  TH2 *hist = getDalitzHist(tree,cut,bins12,bins13);
  hist->SetName("hist_dkbadd");

  // Remove the palette since that gives warnings on reading the root file
  TList *l = hist->GetListOfFunctions();
  l->Remove(l->FindObject("palette"));

  const char* file = "hist_dkbadd.root";
  cout << "Writing DKBadD histogram PDF to "<<file<<endl;
  TFile f(file,"recreate");
  hist->Write();
  f.Close();
}


// Define the DPiX histogram binning
void getDalitzParamsDPiX(TTree *tree, TCut cut)
{
 // Define binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);

  Double_t d12[] = {0.2,0.3,0.4,0.5,0.6,0.8,1.1,1.4,1.7,2.0,
		    2.15,2.30,2.45,2.60,2.75};
  Double_t d13[] = {0.24,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
 
  setBinning(bins12,d12,sizeof(d12)/sizeof(Double_t));
  setBinning(bins13,d13,sizeof(d13)/sizeof(Double_t));

  TH2 *hist = getDalitzHist(tree,cut,bins12,bins13);
  hist->SetName("hist_dpix");

  // Remove the palette since that gives warnings on reading the root file
  TList *l = hist->GetListOfFunctions();
  l->Remove(l->FindObject("palette"));

  const char* file = "hist_dpix.root";
  cout << "Writing DPiX histogram PDF to "<<file<<endl;
  TFile f(file,"recreate");
  hist->Write();
  f.Close();
}


TH2* getDalitzHist(TTree *tree, TCut cut, RooBinning& bins12, RooBinning& bins13)
{
  TString baseFile = TString("getDalitzHist_")+cut.GetName();

  m12->setBinning(bins12);
  m13->setBinning(bins13);

  // Cut the tree
  gROOT->cd();
  TTree *cutTree = tree->CopyTree(cut);
  // Read tree into RooDataSet
  RooDataSet* data = dalitzDataSet(cutTree);

  TH2* h = (TH2*)data->createHistogram("h",*m12,Binning(bins12),
				       YVar(*m13,Binning(bins13)));

  // Weight the bins according to their area inside the Dalitz plot
  ((BdkDalitzBase*)eff.getPdf())->weightBins(h,true,0.0005);
 
  // Plot of original data (at most 5000 data points)
  data->tree().Draw(TString(m13->GetName())+":"+TString(m12->GetName()),"","goff",5000);
  TGraph* g = new TGraph(data->tree().GetSelectedRows(),
                         data->tree().GetV2(), data->tree().GetV1());

  RooDataHist dataHist("dataHist","",RooArgSet(*m12,*m13),h);
  RooHistPdf *histPdf = new RooHistPdf("histPdf","",RooArgSet(*m12,*m13),dataHist);

  bins12.Print();
  bins13.Print(); 

  // Plotting
  gStyle->SetOptStat(0);
  TCanvas *can = new TCanvas("can","can",1200,400);
  can->Divide(3,1);

  can->cd(1);
  h->SetTitle(TString("Events for ")+cut.GetName());
  h->GetYaxis()->SetTitleOffset(1.3);
  h->Draw("colz");
  can->Update();
  TPaletteAxis *pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette"); 
  pal->GetAxis()->SetTitle("");
  pal->GetAxis()->SetTitleSize(0);
  can->Update();


  g->Draw("p same");
  drawBinning(h);
  
  TGraph *g = ((BdkDalitzBase*)eff.getPdf())->drawBoundary();
  g->SetLineWidth(2);
  g->Draw("c same");

  can->cd(2);
  TH1* hx = h->ProjectionX();
  hx->SetTitle("Projection of M(#pi^{+}#pi^{0})^{2}");
  //  hx->GetXaxis()->SetTitle(hfine->GetXaxis()->GetTitle());
  hx->SetLineWidth(2);
  hx->SetLineColor(kBlue);
  hx->Draw();

  can->cd(3);
  TH1* hy = h->ProjectionY();
  hy->SetTitle("Projection of M(#pi^{-}#pi^{0})^{2}");
  hy->GetXaxis()->SetTitleOffset(1);
  //  hy->GetXaxis()->SetTitle(hfine->GetYaxis()->GetTitle());
  hy->SetLineWidth(2);
  hy->SetLineColor(kBlue);
  hy->Draw();

  can->SaveAs(baseFile+".root");
  can->SaveAs(baseFile+".eps");

  return h;
}


void setBinning(RooBinning& b, Double_t bins[], Int_t n)
{
  if (!bins) return;
  for (int i=0; i<n; i++) {
    b.addBoundary(bins[i]);
  }
}


void sigDalitzComp(TTree *tree, TCut cut1, TCut cut2,
                   TTree *effTree = 0) {

  
  gROOT->cd();
  TTree* cutTree1 = tree->CopyTree(cut1);  
  TTree* cutTree2 = tree->CopyTree(cut2);
  RooDataSet* data1 = dalitzDataSet(cutTree1);
  RooDataSet* data2 = dalitzDataSet(cutTree2);
  
  TTree *cutEffTree1 = 0;
  TTree *cutEffTree2 = 0;
  RooDataSet* effData1 = 0;
  RooDataSet* effData2 = 0;
  
  if (effTree) {
    cutEffTree1 = effTree->CopyTree(cut1);
    cutEffTree2 = effTree->CopyTree(cut2);
    effData1 = dalitzDataSet(cutEffTree1);
    effData2 = dalitzDataSet(cutEffTree2);
  }
  

  const Int_t BINS = 50;

  TH2 *h1 = data1->createHistogram(*m12,*m13,BINS,BINS);
  TH2 *h2 = data2->createHistogram(*m12,*m13,BINS,BINS);

  TH2 *e1, *e2;
  if (effTree) {
    e1 = effData1->createHistogram(*m12,*m13,BINS,BINS);
    e2 = effData2->createHistogram(*m12,*m13,BINS,BINS);
    h1->Multiply(e2);
    h2->Multiply(e1);
  }
  
  TCanvas *canvas = new TCanvas("Dalitz","Dalitz",0,0,1200,800);
  canvas->Divide(3,2);

  canvas->cd(1);  
  h1->Draw("colz"); 

  canvas->cd(2);
  h2->Draw("colz");

  canvas->cd(3);
  TH2* hchi2 = h1->Clone("hchi2");
  chi2test2d(h1,h2,hchi2);
  hchi2->Draw("colz");

  canvas->cd(4);
  TH1D *h12_1 = h1->ProjectionX();
  TH1D *h12_2 = h2->ProjectionX();
  h12_1->Sumw2();
  h12_2->Sumw2();
  
  h12_1->Scale(1.0/h12_1->Integral());
  h12_2->Scale(1.0/h12_2->Integral());

  h12_2->SetLineColor(kBlue);
  h12_2->Draw("pe");
  h12_1->Draw("pe same");
  
  canvas->cd(5);  

  TH1D* h13_1 = h1->ProjectionY();
  TH1D* h13_2 = h2->ProjectionY();
  h13_1->Sumw2();
  h13_2->Sumw2();
  
  h13_1->Scale(1.0/h13_1->Integral());
  h13_2->Scale(1.0/h13_2->Integral());

  h13_2->SetLineColor(kBlue);
  h13_2->Draw("pe");
  h13_1->Draw("pe same");
}


void plot3D(TTree *tree, TCut cut, BdkPdfAbsBase *pdf)
{
  gROOT->cd();
  TTree *cutTree = tree->CopyTree(cut);  
  RooDataSet *data = dalitzDataSet(cutTree);
  
  plotDalitz3D(data, pdf->getPdf(),m12,m13,50,50);
}
