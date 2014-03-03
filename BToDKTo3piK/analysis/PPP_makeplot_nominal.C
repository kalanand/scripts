{
  gROOT->Reset();
  gStyle->Reset();
  gROOT->SetStyle("BABAR");

  bool doNorm = true;
  int numBINS = 30;


  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0}) [GeV^{2}/c^{4}]",0,3); //This is pi+ pi0  
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]",0,3); //This is pi- pi0  
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+}) [GeV^{2}/c^{4}]",0,3); //This is pi+ pi-
  fit_S23.setBins(numBINS);
  fit_S31.setBins(numBINS);
  fit_S12.setBins(numBINS);


  //read in the background dataset
  TFile hello("hhPi0.root","READ");
  TTree* suptree = (TTree*)hello.Get("pipipi0_data");
  gROOT->cd();


  TCut remKs("(sqrt(fit_S12)<0.489||sqrt(fit_S12)>0.508)&& abs(hh_SignedDistance)<2.0");
//   TTree* sigtree = suptree->CopyTree(remKs && "abs(Dmass-1.8637)<0.016");
  TTree* sigtree = suptree->CopyTree(remKs && "abs(Dmass-1.864)<0.02");
  TTree* bkgtree = suptree->CopyTree(remKs && "Dmass>1.93 && Dmass<1.99");


  RooArgList myList;
  myList.add(RooArgSet(fit_S23,fit_S31,fit_S12));
  RooDataSet *data = new RooDataSet("data","data",sigtree,myList);
  RooDataSet *d2 =  new RooDataSet("d2","d2",bkgtree,myList); 


  BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
  dalitzCfg->getParameters(RooArgSet())->readFromFile("../BToDKTo3piK/Dstar/dalitzCfg.par");
  dalitzCfg->getParameters(RooArgSet())->Print("v");


  // define signal pdf
  BdkPdf2DpolyDalitz eff("eff","",fit_S23, fit_S31);
  eff.parameters().readFromFile("PPP_eff.par");
  BdkPdfDDalitz dstar("dstar","dstar",fit_S23,fit_S31,BdkDalitzBase::D0);
  dstar.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());


  // define bkg pdf
  TH2D bgHist("bgHist","2d bkgd hist from data sideband",6,0,3,6,0,3);
  bkgtree->Draw("fit_S23:fit_S31>>bgHist","","goff");
  ((BdkDalitzBase*) dstar.getPdf())->weightBins(&bgHist,true,0.001);
  BdkDalitzHist bkgpdf("bkgpdf", "bkgpdf", BdkDalitzBase::D0, 
		       BdkDalitzBase::PPP0, fit_S23, fit_S31, bgHist);



  // define signal+bkg pdf
  RooRealVar fsig("fsig","fsig",0.982);
  RooRealVar fbkg("fbkg","fbkg",1.0-fsig.getVal());
  RooAddPdf mypdf("mypdf","mypdf",RooArgList(*dstar.getPdf(),bkgpdf),
		  RooArgList(fsig,fbkg));


  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(5E-7);
  cfg->setEpsRel(5E-7);
 mypdf.setIntegratorConfig(*cfg);


  RooArgSet allPars = dstar.parameters();
  dstar.parameters().readFromFile("PPP_NominalFit.par");



  // Machinary to plot the third variable
  Double_t md0 = 1.8645;
  Double_t mpi0=0.1349766;
  Double_t mpi = 0.13957;
  Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;
  RooRealVar totalm("totalm","totalm",total);
  RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));


  BdkPdfDDalitz dstar12("dstar12","dstar12",fit_S23,fit_S31a,BdkDalitzBase::D0); 
  dstar12.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());
  RooNumIntConfig* cfg2 = RooAbsReal::defaultIntegratorConfig();
  cfg2->setEpsAbs(2E-6);
  cfg2->setEpsRel(2E-6);
 dstar12.getPdf()->setIntegratorConfig(*cfg2);
  RooArgSet allPars12 = dstar12.parameters();
  TIterator * iter12 = allPars12.createIterator();
  RooRealVar * var;
  while(0!= (var= (RooRealVar*)iter12->Next())) {
    TString vname = TString(var->GetName()).ReplaceAll("dstar12","dstar");
    if(vname.Contains("dstar.")) {
      RooRealVar* arvar = allPars.find(vname);
      var->setVal(arvar->getVal());
    }
  }
  
  delete iter12;



  RooAddPdf mypdf12("mypdf12","mypdf12",RooArgList(*dstar12.getPdf(),bkgpdf),
		    RooArgList(fsig,fbkg));



RooDataSet *fitdata = mypdf.generate(RooArgSet(fit_S23,fit_S31),data->numEntries());
RooFormulaVar fit_S12a("fit_S12","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S31));
fitdata->addColumn(fit_S12a);

  //Make the plots
  //******************************************************
//   RooPlot* xframe = fit_S23.frame(0,3,30);
//   data->plotOn(xframe,MarkerSize(0.1),DrawOption("z"));
//   mypdf.plotOn(xframe);
//   xframe->getAttLine()->SetLineWidth(1);
//   xframe->getAttLine()->SetLineStyle(1);
//   xframe->SetTitleSize(0.08, "Y" );
//   xframe->SetTitleOffset(1.6, "Y");
//   xframe->SetNdivisions(505, "X");
//   xframe->SetNdivisions(503, "Y");
//   xframe->GetYaxis()->SetTitle("Events / 0.1 GeV^{2}/c^{4}");
//   xframe->SetTitle("");

  RooPlot* yframe = fit_S31.frame(0,3,30);
  fitdata->plotOn(yframe,MarkerSize(0.1),DrawOption("z"));
  mypdf.plotOn(yframe); 
  yframe->getAttLine()->SetLineWidth(1);
  yframe->getAttLine()->SetLineStyle(1);
  yframe->SetTitleSize(0.08, "Y" );
  yframe->SetTitleOffset(1.6, "Y");
  yframe->SetNdivisions(505, "X");
  yframe->SetNdivisions(503, "Y");
  yframe->GetYaxis()->SetTitle("Events / 0.1 GeV^{2}/c^{4}");
  yframe->SetTitle("");


  RooBinning mbBinning(30,0,3);
//   mbBinning.addUniform(1,0.24,0.26);
  mbBinning.addUniform(15,0.1,0.4);
  fit_S12.setBinning(mbBinning);
  RooPlot* zframe = fit_S12.frame();
  fitdata->plotOn(zframe,MarkerSize(0.1),DrawOption("z"));
  mypdf12.plotOn(zframe);
  zframe->getAttLine()->SetLineWidth(1);
  zframe->getAttLine()->SetLineStyle(1); 
  zframe->SetTitleSize(0.03, "Y" );
  zframe->SetTitleOffset(1.9, "Y");
  zframe->GetYaxis()->SetTitle("Events / 0.1 GeV/c^{2}");
  zframe->SetTitle("");




//   TCanvas *canvas1 = new TCanvas("canvas1","allevents",400,400);
//   xframe->Draw();
  TCanvas *canvas2 = new TCanvas("canvas2","allevents",400,400);
  yframe->Draw();

  TCanvas *canvas3 = new TCanvas("canvas3","allevents",400,400);
  zframe->Draw();




/*

TH1D  histDat1("histDat1","",60,0,3);
TH1D  histDat2("histDat2","",60,0,3);
histDat1.Sumw2();
histDat2.Sumw2();
histDat1.GetXaxis()->SetTitle("m^{2}_{#pi^{+}#pi^{0}} [GeV^{2}/c^{4}]");
histDat1.GetYaxis()->SetTitle("Events/ 0.05 GeV^{2}/c^{4}");

histDat2.GetXaxis()->SetTitle("m^{2}_{#pi^{-}#pi^{0}} [GeV^{2}/c^{4}]");
histDat2.GetYaxis()->SetTitle("Events/ 0.05 GeV^{2}/c^{4}");


TH1D  histFit1("histFit1","",60,0,3);
TH1D  histFit2("histFit2","",60,0,3);
histFit1.Sumw2();
histFit2.Sumw2();
histFit1.SetLineColor(4);
histFit2.SetLineColor(4);

sigtree->Draw("fit_S23>>histDat1","","goff");
sigtree->Draw("fit_S31>>histDat2","","goff");
fitdata->tree().Draw("fit_S23>>histFit1","","goff");
fitdata->tree().>Draw("fit_S31>>histFit2","","goff");

gStyle->SetOptStat(0);
TCanvas *can1 = new TCanvas("can1","",400,400);
histDat1.Draw();
histFit1.Draw("HIST same");
TCanvas *can2 = new TCanvas("can2","",400,400);
histDat2.Draw();
histFit2.Draw("HIST same");

*/

//   TCanvas *canvas = new TCanvas("canvas","allevents",900,300);
//   canvas->Divide(3,1);
//   canvas->cd(1);xframe->Draw();
//   canvas->cd(2);yframe->Draw();
//   canvas->cd(3);
//   zframe->Draw();
//   canvas->SaveAs(TString("ppp_nominal.gif"));
//   canvas->SaveAs(TString("ppp_nominal.eps"));
//   canvas->SaveAs(TString("ppp_nominal.C"));
//  delete canvas;


  //
  // Plot chi^2
  //

//   TH2D h("h","", numBINS,0,3,numBINS,0,3);
//   sigtree->Draw("fit_S23:fit_S31>>h","","goff");
//   h->GetXaxis()->SetTitle("m^{2}_{#pi^{-}#pi^{0}}");
//   h->GetYaxis()->SetTitle("m^{2}_{#pi^{+}#pi^{0}}");
//   h->SetTitle("Data");
//   RooDataSet *fitdata = mypdf.generate(RooArgSet(fit_S23,fit_S31),data->numEntries());
//   RooFormulaVar fit_S12a("fit_S12","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S31));
//   fitdata->addColumn(fit_S12a);

//   TH2F* hFit = fitdata->createHistogram(fit_S31,fit_S23);
//   h.Sumw2();
//   hFit->Sumw2();
//   Double_t hintegral = h.Integral();
//   Double_t hFitintegral = hFit->Integral();
//   Double_t Normalization = hintegral/hFitintegral;
//   hFit->Scale(Normalization);
//   hFit->GetXaxis()->SetTitle("m^{2}_{#pi^{-}#pi^{0}}");
//   hFit->GetYaxis()->SetTitle("m^{2}_{#pi^{+}#pi^{0}}");
//   hFit->SetTitle("Data");

//   TH2D *hchi2 = new TH2D("hchi2","",numBINS,0,3,numBINS,0,3);  
//   TH1D* pullhist = new TH1D("pullhist","",25,-4,4); 
//   double chi2Total =0.0;
//   int ndof = numBINS*numBINS;
//   ndof--;   // histos are normalized to each other

//   Double_t chi2Total = 0;    // total chi2
//   for (int x=1; x <= numBINS; x++) {
//     for (int y=1; y <= numBINS; y++) {
      
//       Double_t bin1 =  h.GetBinContent(x,y);
//       Double_t bin2 = hFit->GetBinContent(x,y); 

//       if (bin2==0 && bin1==0) ndof--;    // no data -> one less dof
//       else {
// 	Double_t sqrtChi2 = (bin1-bin2);
// 	sqrtChi2 /= sqrt(h.GetBinError(x,y)**2 + hFit->GetBinError(x,y)**2);
// 	hchi2->SetBinContent(x,y,sqrtChi2);
// 	chi2Total += sqrtChi2**2;
// 	pullhist->Fill(sqrtChi2);
//       }
//     }
//   }


//   Double_t chi2Prob = TMath::Prob(chi2Total, ndof);
//   cout << "Chi2          = " << chi2Total << endl;
//   cout << "ndof          = " << ndof << endl;
//   cout << "Chi2/ndof     = " << chi2Total/ndof << endl;
//   cout << "Chi2 prob     = " << chi2Prob << endl;


//   TCanvas* aCanvas = new TCanvas("c3", "c3", 900, 600);
//   gStyle->SetOptStat(0);
//   aCanvas->Divide(2,2);
//   aCanvas->cd(1);
//   h.Draw("colz");
//   aCanvas->cd(2);
//   hFit->Draw("colz");
//   aCanvas->cd(3);
//   hchi2->SetMaximum(4);
//   hchi2->SetMinimum(-4);
//   hchi2->SetTitle("#sqrt{#chi^{2}}");
//   hchi2->GetXaxis()->SetTitle("m^{2}_{#pi^{-}#pi^{0}}");
//   hchi2->GetYaxis()->SetTitle("m^{2}_{#pi^{+}#pi^{0}}");
//   hchi2->Draw("colz");
//   aCanvas->cd(4);  
//   pullhist->SetTitle("#sqrt{#chi^{2}} 1D");
//   pullhist->Fit("gaus");
//   gStyle->SetStatX(0.99);
//   gStyle->SetStatY(0.9);
//   gStyle->SetOptFit(11);
//   pullhist->Draw();
//   aCanvas->Draw();
//   aCanvas->SaveAs(TString("ppp_chi_nominal.gif"));
//   aCanvas->SaveAs(TString("ppp_chi_nominal.eps"));
//   aCanvas->SaveAs(TString("ppp_chi_nominal.C"));
}  //end the macro







