// Does 2^6 Dalitz plot fit for data. Each fit corresponds to a different
// starting point. Each starting point corresponds to either a flip of sign 
// or an addition of 180 in one of the phases. 

void DP_PhaseScan() {
  TFile* resultFile = new TFile("./PhaseScan/DP_PhaseScan.root","update");

  Loop(*resultFile, 0, 0, 0, 0, 0, 0);
  Loop(*resultFile, 1, 1, 1, 1, 1, 1);

  Loop(*resultFile, 1, 0, 0, 0, 0, 0);
  Loop(*resultFile, 0, 1, 0, 0, 0, 0);
  Loop(*resultFile, 0, 0, 1, 0, 0, 0);
  Loop(*resultFile, 0, 0, 0, 1, 0, 0);
  Loop(*resultFile, 0, 0, 0, 0, 1, 0);
  Loop(*resultFile, 0, 0, 0, 0, 0, 1);

  Loop(*resultFile, 1, 1, 0, 0, 0, 0);
  Loop(*resultFile, 0, 1, 1, 0, 0, 0);
  Loop(*resultFile, 0, 0, 1, 1, 0, 0);
  Loop(*resultFile, 0, 0, 0, 1, 1, 0);
  Loop(*resultFile, 0, 0, 0, 0, 1, 1);
  Loop(*resultFile, 1, 0, 0, 0, 0, 1);
  Loop(*resultFile, 1, 0, 1, 0, 0, 0);
  Loop(*resultFile, 0, 1, 0, 1, 0, 0);
  Loop(*resultFile, 0, 0, 1, 0, 1, 0);
  Loop(*resultFile, 0, 0, 0, 1, 0, 1);
  Loop(*resultFile, 1, 0, 0, 0, 1, 0);
  Loop(*resultFile, 0, 1, 0, 0, 0, 1);
  Loop(*resultFile, 1, 0, 0, 1, 0, 0);
  Loop(*resultFile, 0, 1, 0, 0, 1, 0);
  Loop(*resultFile, 0, 0, 1, 0, 0, 1);

  Loop(*resultFile, 1, 1, 1, 0, 0, 0);
  Loop(*resultFile, 0, 1, 1, 1, 0, 0);
  Loop(*resultFile, 0, 0, 1, 1, 1, 0);
  Loop(*resultFile, 0, 0, 0, 1, 1, 1);
  Loop(*resultFile, 1, 0, 0, 0, 1, 1);
  Loop(*resultFile, 1, 1, 0, 0, 0, 1);
  Loop(*resultFile, 1, 0, 1, 1, 0, 0);
  Loop(*resultFile, 0, 1, 0, 1, 1, 0);
  Loop(*resultFile, 0, 0, 1, 0, 1, 1);
  Loop(*resultFile, 1, 0, 0, 1, 0, 1);
  Loop(*resultFile, 1, 1, 0, 0, 1, 0);
  Loop(*resultFile, 0, 1, 1, 0, 0, 1);
  Loop(*resultFile, 1, 0, 0, 1, 1, 0);
  Loop(*resultFile, 0, 1, 0, 0, 1, 1);
  Loop(*resultFile, 1, 0, 1, 0, 0, 1);
  Loop(*resultFile, 1, 1, 0, 1, 0, 0);
  Loop(*resultFile, 0, 1, 1, 0, 1, 0);
  Loop(*resultFile, 0, 0, 1, 1, 0, 1);
  Loop(*resultFile, 0, 0, 0, 1, 1, 1);
  Loop(*resultFile, 1, 1, 0, 0, 0, 1);

  Loop(*resultFile, 1, 1, 1, 1, 0, 0);
  Loop(*resultFile, 0, 1, 1, 1, 1, 0);
  Loop(*resultFile, 0, 0, 1, 1, 1, 1);
  Loop(*resultFile, 1, 0, 0, 1, 1, 1);
  Loop(*resultFile, 1, 1, 0, 0, 1, 1);
  Loop(*resultFile, 1, 1, 1, 0, 0, 1);
  Loop(*resultFile, 1, 0, 1, 1, 1, 0);
  Loop(*resultFile, 0, 1, 0, 1, 1, 1);
  Loop(*resultFile, 1, 0, 1, 0, 1, 1);
  Loop(*resultFile, 1, 1, 0, 1, 0, 1);
  Loop(*resultFile, 1, 1, 1, 0, 1, 0);
  Loop(*resultFile, 0, 1, 1, 1, 0, 1);
  Loop(*resultFile, 1, 0, 1, 1, 0, 1);
  Loop(*resultFile, 1, 1, 0, 1, 1, 0);
  Loop(*resultFile, 0, 1, 1, 0, 1, 1);

  Loop(*resultFile, 0, 1, 1, 1, 1, 1);
  Loop(*resultFile, 1, 0, 1, 1, 1, 1);
  Loop(*resultFile, 1, 1, 0, 1, 1, 1);
  Loop(*resultFile, 1, 1, 1, 0, 1, 1);
  Loop(*resultFile, 1, 1, 1, 1, 0, 1);
  Loop(*resultFile, 1, 1, 1, 1, 1, 0);


  delete resultFile;
} 


void Loop(TFile& resultFile, bool NR_180, bool NR_Sign, bool RhoM_180, 
	  bool RhoM_Sign, bool Rho0_180, bool Rho0_Sign){
  gROOT->Reset();
  bool doNorm = false;

  //calculate the dalitz boundary
  Double_t md0 = 1.8645;
  Double_t mpi0=0.1349766;
  Double_t mpi = 0.13957;

  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",0,3); 
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",0,3); 
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",0,3); 
  fit_S23.setBins(60);
  fit_S31.setBins(60);
  fit_S12.setBins(60);


  //read in the background dataset
  TFile hello("hhPi0.root","READ");
  TTree* suptree = (TTree*)hello.Get("pipipi0_data");
  gROOT->cd();
  TCut remKs("(hh_mass<0.4899||hh_mass>0.5079) && abs(hh_SignedDistance)<2.0");
  TTree* sigtree = suptree->CopyTree(remKs && "abs(Dmass-1.864)<0.016");
  TTree* bkgtree = suptree->CopyTree(remKs && "Dmass>1.93");
  RooArgList myList;
  myList.add(RooArgSet(fit_S23,fit_S31,fit_S12));
  RooDataSet *data = new RooDataSet("data","data",sigtree,myList);
  RooDataSet *d2 =  new RooDataSet("d2","d2",bkgtree,myList); 


  // define signal pdf
  BdkPdf2DpolyDalitz eff("eff","",fit_S31, fit_S23);
  eff.parameters().readFromFile("eff.par");
  BdkPdfDDalitz sigpdf("sigpdf","sigpdf",fit_S31,fit_S23,1); 
  sigpdf.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());
  TH2D bgHist("bgHist","2d bkgd hist from data sideband",9,0,3,9,0,3);
  bkgtree->Draw("fit_S23:fit_S31>>bgHist","","goff");
  ((BdkDalitzBase*) sigpdf.getPdf())->weightBins(&bgHist,true,0.0005);
  BdkDalitzHist bkgpdf("bkgpdf", "bkgpdf", BdkDalitzBase::D0, 
		       BdkDalitzBase::PPP0, fit_S31, fit_S23, bgHist);
  RooRealVar fsig("fsig","fsig",0.98);
  RooRealVar fbkg("fbkg","fbkg",1.0-fsig.getVal());
  RooAddPdf mypdf("mypdf","mypdf",RooArgList(*sigpdf.getPdf(),bkgpdf),
		  RooArgList(fsig,fbkg));

  if(doNorm==true) {               
    sigpdf.parameters(); 
    BdkDDalitzAmp::normalizeAll(); 
    sigpdf.parameters().writeToFile("parsdata.txt");  
  }
  else { sigpdf.parameters().readFromFile("parsdata.txt"); }




  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  mypdf.setIntegratorConfig(*cfg);

  RooArgSet allPars = sigpdf.parameters();
  RooRealVar* var1 = allPars.find("sigpdf.pdf.dalitzAmp.Nonres_amp");
  RooRealVar* var2 = allPars.find("sigpdf.pdf.dalitzAmp.Nonres_phase");
  RooRealVar* var3 = allPars.find("sigpdf.pdf.dalitzAmp.Rho+_amp");
  RooRealVar* var4 = allPars.find("sigpdf.pdf.dalitzAmp.Rho+_phase");
  RooRealVar* var5 = allPars.find("sigpdf.pdf.dalitzAmp.Rho-_amp");
  RooRealVar* var6 = allPars.find("sigpdf.pdf.dalitzAmp.Rho-_phase");
  RooRealVar* var7 = allPars.find("sigpdf.pdf.dalitzAmp.Rho0_amp");
  RooRealVar* var8 = allPars.find("sigpdf.pdf.dalitzAmp.Rho0_phase");

  var1->setConstant(false);
  var2->setConstant(false);
  var5->setConstant(false);
  var6->setConstant(false);
  var7->setConstant(false);
  var8->setConstant(false);

  TString name = "PhaseScan";
  if(NR_180 == true) { var2->setVal(257.0); name += "_NR_180"; }
  if(NR_Sign == true) { var2->setVal(-77.0); name += "_NR_Sign"; }
  if(RhoM_180 == true) { var6->setVal(176.0); name += "_Rho-_180"; }
  if(RhoM_Sign == true) { var6->setVal(4.0); name += "_Rho-_Sign"; }
  if(Rho0_180 == true) { var8->setVal(190.0); name += "_Rho0_180"; }
  if(Rho0_Sign == true) { var8->setVal(-10.0); name += "_Rho0_Sign"; }


  cout << "------------------------------------------------" << endl;
  cout << "$$$$$$$$$$$$$$$$$  now doing " << name << endl;
  cout << "------------------------------------------------" << endl;


  RooFitResult* fitresult = mypdf.fitTo(*data,"mr");
  resultFile.cd();
  fitresult->Write(name);

  Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;
  RooRealVar totalm("totalm","totalm",total);
  RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));
  BdkPdfDDalitz sigpdf12("sigpdf12","sigpdf12",fit_S31a,fit_S23,1); 
  TH2D bgHist12("bgHist12","2d bkgd hist from data sideband",6,0,3,6,0,3);
  bkgtree->Draw("fit_S23:fit_S12>>bgHist12","","goff");
  ((BdkDalitzBase*) sigpdf12.getPdf())->weightBins(&bgHist12,true,0.0005);
  BdkDalitzHist bkgpdf12("bkgpdf12", "bkgpdf12", BdkDalitzBase::D0, 
			 BdkDalitzBase::PPP0, fit_S23, fit_S12, bgHist12);
  RooAddPdf mypdf12("mypdf12","mypdf12",RooArgList(*sigpdf12.getPdf(),bkgpdf12),
		    RooArgList(fsig,fbkg));



  //Make the plots
  //******************************************************
  RooPlot* xframe = fit_S23.frame(0,3);
  data->plotOn(xframe,MarkerSize(0.1),DrawOption("z"));
  mypdf.plotOn(xframe,Project(fit_S31));
  xframe->getAttLine()->SetLineWidth(1);
  xframe->getAttLine()->SetLineStyle(1);
  xframe->SetTitleSize(0.03, "Y" );
  xframe->SetTitleOffset(1.8, "Y");

  RooPlot* yframe = fit_S31.frame(0,3);
  data->plotOn(yframe,MarkerSize(0.1),DrawOption("z"));
  mypdf.plotOn(yframe,Project(fit_S23)); 
  yframe->getAttLine()->SetLineWidth(1);
  yframe->getAttLine()->SetLineStyle(1);
  yframe->SetTitleSize(0.03, "Y" );
  yframe->SetTitleOffset(1.8, "Y");

  RooPlot* zframe = fit_S12.frame(0,3);
  data->plotOn(zframe,MarkerSize(0.1),DrawOption("z"));
  mypdf12.plotOn(zframe,Project(fit_S23)); 
  zframe->getAttLine()->SetLineWidth(1);
  zframe->getAttLine()->SetLineStyle(1);
  zframe->SetTitleSize(0.03, "Y" );
  zframe->SetTitleOffset(1.8, "Y");

  TCanvas *canvas = new TCanvas("canvas","allevents",900,400);
  canvas->Divide(3,1);
  canvas->cd(1);xframe->Draw();
  canvas->cd(2);yframe->Draw();
  canvas->cd(3);
  zframe->Draw();
  canvas->SaveAs(TString("./PhaseScan/")+name+TString(".gif"));
  canvas->SaveAs(TString("./PhaseScan/")+name+TString(".eps"));
  delete canvas;

}  //end the macro







