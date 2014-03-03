// generates toy MC DP, fits it and plots chi distribution

{
  gROOT->Reset();
  bool doNorm = true;
  bool theOther=true;

  Double_t  md0 = 1.8645;
  Double_t  mpi0 = 0.1349766;
  Double_t  mpi = 0.13957;
  Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;

  RooRealVar totalm("totalm","totalm",total);
  RooRealVar Dmass("Dmass","Dmass",1.72,1.99);
  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",0,3); 
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",0,3);
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",0,3); 
//   fit_S23.setBins(24);
//   fit_S31.setBins(24);
//   fit_S12.setBins(24);
  RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));


  //read in the background dataset
  TFile hello("hhPi0.root","READ");
  TTree* tree = (TTree*)hello.Get("pipipi0_ccbar");
  RooArgList myList;
  myList.add(RooArgSet(Dmass,fit_S23,fit_S31,fit_S12));
  gROOT->cd();
  RooDataSet fulldata("fulldata","fulldata",tree, myList);
  TTree* tree1 = tree->CopyTree("Flag==1");
  RooDataSet othdata("othdata","",tree1, myList);
  RooDataSet* d2;
  if(theOther==true) d2 = (RooDataSet*) othdata.reduce("Dmass>1.93"); 
  else d2 = (RooDataSet*) fulldata.reduce("Dmass>1.93"); 
  RooDataHist bkg("bkg","bkg",RooArgSet(fit_S23,fit_S31),*d2);
  RooHistPdf bkgpdf("bkgpdf","bkgpdf",RooArgSet(fit_S23,fit_S31),bkg);


  // define signal pdf
  BdkPdfDDalitz sigpdf("sigpdf","sigpdf",fit_S23,fit_S31,1); 
  RooRealVar fsig("fsig","fsig",0.976);
  RooRealVar fbkg("fbkg","fbkg",1.0-fsig.getVal());
  RooAddPdf mypdf("mypdf","mypdf",RooArgList(*sigpdf.getPdf(),bkgpdf),RooArgList(fsig,fbkg));
  
  if(doNorm==true) {               
    sigpdf.parameters(); 
    BdkDDalitzAmp::normalizeAll(); 
    sigpdf.parameters().writeToFile("pars.txt");  
  }
  else { sigpdf.parameters().readFromFile("pars.txt"); }



                       
  //  signal + background dataset
  RooDataSet *data = sigpdf.generate(20000) ;
  TTree* bgTree = tree->CopyTree("abs(Dmass-1.861)<0.015 && !(Flag==1||Flag==2||Flag==11)","",492,0);
  RooDataSet d3("d3","background in the signal region", bgTree, RooArgSet(fit_S23,fit_S31));
  data->append(d3);

  RooFormulaVar fit_S12a("fit_S12","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S31));
  data->addColumn(fit_S12a);

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
  var3->setVal(1.0);
  var3->setConstant(true);
  var4->setVal(0.0);
  var4->setConstant(true);
  var1->setRange("var1Range", 0., 2.0);
  var2->setRange("var2Range", 60., 100.); 
  var5->setRange("var5Range", 0.0, 1.0);
  var6->setRange("var6Range", -20., 20.);
  var7->setRange("var3Range", 0., 1.0);
  var8->setRange("var4Range", -20., 20.);


                
  mypdf.fitTo(*data,"m");



  BdkPdfDDalitz sigpdf12("sigpdf12","sigpdf12",fit_S23,fit_S31a,1); 
  RooDataHist bkg12("bkg12","bkg12",RooArgSet(fit_S23,fit_S12),*d2);
  RooHistPdf bkgpdf12("bkgpdf12","bkgpdf12",RooArgSet(fit_S23,fit_S12),bkg12);
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

  TCanvas *c1 = new TCanvas("c1","allevents",1000,380);
  c1->Divide(3,1);
  c1->cd(1);xframe->Draw();
  c1->cd(2);yframe->Draw();
  c1->cd(3);
  zframe->Draw();

// if(theOther==true) {
//    c1->SaveAs("DPSigBack_2.gif");
//    c1->SaveAs("DPSigBack_2.eps");
// }
//   else {
//    c1->SaveAs("DPSigBack.gif");
//    c1->SaveAs("DPSigBack.eps");
//  }
//   delete c1;

}  //end the macro

