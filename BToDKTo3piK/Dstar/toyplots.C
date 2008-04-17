// generate toy MC for (1.) signal only,  (2.) signal + bkg, and
// (3.) signal + bkg + K_short
// then fit them using the same pdf


void toyplots() {
  fitToyMc();
  //  fitToyMcWithBackground();
  // fitToyMcWithBackgroundAndKs();
}






void fitToyMcWithBackgroundAndKs() {
gROOT->Reset();
bool doNorm = true;

Double_t  md0 = 1.8645;
Double_t  mpi0 = 0.1349766;
Double_t  mpi = 0.13957;
Double_t upper=pow((md0-mpi),2);
Double_t lower=pow((mpi0+mpi),2);
Double_t s12lower=pow((mpi+mpi),2);
Double_t s12upper=pow((md0-mpi0),2);
Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;

RooRealVar totalm("totalm","totalm",total);
RooRealVar Dmass("Dmass","Dmass",1.72,1.99);
RooRealVar hh_mass("hh_mass","hh_mass",0.0,3.0);
RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",lower,upper); //This is pi+ pi0  rho+
RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",lower,upper); //This is pi- pi0  rho-
RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",s12lower,s12upper); //This is pi+ pi-

RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));


//read in the background dataset
TFile hello("hhPi0.root","READ");
TTree* tree = (TTree*)hello.Get("pipipi0_ccbar");
RooArgList myList;
myList.add(RooArgSet(Dmass,hh_mass,fit_S23,fit_S31,fit_S12));
gROOT->cd();
RooDataSet fulldata("fulldata","fulldata",tree, myList);
RooDataSet *d2 = fulldata.reduce("Dmass>1.93"); 
RooDataHist bkg("bkg","bkg",RooArgSet(fit_S23,fit_S31),*d2);
RooHistPdf bkgpdf("bkgpdf","bkgpdf",RooArgSet(fit_S23,fit_S31),bkg);
RooDataSet *d3 = fulldata.reduce("hh_mass>0.489 && hh_mass<0.508"); 
RooDataHist ks("ks","ks",RooArgSet(fit_S23,fit_S31),*d3);
RooHistPdf kspdf("kspdf","kspdf",RooArgSet(fit_S23,fit_S31),ks);


// define signal pdf
BdkPdfDDalitz sigpdf("sigpdf","sigpdf",fit_S23,fit_S31,1); 
RooRealVar fsig("fsig","fsig",0.956, 0.94, 0.97);
RooRealVar fks("fks","fks",0.019, 0.014, 0.024);
RooRealVar fbkg("fbkg","fbkg", 0.025, 0.02,0.035);
RooAddPdf mypdf("mypdf","mypdf",RooArgList(*sigpdf.getPdf(),bkgpdf,kspdf),RooArgList(fsig,fbkg,fks));
  
if(doNorm==true) {               
  sigpdf.parameters(); 
  BdkDDalitzAmp::normalizeAll(); 
  sigpdf.parameters().writeToFile("pars.txt");  
}
else { sigpdf.parameters().readFromFile("pars.txt"); }



                       
// Generate dataset from pdf
RooDataSet *data = mypdf.generate(RooArgSet(fit_S23,fit_S31),20000) ;
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

RooDataHist ks12("ks12","ks12",RooArgSet(fit_S23,fit_S12),*d3);
RooHistPdf kspdf12("kspdf12","kspdf12",RooArgSet(fit_S23,fit_S12),ks12);

RooAddPdf mypdf12("mypdf12","mypdf12",RooArgList(*sigpdf12.getPdf(),bkgpdf12,kspdf12),
		  RooArgList(fsig,fbkg,fks));
 

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

TCanvas *c1 = new TCanvas("c1","allevents",1200,400);
c1->Divide(3,1);
c1->cd(1);xframe->Draw();
c1->cd(2);yframe->Draw();
c1->cd(3);
zframe->Draw();
c1->SaveAs("toydalitz_withKs.gif");
c1->SaveAs("toydalitz_withKs.eps");
delete c1;
}



void fitToyMcWithBackground() {

  gROOT->Reset();
  bool doNorm = true;

  Double_t  md0 = 1.8645;
  Double_t  mpi0 = 0.1349766;
  Double_t  mpi = 0.13957;
  Double_t upper=pow((md0-mpi),2);
  Double_t lower=pow((mpi0+mpi),2);
  Double_t s12lower=pow((mpi+mpi),2);
  Double_t s12upper=pow((md0-mpi0),2);
  Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;

  RooRealVar totalm("totalm","totalm",total);
  RooRealVar Dmass("Dmass","Dmass",1.72,1.99);
  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",lower,upper); //This is pi+ pi0  rho+
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",lower,upper); //This is pi- pi0  rho-
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",s12lower,s12upper); //This is pi+ pi-
  RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));


  //read in the background dataset
  TFile hello("hhPi0.root","READ");
  TTree* tree = (TTree*)hello.Get("pipipi0_ccbar");
  RooArgList myList;
  myList.add(RooArgSet(Dmass,fit_S23,fit_S31,fit_S12));
  gROOT->cd();
  RooDataSet fulldata("fulldata","fulldata",tree, myList);
  RooDataSet *d2 = fulldata.reduce("Dmass>1.93"); 
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



                       
  // Generate dataset from pdf
  RooDataSet *data = mypdf.generate(RooArgSet(fit_S23,fit_S31),20000) ;
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

  TCanvas *c1 = new TCanvas("c1","allevents",1200,400);
  c1->Divide(3,1);
  c1->cd(1);xframe->Draw();
  c1->cd(2);yframe->Draw();
  c1->cd(3);
  zframe->Draw();
  c1->SaveAs("toydalitz_withBkg.gif");
  c1->SaveAs("toydalitz_withBkg.eps");
  delete c1;

}  //end the macro











void fitToyMc() {

  gROOT->Reset();
  bool doNorm = true;

  Double_t  md0 = 1.8645;
  Double_t  mpi0 = 0.1349766;
  Double_t  mpi = 0.13957;
  Double_t upper=pow((md0-mpi),2);
  Double_t lower=pow((mpi0+mpi),2);
  Double_t s12lower=pow((mpi+mpi),2);
  Double_t s12upper=pow((md0-mpi0),2);
  Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;

  RooRealVar totalm("totalm","totalm",total);
  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",lower,upper); //This is pi+ pi0  rho+
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",lower,upper); //This is pi- pi0  rho-
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",s12lower,s12upper); //This is pi+ pi-
  RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));



  // define signal pdf
  BdkPdfDDalitz sigpdf("sigpdf","sigpdf",fit_S23,fit_S31,1);   
  if(doNorm==true) {               
    sigpdf.parameters(); 
    BdkDDalitzAmp::normalizeAll(); 
    sigpdf.parameters().writeToFile("pars.txt");  
  }
  else { sigpdf.parameters().readFromFile("pars.txt"); }



                       
  // Generate dataset from pdf
  RooDataSet *data = sigpdf.generate(20000) ;
  RooFormulaVar fit_S12a("fit_S12","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S31));
  data->addColumn(fit_S12a);

  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  sigpdf.getPdf()->setIntegratorConfig(*cfg);

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


                
  sigpdf.getPdf()->fitTo(*data,"m");
  BdkPdfDDalitz sigpdf12("sigpdf12","sigpdf12",fit_S23,fit_S31a,1); 


  //Make the plots
  //******************************************************
  RooPlot* xframe = fit_S23.frame(0,3);
  data->plotOn(xframe,MarkerSize(0.1),DrawOption("z"));
  sigpdf.getPdf()->plotOn(xframe,Project(fit_S31));
  xframe->getAttLine()->SetLineWidth(1);
  xframe->getAttLine()->SetLineStyle(1);
  xframe->SetTitleSize(0.03, "Y" );
  xframe->SetTitleOffset(1.8, "Y");

  RooPlot* yframe = fit_S31.frame(0,3);
  data->plotOn(yframe,MarkerSize(0.1),DrawOption("z"));
  sigpdf.getPdf()->plotOn(yframe,Project(fit_S23)); 
  yframe->getAttLine()->SetLineWidth(1);
  yframe->getAttLine()->SetLineStyle(1);
  yframe->SetTitleSize(0.03, "Y" );
  yframe->SetTitleOffset(1.8, "Y");

  RooPlot* zframe = fit_S12.frame(0,3);
  data->plotOn(zframe,MarkerSize(0.1),DrawOption("z"));
  sigpdf12.getPdf()->plotOn(zframe,Project(fit_S23)); 
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
   c1->SaveAs("toydalitz.gif");
  c1->SaveAs("toydalitz.eps");
  delete c1;

}  //end the macro











