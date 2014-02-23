// does data Dalitz plot fit and also plots 
// Chi distribution of fitted DP - data DP

{
gROOT->Reset();
bool doNorm = true;

//calculate the dalitz boundary
Double_t md0 = 1.8645;
Double_t mpi0=0.1349766;
Double_t mpi = 0.13957;


RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",0,3); //This is pi+ pi0  rho+
RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",0,3); //This is pi- pi0  rho-
RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",0,3); //This is pi+ pi-
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
RooAddPdf mypdf("mypdf","mypdf",RooArgList(*sigpdf.getPdf(),bkgpdf),RooArgList(fsig,fbkg));

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

// var1->setRange("var1Range", 0.7, 1.3);
// var2->setRange("var2Range", 55., 95.); 
// var5->setRange("var5Range", 0.5, 0.8);
//var6->setVal(0.0);
//var6->setRange("var6Range", -10., 10.);
// var7->setRange("var3Range", 0.2, 0.4);
// var8->setRange("var4Range", 0., 20.);



mypdf.fitTo(*data,"m");


 //
  // Plot chi^2
  //

TH2D h("h","m^{2}(#pi^{-}#pi^{0}) : m^{2}(#pi^{+}#pi^{0})",60,0,3,60,0,3);
sigtree->Draw("fit_S23:fit_S31>>h","","goff");
TH2F* hFit = fit_S31.createHistogram("hFit",fit_S23,0,0,0,0);
mypdf->fillHistogram(hFit,RooArgList(fit_S31,fit_S23));
Double_t hintegral = h.Integral();
Double_t hFitintegral = hFit->Integral();
Double_t Normalization = hintegral/hFitintegral;
hFit->Scale(Normalization);

  TH2D *hchi2 = new TH2D("hchi2","",60,0,3,60,0,3);  
  TH1D* pullhist = new TH1D("pullhist","",25,-4,4); 
  int ndof = 60*60;
  ndof--;   // histos are normalized to each other

  Double_t chi2Total = 0;    // total chi2
  for (int x=1; x <= 60; x++) {
    for (int y=1; y <= 60; y++) {
      
      Double_t bin1 =  h.GetBinContent(x,y);
      Double_t bin2 = hFit->GetBinContent(x,y); 

      if (bin1==0 || bin2==0) ndof--;    // no data -> one less dof
      else {
        Double_t sqrtChi2 = (bin1-bin2);
        sqrtChi2 /= sqrt(h->GetBinError(x,y)**2 
			 + hFit->GetBinError(x,y)**2);        
        hchi2->SetBinContent(x,y,sqrtChi2);
	pullhist->Fill(sqrtChi2);
      }
    }
  }


  gStyle->SetOptStat(0);
  TCanvas* aCanvas = new TCanvas("c3", "c3", 900, 600);
  aCanvas->Divide(2,2);
  aCanvas->cd(1);
  h.Draw("colz");
  aCanvas->cd(2);
  hFit->Draw("colz");
  aCanvas->cd(3);
  hchi2->SetMaximum(4);
  hchi2->SetMinimum(-4);
  hchi2->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hchi2->Draw("colz");
  aCanvas->cd(4);  
  pullhist->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  pullhist->Fit("gaus");
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  pullhist->Draw();
  aCanvas->Draw();
  aCanvas->SaveAs("DP_chi_NR249.eps");
  aCanvas->SaveAs("DP_chi_NR249.gif");
//   aCanvas->SaveAs("DP_chi_NR77.eps");
//   aCanvas->SaveAs("DP_chi_NR77.gif");





/*

Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;
RooRealVar totalm("totalm","totalm",total);
RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));
BdkPdfDDalitz sigpdf12("sigpdf12","sigpdf12",fit_S31a,fit_S23,1); 
TH2D bgHist12("bgHist12","2d bkgd hist from data sideband",9,0,3,9,0,3);
bkgtree->Draw("fit_S23:fit_S12>>bgHist12","","goff");
((BdkDalitzBase*) sigpdf12.getPdf())->weightBins(&bgHist12,true,0.0005);
BdkDalitzHist bkgpdf12("bkgpdf12", "bkgpdf12", BdkDalitzBase::D0, 
		     BdkDalitzBase::PPP0, fit_S23, fit_S12, bgHist12);
RooAddPdf mypdf12("mypdf12","mypdf12",RooArgList(*sigpdf12.getPdf(),bkgpdf12),
		  RooArgList(fsig,fbkg));



//Make the plots
// ******************************************************
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
canvas->SaveAs("newdalitz.gif");
canvas->SaveAs("newdalitz.eps");
delete canvas;

*/

}  //end the macro







