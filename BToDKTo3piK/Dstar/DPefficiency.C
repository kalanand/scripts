// Binned fit for efficiency on Dalitz plot

DPefficiency(){
  bool doSmooth = false;
TFile inputFile("hhPi0Dalitz.root","READ");
TTree* recoTree = (TTree* )inputFile.Get("pipipi0_signal"));
TFile genFile("pipipi0_truth.root","READ");
TTree* genTree = (TTree* )genFile.Get("ntp1");
TH2D genHist("genHist","",60,0,3,60,0,3);
TH2D effHist("effHist","",60,0,3,60,0,3);
genHist.Sumw2();
effHist.Sumw2();
TCut genCut = TCut("!(Flag==2) && truD0_Pcm>2.77");
TCut recoCut = TCut("!(Flag==2) && abs(Dmass-1.861)<0.015 && abs(hh_SignedDistance)<2.0");
genTree->Draw("tru_S23:tru_S31>>genHist", genCut, "goff");
recoTree->Draw("fit_S23:fit_S31>>effHist", recoCut, "goff");
if(doSmooth==true) Smooth( effHist, genHist);
effHist.Divide(&genHist); 
effHist.Scale(100);


BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
const char* cfgFile = "../BToDKTo3piK/params/dalitzCfg.par";
cout << "Reading Daliz configuration from "<<cfgFile<<endl;
dalitzCfg->getParameters(RooArgSet())->readFromFile(cfgFile);
dalitzCfg->getParameters(RooArgSet())->Print("v");


RooRealVar S31("S31","",0,3);
RooRealVar S23("S23","",0,3);
BdkPdf2DpolyDalitz eff("eff","",S31,S23);
Bdk2DpolyDalitz* mypdf = eff.getPdf();
RooArgSet allPars = eff.parameters();
RooRealVar* var1 = allPars.find("eff.c0");
var1->setVal(1.0);
var1->setConstant(true);


RooDataHist H("H", "H", RooArgList(S31,S23), &effHist); 
RooChi2Var chi2("chi2","chi2",*mypdf,H);
RooMinuit minuit(chi2);
minuit.fit("m");


Double_t Norm = mypdf->getNorm(RooArgSet(S31,S23)); 
TH2F* fittedHist = S31.createHistogram("fittedHist",S23,0,0,0,0);
mypdf->fillHistogram(fittedHist,RooArgList(S31,S23));

Double_t Normalization = 0.0;
Double_t effintegral = effHist.Integral();
Double_t fittedintegral = fittedHist->Integral();
int bins = fittedHist->GetNbinsX()*fittedHist->GetNbinsY();
if(fittedintegral != 0.0) Normalization = effintegral*bins/(3600*fittedintegral);
fittedHist->Scale(Normalization);
if(Norm != 0.0)  Normalization = Normalization/Norm;
cout<< "Normalization = " << Normalization << endl;

TCanvas aCanvas("aCanvas","",900,400);
gStyle->SetOptStat(0);
aCanvas.Divide(2,1);
aCanvas.cd(1);
effHist.Draw("colz");
aCanvas.cd(2);
fittedHist->Draw("colz");
aCanvas.SaveAs("DPeff.gif");
}


void Smooth(TH2& recoHist, TH2& genHist) {

  for(int XBin=0; XBin<60; XBin++) {
    for(int YBin=0; YBin<60; YBin++) {
      Int_t count = 0;
      Double_t reco =0.0, erreco = recoHist.GetBinError( XBin, YBin);
      Double_t gen =0.0, ergen = genHist.GetBinError( XBin, YBin);
      for(int i=-1; i<=1; i++) {
	for(int j=-1; j<=1; j++) {
	  Double_t ent1 = recoHist.GetBinContent( XBin+i, YBin+j);
	  Double_t ent2 = genHist.GetBinContent( XBin+i, YBin+j);
	  if( !(ent1<1 || ent2<1)){ reco += ent1; gen += ent2; count++;}
	}
      } 
      if(!(count==0 || genHist.GetBinContent( XBin, YBin)<1.0) ) {
	recoHist.SetBinContent( XBin, YBin, reco/count );
	recoHist.SetBinError( XBin, YBin, erreco/3.0);
	genHist.SetBinContent( XBin, YBin, gen/count );
	genHist.SetBinError( XBin, YBin, ergen/3.0);
      }  
    }
  }
}



