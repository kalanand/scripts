{
gROOT->Reset();
gStyle->Reset();
gROOT->SetStyle("BABAR");
int numBINS = 30;
bool doNorm = true;

RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0}) [GeV^{2}/c^{4}]",0,3);
RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0}) [GeV^{2}/c^{4}]",0,3);
RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+}) [GeV^{2}/c^{4}]",0,3); 
fit_S23.setBins(numBINS);
fit_S31.setBins(numBINS);
fit_S12.setBins(numBINS);

BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
dalitzCfg->getParameters(RooArgSet())->readFromFile("../BToDKTo3piK/Dstar/dalitzCfg.par");
dalitzCfg->getParameters(RooArgSet())->Print("v");


// define signal pdf
BdkPdfDDalitz dstar("dstar","dstar",fit_S23,fit_S31,BdkDalitzBase::D0);
RooArgSet allPars = dstar.parameters();
dstar.parameters().readFromFile("PPP_NominalFit.par");
if(doNorm==true) {               
  dstar.parameters(); 
  BdkDDalitzAmp::normalizeAll(); 
}


RooArgSet fractions = (RooArgSet&) dstar.BreitWignerNormalizationCoefficients();
fractions.Print("v");
}  







