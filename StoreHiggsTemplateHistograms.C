

const double MINRange = 140.0;
const double MAXRange = 620.0;
const int BINWIDTH = 15;
const int NBINSFORPDF = (int)((MAXRange-MINRange)/BINWIDTH);
const char* selection = "Mass2j_PFCor>65. && Mass2j_PFCor<95. && abs(JetPFCor_Pt[1]/Mass2j_PFCor-0.5)<0.2 && gdevtt";


void StoreHiggsTemplateHistograms() {
   grabDataSubtractedHistograms();
   getHiggsMasspoint(160);
   getHiggsMasspoint(180);
   getHiggsMasspoint(200);
   getHiggsMasspoint(250);
   getHiggsMasspoint(300);
   getHiggsMasspoint(350);
   getHiggsMasspoint(400);
   getHiggsMasspoint(450);
   getHiggsMasspoint(500);
   getHiggsMasspoint(550);
   getHiggsMasspoint(600);
}




void grabDataSubtractedHistograms() {
  TFile f("Histograms_data_and_template.root", "update");



  TFile* fitFile = new TFile( "mLnuJJ-combined-fit.root", "read");
  TCanvas* fitCan = (TCanvas*) fitFile->Get( "mLnuJJ-combined-fit" );
  RooHist* data = (RooHist*) fitCan->FindObject( "h_data" );
  RooCurve* fit = (RooCurve*) fitCan->FindObject( "h_total" );
  RooCurve* wjj = (RooCurve*) fitCan->FindObject( "h_Wjets" );


  TFile* subtrFile = new TFile( "mLnuJJ-combined-fit-subtracted.root", "read");
  TCanvas* subtrCan = (TCanvas*) subtrFile->Get( "mLnuJJ-combined-fit-subtracted" );
  RooHist* subtrHist = (RooHist*) subtrCan->FindObject( "resid_h_data_h_SM" );


  TFile* subtrWjetsFile = new TFile( "mLnuJJ-combined-fit-subtractedWjets.root", "read");
  TCanvas* subtrWjetsCan = (TCanvas*) subtrWjetsFile->Get( "mLnuJJ-combined-fit-subtractedWjets" );
  RooHist* subtrWjetsHist = (RooHist*) subtrWjetsCan->FindObject( "resid_h_data_h_Background" );
  RooCurve* Diboson = (RooCurve*) subtrWjetsCan->FindObject( "h_diboson" );


 f.cd();
 data->Write("hist_data");
 fit->Write("curve_fit");
 wjj->Write("curve_WJets");
 subtrHist->Write("hist_data_AllSubtracted");
 subtrWjetsHist->Write("hist_data_WJetsSubtracted");
 Diboson->Write("curve_diboson");
 f.Close();
}



void getHiggsMasspoint(int mHiggs) {

   TFile f("Histograms_data_and_template.root", "update");
   char temp[100];
   sprintf(temp, "data/ReducedTree/RD_WmunuJets_CMSSW415-Spring11MC_WWToLNuQQ_M-%d.root", mHiggs);
   TFile* mhfile = new TFile( temp, "read");
   TTree* treeTemp = (TTree*) mhfile->Get("WJet");
   sprintf(temp, "Higgs_Mass_%d", mHiggs);

   TH1* th1H = new TH1D(temp, temp, NBINSFORPDF, MINRange, MAXRange);
   treeTemp->Draw( TString("fit_mlvjj>>")+TString(temp), selection,"goff");
   th1H->Scale( 1. / th1H->Integral() );

   f.cd();
   th1H->Write();
   f.Close();
}
