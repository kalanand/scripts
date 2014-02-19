
void makeplots ()
{

   //const double intLUMI = 135.;
   const double intLUMI = 250.;
  const double WJets_scale = 31500.* intLUMI/15161497;
  // const double WW_scale = (27.79 * intLUMI/2050240) / 0.1075 / 0.676;
  const double WW_scale = (43.* intLUMI/2050240) / 0.1075 / 0.676;
  // correcting for BR
  // const double WZ_scale = (10.4 * intLUMI/2194752) / 0.1075 / 0.6991;
  const double WZ_scale = (18. * intLUMI/2194752) / 0.1075 / 0.6991;
  // correcting for BR
  const double TT_scale = 65.83 * intLUMI/4816224;


  const double normFactor = 1.6;

  TFile *file0 = TFile::Open("Data_WenuJets_2011_3MayDCS.root");
  TTree* tree = (TTree*) file0->Get("WJet");

  TFile *file1 = TFile::Open("MC_WJets_enujj_Spring11.root");
  TTree* treeMC_WJets = (TTree*) file1->Get("WJet");

  TFile *file2 = TFile::Open("MC_WW_enujj_Spring11.root");
  TTree* treeMC_WW_enujj = (TTree*) file2->Get("WJet");

  TFile *file3 = TFile::Open("MC_WZ_enujj_Spring11.root");
  TTree* treeMC_WZ_enujj = (TTree*) file3->Get("WJet");

  TFile *file4 = TFile::Open("MC_TT_enujj_Spring2011.root");
  TTree* treeMC_TT = (TTree*) file4->Get("WJet");


//////////// binning
  const int jjmass_NBINS = 25;
  const double jjmass_MIN = 30.;
  const double jjmass_MAX = 230.;


//   const int jjmass_NBINS = 15;
//   const double jjmass_MIN = 70.;
//   const double jjmass_MAX = 100.;


  const int jjdeta_NBINS = 40;
  const double jjdeta_MIN = -5.;
  const double jjdeta_MAX = 5.;

  const int jjdphi_NBINS = 64;
  const double jjdphi_MIN = 0;
  const double jjdphi_MAX = 3.14159265358979312;


///////////// cuts 

  //TCut goodMET("W_electron_et>25. && event_met_pfmetsignificance>3.");
  TCut goodMET("W_electron_et>30. && (event_met_pfmet-0.5*event_nPV)>25. && W_electron_isWP80==1");
  TCut goodW(goodMET && "W_mt>50.");
  TCut twojets("numPFCorJets==2");
  TCut costheta("cosJacksonAngleV2j_PFCor>0.9");
  TCut cosjj("cosJacksonAngle2j_PFCor>0.6");

  TCut twoBJets("numPFCorJetBTags==2");
  TCut noBJets("numPFCorJetBTags==0");

  TCut dphicutWjj("abs(abs(abs(W_phi-atan2((JetPFCor_Py[0]+JetPFCor_Py[1]),(JetPFCor_Px[0]+JetPFCor_Px[1])))-TMath::Pi())-TMath::Pi())>2.94");

  TCut allCuts( goodW && twojets && costheta && cosjj && noBJets && dphicutWjj );

// define n-1 cuts
  TCut nmoc_met("W_electron_et>30. && W_electron_isWP80==1 && W_mt>50."  && twojets && cosjj && noBJets && dphicutWjj );
  TCut nmoc_mt( goodMET && twojets && cosjj && noBJets && dphicutWjj );
  TCut nmoc_cos( goodW && twojets && cosjj && noBJets && dphicutWjj );
  TCut nmoc_cosjj( goodW && twojets && costheta && noBJets && dphicutWjj );
  TCut nmoc_dphicutWjj( goodW && twojets && costheta && cosjj && noBJets);
  TCut nmoc_noBJets( goodW && twojets && costheta && cosjj && dphicutWjj );



///////////// variables
  char* mjj = "Mass2j_PFCor";
  char* jjdeta = "JetPFCor_Eta[0]-JetPFCor_Eta[1]";
  char* jjcostheta = "tanh(0.5*(JetPFCor_Y[0]-JetPFCor_Y[1]))";
  char* jjdphi = "abs(abs(abs(JetPFCor_Phi[0]-JetPFCor_Phi[1])-TMath::Pi())-TMath::Pi())";
  char* eMETdphi = "abs(abs(abs(W_electron_phi-event_met_pfmetPhi)-TMath::Pi())-TMath::Pi())";
  char* ej1dphi = "abs(abs(abs(W_electron_phi-JetPFCor_Phi[0])-TMath::Pi())-TMath::Pi())";
  char* ej2dphi = "abs(abs(abs(W_electron_phi-JetPFCor_Phi[1])-TMath::Pi())-TMath::Pi())";

  char* dphiWjj = "abs(abs(abs(W_phi-atan2((JetPFCor_Py[0]+JetPFCor_Py[1]),(JetPFCor_Px[0]+JetPFCor_Px[1])))-TMath::Pi())-TMath::Pi())";
// ################# How many events we started with ?
  int nEvents = tree->Draw( "W_mt", goodW, "goff");;
  cout << "# events in the W+jj skim = " << nEvents << endl;


// ################# How many good W+jj are there ?
  nEvents = tree->Draw( "W_mt", goodW, "goff");
  cout << "Number of W+jj events = " << nEvents << endl;



  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadRightMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  gStyle->SetOptFit(1111);
  // tdrStyle->SetOptStat("mer"); 
  char temp[1000];



  TCanvas* c_wmt = new TCanvas("c_wmt","W mT", 500, 500);
  TH1D* h_wmt = new TH1D("h_wmt","", 40, 40, 120);
  tree->Draw("W_mt>>h_wmt", nmoc_mt, "goff");
  h_wmt->GetXaxis()->SetTitle("W m_{T} [GeV]");
  h_wmt->GetYaxis()->SetTitle("Events / 2 GeV");
  TH1D* h_wmt_WJets = new TH1D("h_wmt_WJets","", 40, 40, 120);
  h_wmt_WJets->SetLineColor(4);
  h_wmt_WJets->SetLineWidth(2);
  treeMC_WJets->Draw( "W_mt>>h_wmt_WJets", nmoc_mt, "goff");
  h_wmt_WJets->Scale( WJets_scale );
  h_wmt->SetMaximum( 1.3 * h_wmt_WJets->GetMaximum() );
  h_wmt->Draw("e");
  h_wmt_WJets->Draw("same");
  c_wmt->SaveAs("W-mt.gif");





  TCanvas* c_wmet = new TCanvas("c_wmet","W MET", 500, 500);
  TH1D* h_wmet = new TH1D("h_wmet","", 60, 0, 120);
  tree->Draw( "(event_met_pfmet-0.5*event_nPV)>>h_wmet", nmoc_met, "e");
  h_wmet->GetXaxis()->SetTitle("MET [GeV]");
  h_wmet->GetYaxis()->SetTitle("Events / 2 GeV");
  TH1D* h_wmet_WJets = new TH1D("h_wmet_WJets","", 60, 0, 120);
  h_wmet_WJets->SetLineColor(4);
  h_wmet_WJets->SetLineWidth(2);
  treeMC_WJets->Draw( "(event_met_pfmet-0.5*event_nPV)>>h_wmet_WJets", nmoc_met, "goff");
  h_wmet_WJets->Scale( WJets_scale );
  h_wmet_WJets->Draw("same");
  TArrow* ar2 = new TArrow(25,30,25,200,0.03,"<|");
  ar2->SetLineColor(2);
  ar2->SetFillColor(2);
  ar2->Draw();
  c_wmet->SaveAs("W-met.gif");



//   TCanvas* c_wy = new TCanvas("c_wy","y of neutrino", 500, 500);
//   TH1D* h_wy = new TH1D("h_wy","", 40, -4, 4);
//   tree->Draw( "W_y>>h_wy", goodWn, "e");
//   h_wy->GetXaxis()->SetTitle("y of neutrino");
//   h_wy->GetYaxis()->SetTitle("Events / 0.2");
//   TArrow* ar3 = new TArrow(-2.6,80,-2.6,140,0.03,"<|");
//   ar3->SetLineColor(2);
//   ar3->SetFillColor(2);
//   ar3->Draw();
//   TArrow* ar32 = new TArrow(2.6,80,2.6,140,0.03,"<|");
//   ar32->SetLineColor(2);
//   ar32->SetFillColor(2);
//   ar32->Draw();
//   c_wy->SaveAs("W-rapidity.gif");




  TCanvas* c_jjmass = new TCanvas("c_jjmass","jj mass", 500, 500);
  TH1D* h_jjmass = new TH1D("h_jjmass","", jjmass_NBINS, jjmass_MIN, jjmass_MAX);
  TH1D* h_jjmass_WJets = new TH1D("h_jjmass_WJets","", jjmass_NBINS, jjmass_MIN, jjmass_MAX);
//  h_jjmass_WJets->SetLineColor(4);
  h_jjmass_WJets->SetFillColor(9);

  h_jjmass_WJets->SetLineWidth(2);
  TH1D* h_jjmass_WW_enujj = new TH1D("h_jjmass_WW_enujj","", jjmass_NBINS, jjmass_MIN, jjmass_MAX);
  h_jjmass_WW_enujj->SetFillColor(kOrange);
  TH1D* h_jjmass_WZ_enujj = h_jjmass_WW_enujj->Clone("h_jjmass_WZ_enujj");

  tree->Draw("1.1*Mass2j_PFCor>>h_jjmass", allCuts, "goff");
  treeMC_WJets->Draw( "Mass2j_PFCor>>h_jjmass_WJets", allCuts, "goff");
  h_jjmass_WJets->Scale( WJets_scale );

  treeMC_WW_enujj->Draw("Mass2j_PFCor>>h_jjmass_WW_enujj", allCuts, "goff");
  h_jjmass_WW_enujj->Scale( WW_scale );
  treeMC_WZ_enujj->Draw("Mass2j_PFCor>>h_jjmass_WZ_enujj", allCuts, "goff");
  h_jjmass_WZ_enujj->Scale( WZ_scale );
  h_jjmass_WW_enujj->Add(h_jjmass_WZ_enujj);
  TH1D* h_jjmass_TT = new TH1D("h_jjmass_TT","", jjmass_NBINS, jjmass_MIN, jjmass_MAX);
  h_jjmass_TT->SetFillColor(6);
  treeMC_TT->Draw( "Mass2j_PFCor>>h_jjmass_TT", allCuts, "goff");
  h_jjmass_TT->Scale( TT_scale );

// #################
  h_jjmass_WJets->Add(h_jjmass_TT);
  TH1D* h_jjmass_WW_clone = h_jjmass_WW_enujj->Clone("h_jjmass_WW_clone");
  h_jjmass_WW_clone->SetLineColor(kOrange);
  h_jjmass_WW_clone->SetLineWidth(3);
  h_jjmass_WW_clone->SetFillColor(0);
  h_jjmass_WW_enujj->Add(h_jjmass_WJets);
// #################

  h_jjmass->SetMaximum( 1.6 * h_jjmass->GetMaximum() );
  h_jjmass->Draw("e");
  h_jjmass_WW_enujj->Draw("same");
  h_jjmass_WJets->Draw("same");
  h_jjmass_TT->Draw("same");
  h_jjmass_WW_clone->Draw("hist same");
  h_jjmass->Draw("esame");
  c_jjmass->RedrawAxis();

  //h_jjmass_WJets->Draw("sames");
  //h_jjmass->SetMaximum( 1.5 * h_jjmass->GetMaximum() );
  h_jjmass->GetXaxis()->SetTitle("m_{jj} [GeV]");
  sprintf(temp, "Events / %d GeV", (jjmass_MAX-jjmass_MIN)/jjmass_NBINS);
  h_jjmass->GetYaxis()->SetTitle(temp);
  h_jjmass->GetXaxis()->SetNdivisions(505);
  c_jjmass->RedrawAxis();
  c_jjmass->SaveAs("jj-mass.gif");


// ################# How many good W+jj are there ?
  cout << "Number of observed events in data = " <<  h_jjmass->Integral() << endl;
  cout << "Number of W+jj events predicted from MC = " <<  h_jjmass_WJets->Integral() << endl;
  cout << "Number of WW events predicted from MC = " <<  h_jjmass_WW_enujj->Integral()/10. << endl;





// // ############## Print out how many events survive W mass cut:70--90
//   h_jjmass->GetXaxis()->SetRangeUser(60., 100.);
//   cout << "Data events in 60--100 GeV mass window = " << h_jjmass->Integral() << endl;
//   h_jjmass_WJets->GetXaxis()->SetRangeUser(60., 100.);
//   cout << "W+jets events in 60--100 GeV mass window = " << h_jjmass_WJets->Integral() << endl;
//   h_jjmass_WW_enujj->GetXaxis()->SetRangeUser(60., 100.);
//   cout << "Signal events in 60--100 GeV mass window = " << h_jjmass_WW_enujj->Integral() << endl;  





  TCanvas* c_jjcostheta = new TCanvas("c_jjcostheta","jj costhetastar", 500, 500);
  TH1D* h_jjcostheta = new TH1D("h_jjcostheta","", 25, -1, 1);
  h_jjcostheta->SetMinimum(0);
  h_jjcostheta->GetXaxis()->SetTitle("cos#theta* (jj, j)");
  h_jjcostheta->GetYaxis()->SetTitle("Events / 0.08");
  TH1D* h_jjcostheta_WJets = new TH1D("h_jjcostheta_WJets","", 25, -1, 1);
  h_jjcostheta_WJets->SetFillColor(9);
  h_jjcostheta_WJets->SetLineWidth(2);

//////////
  TH1D* h_jjcostheta_WW_enujj = new TH1D("h_jjcostheta_WW_enujj","", 25, -1, 1);
  h_jjcostheta_WW_enujj->SetFillColor(kOrange);
  TH1D* h_jjcostheta_WZ_enujj = h_jjcostheta_WW_enujj->Clone("h_jjcostheta_WZ_enujj");

   tree->Draw("cosJacksonAngle2j_PFCor>>h_jjcostheta", nmoc_cosjj, "goff");
  treeMC_WJets->Draw( "cosJacksonAngle2j_PFCor>>h_jjcostheta_WJets", nmoc_cosjj, "goff");
  h_jjcostheta_WJets->Scale( WJets_scale );

  treeMC_WW_enujj->Draw("cosJacksonAngle2j_PFCor>>h_jjcostheta_WW_enujj", nmoc_cosjj, "goff");
  h_jjcostheta_WW_enujj->Scale( WW_scale );
  treeMC_WZ_enujj->Draw("cosJacksonAngle2j_PFCor>>h_jjcostheta_WZ_enujj", nmoc_cosjj, "goff");
  h_jjcostheta_WZ_enujj->Scale( WZ_scale );
  h_jjcostheta_WW_enujj->Add(h_jjcostheta_WZ_enujj);
  TH1D* h_jjcostheta_TT = new TH1D("h_jjcostheta_TT","", 25, -1, 1);
  h_jjcostheta_TT->SetFillColor(6);
  treeMC_TT->Draw( "cosJacksonAngle2j_PFCor>>h_jjcostheta_TT", nmoc_cosjj, "goff");
  h_jjcostheta_TT->Scale( TT_scale );

// #################
  h_jjcostheta_WJets->Add(h_jjcostheta_TT);
  h_jjcostheta_WW_enujj->Add(h_jjcostheta_WJets);
// #################

  h_jjcostheta->SetMaximum( 1.6 * h_jjcostheta->GetMaximum() );
  h_jjcostheta->Draw("e");
  h_jjcostheta_WW_enujj->Draw("same");
  //h_jjcostheta_WJets->Draw("same");
  h_jjcostheta_TT->Draw("same");
  h_jjcostheta->Draw("esame");
  c_jjcostheta->RedrawAxis();
//////////



  h_jjcostheta->Draw("e");
  h_jjcostheta_WJets->Draw("same");
  h_jjcostheta->Draw("esame");
  TArrow* ar25 = new TArrow(0.6, 5, 0.6, 140,0.03,"<|");
  ar25->SetLineColor(2);
  ar25->SetFillColor(2);
  ar25->Draw();
  c_jjcostheta->SaveAs("jj-costhetastar.gif");






  TCanvas* c_costheta = new TCanvas("c_costheta","costhetastar", 500, 500);
  TH1D* h_costheta = new TH1D("h_costheta","", 50, 0.5, 1);
  h_costheta->SetMinimum(0);
  tree->Draw("cosJacksonAngleV2j_PFCor>>h_costheta", nmoc_cos, "goff");
  h_costheta->GetXaxis()->SetTitle("cos#theta* (W+jj, W)");
  h_costheta->GetYaxis()->SetTitle("Events / 0.01");
  TH1D* h_costheta_WJets = new TH1D("h_costheta_WJets","", 50, 0.5, 1);
  h_costheta_WJets->SetFillColor(9);
  h_costheta_WJets->SetLineWidth(2);
  treeMC_WJets->Draw("cosJacksonAngleV2j_PFCor>>h_costheta_WJets", nmoc_cos, "goff");
  h_costheta_WJets->Scale( WJets_scale );
  h_costheta->Draw("e");
  h_costheta_WJets->Draw("same");
  h_costheta->Draw("esame");
  TArrow* ar26 = new TArrow(0.9, 5, 0.9, 140,0.03,"<|");
  ar26->SetLineColor(2);
  ar26->SetFillColor(2);
  ar26->Draw();
  c_costheta->SaveAs("Wjj-costhetastar.gif");








  TCanvas* c_jjdphi = new TCanvas("c_jjdphi","jj #Delta#phi", 500, 500);
  TH1D* h_jjdphi = new TH1D("h_jjdphi","", jjdphi_NBINS/4, jjdphi_MIN, jjdphi_MAX);
  tree->Draw(TString(jjdphi)+TString(">>h_jjdphi"), allCuts, "goff");
  h_jjdphi->GetXaxis()->SetTitle("#Delta#phi of leading two jets");
  sprintf(temp, "Events / %.1f", 4*(jjdphi_MAX-jjdphi_MIN)/jjdphi_NBINS);
  h_jjdphi->GetYaxis()->SetTitle(temp);
  h_jjdphi->SetMinimum(3);
  TH1D* h_jjdphi_WJets = new TH1D("h_jjdphi_WJets","", jjdphi_NBINS/4, jjdphi_MIN, jjdphi_MAX);
  h_jjdphi_WJets->SetFillColor(9);
  h_jjdphi_WJets->SetLineWidth(2);
  treeMC_WJets->Draw( TString(jjdphi)+TString(">>h_jjdphi_WJets"), allCuts, "goff");
  h_jjdphi_WJets->Scale( WJets_scale );
  h_jjdphi->SetMinimum(0);
  h_jjdphi->Draw("e");
  h_jjdphi_WJets->Draw("same");
  h_jjdphi->Draw("esame");
  c_jjdphi->RedrawAxis();
  c_jjdphi->SaveAs("jj-dphi.gif");







  TCanvas* c_j1dphiMET = new TCanvas("c_j1dphiMET","dPhi between leading jet and MET", 500, 500);
  TH1D* h_j1dphiMET = new TH1D("h_j1dphiMET","", jjdphi_NBINS/4, jjdphi_MIN, jjdphi_MAX);
  tree->Draw("abs(JetPFCor_dphiMET[0])>>h_j1dphiMET",allCuts, "goff");
  h_j1dphiMET->GetXaxis()->SetTitle("#Delta#phi between leading jet and MET");
  sprintf(temp, "Events / %.1f", 4*(jjdphi_MAX-jjdphi_MIN)/jjdphi_NBINS);
  h_j1dphiMET->GetYaxis()->SetTitle(temp);
  h_j1dphiMET->SetMinimum(3);
  TH1D* h_j1dphiMET_WJets = new TH1D("h_j1dphiMET_WJets","", jjdphi_NBINS/4, jjdphi_MIN, jjdphi_MAX);
  h_j1dphiMET_WJets->SetFillColor(4);
  h_j1dphiMET_WJets->SetLineWidth(2);
  treeMC_WJets->Draw( "abs(JetPFCor_dphiMET[0])>>h_j1dphiMET_WJets",allCuts, "goff");
  h_j1dphiMET_WJets->Scale( WJets_scale );
  h_j1dphiMET->SetMinimum(0);
  h_j1dphiMET->Draw("e");
  h_j1dphiMET_WJets->Draw("same");
  c_j1dphiMET->RedrawAxis();
  c_j1dphiMET->SaveAs("j1dphiMET.gif");






  TCanvas* c_dphiWjj = new TCanvas("c_dphiWjj","dPhi between two W's", 500, 500);
  TH1D* h_dphiWjj = new TH1D("h_dphiWjj","", jjdphi_NBINS, jjdphi_MIN, jjdphi_MAX);
  tree->Draw( TString(dphiWjj)+TString(">>h_dphiWjj"), nmoc_dphicutWjj, "goff");
  h_dphiWjj->GetXaxis()->SetTitle("#Delta#phi between two W's");
  sprintf(temp, "Events / %.2f", (jjdphi_MAX-jjdphi_MIN)/jjdphi_NBINS);
  h_dphiWjj->GetYaxis()->SetTitle(temp);
  TH1D* h_dphiWjj_WJets = new TH1D("h_dphiWjj_WJets","", jjdphi_NBINS, jjdphi_MIN, jjdphi_MAX);
  h_dphiWjj_WJets->SetFillColor(9);
  h_dphiWjj_WJets->SetLineWidth(2);
  treeMC_WJets->Draw( TString(dphiWjj)+TString(">>h_dphiWjj_WJets"), nmoc_dphicutWjj, "goff");
  h_dphiWjj_WJets->Scale( WJets_scale );
  h_dphiWjj->GetXaxis()->SetRangeUser(1.5,3.2);
  h_dphiWjj->SetMinimum(0);
  h_dphiWjj->Draw("e");
  h_dphiWjj_WJets->Draw("same");
  c_dphiWjj->RedrawAxis();
  TArrow* ar5 = new TArrow(2.94, 5, 2.94, 160,0.03,"<|");
  ar5->SetLineColor(2);
  ar5->SetFillColor(2);
  ar5->Draw();
  c_dphiWjj->SaveAs("dphiWjj.gif");












  TCanvas* c_numTags = new TCanvas("c_numTags","number of b-tags", 500, 500);
  TH1D* h_numTags = new TH1D("h_numTags","", 3, -0.5, 2.5);
  tree->Draw("numPFCorJetBTags>>h_numTags", nmoc_noBJets, "goff");
  h_numTags->GetXaxis()->SetTitle("Number of B-tags");
  h_numTags->SetMinimum(0.6);
  TH1D* h_numTags_WJets = new TH1D("h_numTags_WJets","", 3, -0.5, 2.5);
  h_numTags_WJets->SetLineColor(4);
  h_numTags_WJets->SetLineWidth(2);
  treeMC_WJets->Draw("numPFCorJetBTags>>h_numTags_WJets", nmoc_noBJets, "goff");
  h_numTags_WJets->Scale( WJets_scale );
  h_numTags->SetMaximum( 10 * h_numTags_WJets->GetMaximum() );
  h_numTags->Draw("e");
  h_numTags_WJets->Draw("same");
  c_numTags->SetLogy();
  c_numTags->SaveAs("nBTags.gif");

// // ################# Print the three bin contents
//   h_numTags->Print("all");


//   tdrStyle->SetPadRightMargin(0.13);
//   TCanvas* c_mTVsMjj = new TCanvas("c_mTVsMjj","W mT vs mjj", 500, 500);
//   TH2D* h_mTVsMjj = new TH2D("h_mTVsMjj","", 30, 0., 300.,  20, 40., 200.);
//   h_mTVsMjj->GetXaxis()->SetTitle("m_{jj}    ");
//   h_mTVsMjj->GetXaxis()->SetNdivisions(505);
//   h_mTVsMjj->GetYaxis()->SetTitle("    m_{T}");
//   h_mTVsMjj->GetYaxis()->SetNdivisions(505);
//   sprintf(temp, "W_mt:%s>>h_mTVsMjj", mjj);
//   gStyle->SetPalette(1);
//   tree->Draw(temp, goodW && twojets && noBJets && dphicut, "lego2z");
//   // c_mTVsMjj->RedrawAxis();
//   c_mTVsMjj->SaveAs("mTVsMjj-lego.gif");


}


