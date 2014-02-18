#include "/uscms/home/kalanand/Utility/utilities.cc"



void ProdValidation(float min=0.0, float max = 500.0, int bins=50) {


  TString ZmassWindow = "abs(mZee-91.2)<10.0";
  TString etaCut = " && abs(JetRecoEta[5][0])<1.3";
  TString secondJetCut = "&& JetRecoDphi[5][0] > 2.94 && JetRecoPt[5][1]/Z_Pt < 0.1";


  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetNdivisions ( 505 ,"X");

  TFile* file = TFile::Open("SD_Ele15_Zee.root");
  TTree* tree = (TTree*)file->Get("ZJet");
  TFile* file2 = TFile::Open("SD_Ele15_QCD.root");
  TTree* tree2 = (TTree*)file2->Get("ZJet");
  TFile* file3 = TFile::Open("SD_Ele15_Wenu.root");
  TTree* tree3 = (TTree*)file3->Get("ZJet");
  TFile* file4 = TFile::Open("SD_Ele15_TTbar.root");
  TTree* tree4 = (TTree*)file4->Get("ZJet");


  TH1F* ECALIsoZee=AddTwoHistos( *file, "EmECALIsolation", "EPECALIsolation");
  TH1F* HCALIsoZee =AddTwoHistos( *file, "EmHCALIsolation", "EPHCALIsolation");
  TH1F* TrackIsoZee=AddTwoHistos(*file,"EmTrackIsolation","EPTrackIsolation");
  TH1F* EoverPZee = AddTwoHistos(*file, "EmEoverP", "EpEoverP");

  TH1F* ECALIsoQCD=AddTwoHistos( *file2, "EmECALIsolation", "EPECALIsolation");
  TH1F* HCALIsoQCD=AddTwoHistos( *file2, "EmHCALIsolation", "EPHCALIsolation");
  TH1F* TrackIsoQCD=AddTwoHistos(*file2,"EmTrackIsolation","EPTrackIsolation");
  TH1F* EoverPQCD = AddTwoHistos(*file2, "EmEoverP", "EpEoverP");

  TH1F* ECALIsoW=AddTwoHistos( *file3, "EmECALIsolation", "EPECALIsolation");
  TH1F* HCALIsoW=AddTwoHistos( *file3, "EmHCALIsolation", "EPHCALIsolation");
  TH1F* TrackIsoW=AddTwoHistos(*file3,"EmTrackIsolation","EPTrackIsolation");
  TH1F* EoverPW = AddTwoHistos(*file3, "EmEoverP", "EpEoverP");

  TH1F* ECALIsoTT=AddTwoHistos( *file4, "EmECALIsolation", "EPECALIsolation");
  TH1F* HCALIsoTT=AddTwoHistos( *file4, "EmHCALIsolation", "EPHCALIsolation");
  TH1F* TrackIsoTT=AddTwoHistos(*file4,"EmTrackIsolation","EPTrackIsolation");
  TH1F* EoverPTT = AddTwoHistos(*file4, "EmEoverP", "EpEoverP");

  SetHistAttributes( *ECALIsoZee, 1, "Electron ECAL isolation");
  SetHistAttributes( *HCALIsoZee, 1, "Electron HCAL isolation");
  SetHistAttributes( *TrackIsoZee, 1, "Electron track-based isolation");
  SetHistAttributes( *EoverPZee, 1, "Electron E/p");
  SetHistAttributes( *ECALIsoQCD, 2);
  SetHistAttributes( *HCALIsoQCD, 2);
  SetHistAttributes( *TrackIsoQCD, 2);
  SetHistAttributes( *EoverPQCD, 2);
  SetHistAttributes( *ECALIsoW, 4);
  SetHistAttributes( *HCALIsoW, 4);
  SetHistAttributes( *TrackIsoW, 4);
  SetHistAttributes( *EoverPW, 4);
  SetHistAttributes( *ECALIsoTT, 6);
  SetHistAttributes( *HCALIsoTT, 6);
  SetHistAttributes( *TrackIsoTT, 6);
  SetHistAttributes( *EoverPTT, 6);


  TH1F* Zmass_noIso      = (TH1F* ) file->Get("Zmass_noIso");
  TH1F* htemp            = (TH1F* ) file2->Get("Zmass_noIso");
  Zmass_noIso->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_noIso");
  Zmass_noIso->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_noIso");
  Zmass_noIso->Add(htemp);

  TH1F* Zmass_TkIso      = (TH1F* ) file->Get("Zmass_TkIso");
  htemp                  = (TH1F* ) file2->Get("Zmass_TkIso");
  Zmass_TkIso->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_TkIso");
  Zmass_TkIso->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_TkIso");
  Zmass_TkIso->Add(htemp);


  TH1F* Zmass_TkECALIso  = (TH1F* ) file->Get("Zmass_TkECALIso");
  htemp                  = (TH1F* ) file2->Get("Zmass_TkECALIso");
  Zmass_TkECALIso->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_TkECALIso");
  Zmass_TkECALIso->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_TkECALIso");
  Zmass_TkECALIso->Add(htemp);


  TH1F* Zmass_TkECALHCALIso = (TH1F* ) file->Get("Zmass_TkECALHCALIso");
  htemp                  = (TH1F* ) file2->Get("Zmass_TkECALHCALIso");
  Zmass_TkECALHCALIso->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_TkECALHCALIso");
  Zmass_TkECALHCALIso->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_TkECALHCALIso");
  Zmass_TkECALHCALIso->Add(htemp);


  TH1F* Zmass_loose      = (TH1F* ) file->Get("Zmass_loose");
  htemp                  = (TH1F* ) file2->Get("Zmass_loose");
  Zmass_loose->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_loose");
  Zmass_loose->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_loose");
  Zmass_loose->Add(htemp);


  TH1F* Zmass_tight      = (TH1F* ) file->Get("Zmass_tight");    
  htemp                  = (TH1F* ) file2->Get("Zmass_tight");
  Zmass_tight->Add(htemp);
  htemp                  = (TH1F* ) file3->Get("Zmass_tight");
  Zmass_tight->Add(htemp);
  htemp                  = (TH1F* ) file4->Get("Zmass_tight");
  Zmass_tight->Add(htemp);



  SetHistAttributes( *Zmass_noIso, 1, "Reconstructed Z mass");
  SetHistAttributes( *Zmass_TkIso, 2);
  SetHistAttributes( *Zmass_TkECALIso, 4);
  SetHistAttributes( *Zmass_TkECALHCALIso, 3);
  SetHistAttributes( *Zmass_loose, 2);
  SetHistAttributes( *Zmass_tight, 3);
  Zmass_TkIso->SetFillColor(2);
  Zmass_TkECALIso->SetFillColor(4);
  Zmass_TkECALHCALIso->SetFillColor(3);
  Zmass_loose->SetFillColor(2);
  Zmass_tight->SetFillColor(3);


  TH1D* EptZee   = CreateHistogram("EptZee",20, 0.0, 100.0, 1,
				 "Electron p_{T} (GeV/c)");
  TH1D* EptQCD   = CreateHistogram("EptQCD",20, 0.0, 100.0, 2);
  TH1D* EptW     = CreateHistogram("EptW",20, 0.0, 100.0, 4);
  TH1D* EptTT    = CreateHistogram("EptTT",20, 0.0, 100.0, 6);

  TH1D* EetaZee  = CreateHistogram("EetaZee", 30, -3, 3, 1,"Electron #eta");
  TH1D* EetaQCD  = CreateHistogram("EetaQCD", 30, -3, 3, 2);
  TH1D* EetaW    = CreateHistogram("EetaW", 30, -3, 3, 4);
  TH1D* EetaTT   = CreateHistogram("EetaTT", 30, -3, 3, 6);

  TH1D* EphiZee  = CreateHistogram("EphiZee", 30, -3, 3, 1,"Electron #phi");
  TH1D* EphiQCD  = CreateHistogram("EphiQCD", 30, -3, 3, 2);
  TH1D* EphiW    = CreateHistogram("EphiW", 30, -3, 3, 4);
  TH1D* EphiTT   = CreateHistogram("EphiTT", 30, -3, 3, 6);

  TH1D* EdhZee   = CreateHistogram("EdhZee", 100, -0.01, 0.01, 1,
				 "Electron #Delta#eta");
  TH1D* EdhQCD   = CreateHistogram("EdhQCD", 100, -0.01, 0.01, 2);
  TH1D* EdhW     = CreateHistogram("EdhW", 100, -0.01, 0.01, 4);
  TH1D* EdhTT    = CreateHistogram("EdhTT", 100, -0.01, 0.01, 6);

  TH1D* EdfiZee  = CreateHistogram("EdfiZee", 100, -0.03, 0.03, 1, 
				"Electron #Delta#phi");
  TH1D* EdfiQCD  = CreateHistogram("EdfiQCD", 100,-0.03, 0.03, 2);
  TH1D* EdfiW    = CreateHistogram("EdfiW", 100,-0.03, 0.03, 4);
  TH1D* EdfiTT   = CreateHistogram("EdfiTT", 100,-0.03, 0.03, 6);

  TH1D* EshhZee  = CreateHistogram("EshhZee", 350, 0, 0.035, 1, 
				"Electron #sigma_{i#eta i#eta}");
  TH1D* EshhQCD  = CreateHistogram("EshhQCD", 350, 0, 0.035, 2);
  TH1D* EshhW  = CreateHistogram("EshhW", 350, 0, 0.035, 4);
  TH1D* EshhTT  = CreateHistogram("EshhTT", 350, 0, 0.035, 6);


  tree->Draw("eMinusEt>>EptZee", ZmassWindow, "goff");
  tree->Draw("ePlusEt>>+EptZee", ZmassWindow, "goff");
  tree2->Draw("eMinusEt>>EptQCD", ZmassWindow, "goff");
  tree2->Draw("ePlusEt>>+EptQCD", ZmassWindow, "goff");
  tree3->Draw("eMinusEt>>EptW", ZmassWindow, "goff");
  tree3->Draw("ePlusEt>>+EptW", ZmassWindow, "goff");
  tree4->Draw("eMinusEt>>EptTT", ZmassWindow, "goff");
  tree4->Draw("ePlusEt>>+EptTT", ZmassWindow, "goff");


  tree->Draw("eMinusEta>>EetaZee", ZmassWindow, "goff");
  tree->Draw("ePlusEta>>+EetaZee", ZmassWindow, "goff");
  tree2->Draw("eMinusEta>>EetaQCD", ZmassWindow, "goff");
  tree2->Draw("ePlusEta>>+EetaQCD", ZmassWindow, "goff");
  tree3->Draw("eMinusEta>>EetaW", ZmassWindow, "goff");
  tree3->Draw("ePlusEta>>+EetaW", ZmassWindow, "goff");
  tree4->Draw("eMinusEta>>EetaTT", ZmassWindow, "goff");
  tree4->Draw("ePlusEta>>+EetaTT", ZmassWindow, "goff");


  tree->Draw("eMinusPhi>>EphiZee", ZmassWindow, "goff");
  tree->Draw("ePlusPhi>>+EphiZee", ZmassWindow, "goff");
  tree2->Draw("eMinusPhi>>EphiQCD", ZmassWindow, "goff");
  tree2->Draw("ePlusPhi>>+EphiQCD", ZmassWindow, "goff");
  tree3->Draw("eMinusPhi>>EphiW", ZmassWindow, "goff");
  tree3->Draw("ePlusPhi>>+EphiW", ZmassWindow, "goff");
  tree4->Draw("eMinusPhi>>EphiTT", ZmassWindow, "goff");
  tree4->Draw("ePlusPhi>>+EphiTT", ZmassWindow, "goff");


  tree->Draw("eMinus_DeltaEtaIn>>EdhZee", ZmassWindow, "goff");
  tree->Draw("ePlus_DeltaEtaIn>>+EdhZee", ZmassWindow, "goff");
  tree2->Draw("eMinus_DeltaEtaIn>>EdhQCD", ZmassWindow, "goff");
  tree2->Draw("ePlus_DeltaEtaIn>>+EdhQCD", ZmassWindow, "goff");
  tree3->Draw("eMinus_DeltaEtaIn>>EdhW", ZmassWindow, "goff");
  tree3->Draw("ePlus_DeltaEtaIn>>+EdhW", ZmassWindow, "goff");
  tree4->Draw("eMinus_DeltaEtaIn>>EdhTT", ZmassWindow, "goff");
  tree4->Draw("ePlus_DeltaEtaIn>>+EdhTT", ZmassWindow, "goff");


  tree->Draw("eMinus_DeltaPhiIn>>EdfiZee", ZmassWindow, "goff");
  tree->Draw("ePlus_DeltaPhiIn>>+EdfiZee", ZmassWindow, "goff");
  tree2->Draw("eMinus_DeltaPhiIn>>EdfiQCD", ZmassWindow, "goff");
  tree2->Draw("ePlus_DeltaPhiIn>>+EdfiQCD", ZmassWindow, "goff");
  tree3->Draw("eMinus_DeltaPhiIn>>EdfiW", ZmassWindow, "goff");
  tree3->Draw("ePlus_DeltaPhiIn>>+EdfiW", ZmassWindow, "goff");
  tree4->Draw("eMinus_DeltaPhiIn>>EdfiTT", ZmassWindow, "goff");
  tree4->Draw("ePlus_DeltaPhiIn>>+EdfiTT", ZmassWindow, "goff");


  tree->Draw("eMinus_SigmaEtaEta>>EshhZee", ZmassWindow, "goff");
  tree->Draw("ePlus_SigmaEtaEta>>+EshhZee", ZmassWindow, "goff");
  tree2->Draw("eMinus_SigmaEtaEta>>EshhQCD", ZmassWindow, "goff");
  tree2->Draw("ePlus_SigmaEtaEta>>+EshhQCD", ZmassWindow, "goff");
  tree3->Draw("eMinus_SigmaEtaEta>>EshhW", ZmassWindow, "goff");
  tree3->Draw("ePlus_SigmaEtaEta>>+EshhW", ZmassWindow, "goff");
  tree4->Draw("eMinus_SigmaEtaEta>>EshhTT", ZmassWindow, "goff");
  tree4->Draw("ePlus_SigmaEtaEta>>+EshhTT", ZmassWindow, "goff");



  TH1D* responseZee = CreateHistogram("responseZee",40,0.0,2.0,1,
		      "p_{T}^{jet} / p_{T}^{Z}");
  TH1D* responseQCD = CreateHistogram("responseQCD", 40, 0.0, 2.0, 2);
  TH1D* responseW = CreateHistogram("responseW", 40, 0.0, 2.0, 4);
  TH1D* responseTT = CreateHistogram("responseTT", 40, 0.0, 2.0, 6);



  TH1D* JetptZee = CreateHistogram("JetptZee", bins, min, max, 1, 
				   "jet p_{T} (GeV/c)");
  TH1D* JetptQCD = CreateHistogram("JetptQCD", bins, min, max, 2);
  TH1D* JetptW = CreateHistogram("JetptW", bins, min, max, 4);
  TH1D* JetptTT = CreateHistogram("JetptTT", bins, min, max, 6);


  TH1D* ZptZee = CreateHistogram("ZptZee", bins, min, max,1,"Z p_{T} (GeV/c)");
  TH1D* ZptQCD = CreateHistogram("ZptQCD",bins, min, max, 2);
  TH1D* ZptW = CreateHistogram("ZptW",bins, min, max, 4);
  TH1D* ZptTT = CreateHistogram("ZptTT",bins, min, max, 6);


  TH1D* dPhiZee = CreateHistogram("dPhiZee", 35, 0, 3.5, 1, 
				  "#Delta#phi (Z, jet)");
  TH1D* dPhiQCD = CreateHistogram("dPhiQCD", 35, 0, 3.5, 2);
  TH1D* dPhiW = CreateHistogram("dPhiW", 35, 0, 3.5, 4);
  TH1D* dPhiTT = CreateHistogram("dPhiTT", 35, 0, 3.5, 6);


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  tree->Draw("JetRecoDphi[5][0]>>dPhiZee", ZmassWindow + etaCut, "goff");
  tree2->Draw("JetRecoDphi[5][0]>>dPhiQCD", ZmassWindow + etaCut, "goff");
  tree3->Draw("JetRecoDphi[5][0]>>dPhiW", ZmassWindow + etaCut, "goff");
  tree4->Draw("JetRecoDphi[5][0]>>dPhiTT", ZmassWindow + etaCut, "goff");


  tree->Draw("JetRecoPt[5][0]>>JetptZee", ZmassWindow + etaCut + secondJetCut, "goff");
  tree2->Draw("JetRecoPt[5][0]>>JetptQCD", ZmassWindow + etaCut + secondJetCut, "goff");
  tree3->Draw("JetRecoPt[5][0]>>JetptW", ZmassWindow + etaCut + secondJetCut, "goff");
  tree4->Draw("JetRecoPt[5][0]>>JetptTT", ZmassWindow + etaCut + secondJetCut, "goff");


  tree->Draw("Z_Pt>>ZptZee", ZmassWindow + etaCut + secondJetCut, "goff");
  tree2->Draw("Z_Pt>>ZptQCD", ZmassWindow + etaCut + secondJetCut, "goff");
  tree3->Draw("Z_Pt>>ZptW", ZmassWindow + etaCut + secondJetCut, "goff");
  tree4->Draw("Z_Pt>>ZptTT", ZmassWindow + etaCut + secondJetCut, "goff");


  tree->Draw("JetRecoPt[5][0]/Z_Pt>>responseZee", ZmassWindow + etaCut + secondJetCut, "goff");
  tree2->Draw("JetRecoPt[5][0]/Z_Pt>>responseQCD", ZmassWindow + etaCut + secondJetCut, "goff");
  tree3->Draw("JetRecoPt[5][0]/Z_Pt>>responseW", ZmassWindow + etaCut + secondJetCut, "goff");
  tree4->Draw("JetRecoPt[5][0]/Z_Pt>>responseTT", ZmassWindow + etaCut + secondJetCut, "goff");



  // Now make all the plots

  ECALIsoZee->GetXaxis()->SetRangeUser(0, 10);
  PlotFourOnCanvas("tempPlots/ecaliso", *ECALIsoZee,"Z#rightarrow ee", 
		   *ECALIsoQCD,"QCD Bkg.", *ECALIsoW,"W#rightarrow e#nu",
		   *ECALIsoTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

  HCALIsoZee->GetXaxis()->SetRangeUser(0, 10);
  TCanvas* c = 
    PlotFourOnCanvas("tempPlots/hcaliso", *HCALIsoZee,"Z#rightarrow ee", 
		     *HCALIsoQCD,"QCD Bkg.", *HCALIsoW,"W#rightarrow e#nu",
		     *HCALIsoTT,"ttbar", 0.55, 0.7, 0.9, 0.9);
  c->SetLogy(1);

  TrackIsoZee->GetXaxis()->SetRangeUser(0.5, 10);
  PlotFourOnCanvas("tempPlots/trackiso", *TrackIsoZee,"Z#rightarrow ee", 
		   *TrackIsoQCD,"QCD Bkg.", *TrackIsoW,"W#rightarrow e#nu",
		   *TrackIsoTT,"ttbar", 0.55, 0.7, 0.9, 0.9);


  TCanvas* can4 = new TCanvas("can4","",500,500);
  gStyle->SetOptStat(0);
  Zmass_noIso->Draw("e");
  Zmass_TkIso->Draw("HIST same");
  Zmass_TkECALIso->Draw("HIST same");
  Zmass_TkECALHCALIso->Draw("HIST same");
  leg_hist = new TLegend(0.6,0.6,0.85,0.85);
  leg_hist->SetHeader("Isolation requirements");
  leg_hist->AddEntry( Zmass_noIso,"No isolation","le");
  leg_hist->AddEntry( Zmass_TkIso,"Track iso","F");
  leg_hist->AddEntry( Zmass_TkECALIso,"Track+ECAL iso","F");
  leg_hist->AddEntry( Zmass_TkECALHCALIso,"Track+ECAL+HCAL iso","F");
  leg_hist->SetMargin(0.12);
  leg_hist->SetFillColor(0);
  leg_hist->Draw();
  can4->SaveAs("tempPlots/Zmass_Isolation.eps");
  can4->SaveAs("tempPlots/Zmass_Isolation.gif");
  can4->SaveAs("tempPlots/Zmass_Isolation.root");


  TCanvas* can5 = new TCanvas("can5","",500,500);
  gStyle->SetOptStat(0);
  Zmass_noIso->Draw("e");
  Zmass_loose->Draw("HIST same");
  Zmass_tight->Draw("HIST same");
  leg_hist = new TLegend(0.6,0.6,0.85,0.85);
  leg_hist->AddEntry( Zmass_noIso,"No Id","le");
  leg_hist->AddEntry( Zmass_loose,"Loose Id","F");
  leg_hist->AddEntry( Zmass_tight,"Tight Id","F");
  leg_hist->SetMargin(0.12);
  leg_hist->SetFillColor(0);
  leg_hist->Draw();
  can5->SaveAs("tempPlots/Zmass_Id.eps");
  can5->SaveAs("tempPlots/Zmass_Id.gif");
  can5->SaveAs("tempPlots/Zmass_Id.root");


  EoverPZee->GetXaxis()->SetRangeUser(0.5, 2);
  PlotFourOnCanvas("tempPlots/EoverP", *EoverPZee,"Z#rightarrow ee", 
		   *EoverPQCD,"QCD Bkg.", *EoverPW,"W#rightarrow e#nu",
		   *EoverPTT,"ttbar", 0.55, 0.7, 0.9, 0.9);
  PlotFourOnCanvas("tempPlots/deta", *EdhZee,"Z#rightarrow ee", 
		   *EdhQCD,"QCD Bkg.", *EdhW,"W#rightarrow e#nu",
		   *EdhTT,"ttbar", 0.6, 0.7, 0.95, 0.9);

  PlotFourOnCanvas("tempPlots/dphi", *EdfiZee,"Z#rightarrow ee", 
		   *EdfiQCD,"QCD Bkg.", *EdfiW,"W#rightarrow e#nu",
		   *EdfiTT,"ttbar", 0.65, 0.7, 0.95, 0.9);

  PlotFourOnCanvas("tempPlots/shh", *EshhZee,"Z#rightarrow ee", 
		   *EshhQCD,"QCD Bkg.", *EshhW,"W#rightarrow e#nu",
		   *EshhTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

 
  PlotFourOnCanvas("tempPlots/electronpt", *EptZee,"Z#rightarrow ee", 
		   *EptQCD,"QCD Bkg.", *EptW,"W#rightarrow e#nu",
		   *EptTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

  PlotFourOnCanvas("tempPlots/electroneta", *EetaZee,"Z#rightarrow ee", 
		   *EetaQCD,"QCD Bkg.", *EetaW,"W#rightarrow e#nu",
		   *EetaTT,"ttbar", 0.65, 0.75, 0.95, 0.95);


  EphiZee->GetYaxis()->SetMoreLogLabels();
  EphiZee->GetYaxis()->SetNoExponent();
  TCanvas* c2 = 
    PlotFourOnCanvas("tempPlots/electronphi", *EphiZee,"Z#rightarrow ee", 
		     *EphiQCD,"QCD Bkg.", *EphiW,"W#rightarrow e#nu",
		     *EphiTT,"ttbar", 0.55, 0.7, 0.9, 0.9);
  c2->SetLogy(1);
  c2->SetLogy(1);


  PlotFourOnCanvas("tempPlots/jetpt", *JetptZee,"Z#rightarrow ee", 
		   *JetptQCD,"QCD Bkg.", *JetptW,"W#rightarrow e#nu",
		   *JetptTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

  PlotFourOnCanvas("tempPlots/Zpt", *ZptZee,"Z#rightarrow ee", 
		   *ZptQCD,"QCD Bkg.", *ZptW,"W#rightarrow e#nu",
		   *ZptTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

  PlotFourOnCanvas("tempPlots/response", *responseZee,"Z#rightarrow ee", 
		   *responseQCD,"QCD Bkg.", *responseW,"W#rightarrow e#nu",
		   *responseTT,"ttbar", 0.55, 0.7, 0.9, 0.9);

  PlotFourOnCanvas("tempPlots/dphiZjet", *dPhiZee,"Z#rightarrow ee", 
		   *dPhiQCD,"QCD Bkg.", *dPhiW,"W#rightarrow e#nu",
		   *dPhiTT,"ttbar", 0.35, 0.7, 0.7, 0.9);
}






void SetHistAttributes(TH1& h, int color=1, 
		       TString xtit="", TString ytit="") {
  h.SetLineColor(color);
  h.SetMarkerColor(color);
  h.SetLineWidth(2);
  SetTitles( h, xtit, ytit);
} 





TH1F* AddTwoHistos(TFile& file, char* hist1, char* hist2) {
  TH1F* h1  = (TH1F* ) file.Get(hist1);
  TH1F* h2  = (TH1F* ) file.Get(hist2);
  h1->Add(h2);
  return h1;
}
