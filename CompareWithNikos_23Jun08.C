 void CompareWithNikos_23Jun08(){


   /*
     {
     gROOT->ProcessLine(".L ~/tdrstyle.C");
     setTDRStyle();
     TFile f2("histos.root");
     TH1F *eff = (TH1F*)f2.Get("gsfelectron_gsfelectronIso_Eta_reReco");
     //
     TFile f("./Nikos-reReco/sc_eta_iso.root");

     // common for all sequences ...........................
     TH1F *base = (TH1F*)f.Get("base");
     TH1F *destination = (TH1F*)f.Get("destination");

     TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
     g1->BayesDivide(destination, base, "");
     g1->GetYaxis()->SetRangeUser(0.5, 1.05);
     eff->SetLineColor(2);

     eff->SetMarkerStyle(22);
     eff->SetMarkerSize(1.4);
     eff->SetMarkerColor(2);


     //
     g1->GetYaxis()->SetTitle("Efficiency");
     g1->GetXaxis()->SetTitle("#eta");
     g1->Draw("APE");
     eff->Draw("esame");
     gPad->Update();
   */



  TString stringname1 = "supercluster_gsfelectron";
  TString stringname2 = "gsfelectron_gsfelectronIso";
  TString stringname3 = "gsfelectronIso_gsfelectronIsoLooseID";
  TString stringname4 = "gsfelectronIsoLooseID_HLT";


  // TFile myhistos("histos_sbs.root"); 
  // TFile myhistos("histos_truth.root"); 
  TFile myhistos("histos.root"); 
  TString name11 = stringname1 + TString("_Eta_reReco");
  TString name12 = stringname1 + TString("_Eta_Reco");
  TString name13 = stringname1 + TString("_Phi_reReco");
  TString name14 = stringname1 + TString("_Phi_Reco");
  TString name15 = stringname1 + TString("_Pt_reReco");
  TString name16 = stringname1 + TString("_Pt_Reco");
  TString name17 = stringname1 + TString("_Pt_Eta_reReco");
  TString name18 = stringname1 + TString("_Pt_Eta_Reco");
  TString name19 = stringname1 + TString("_Phi_Eta_reReco");
  TString name10 = stringname1 + TString("_Phi_Eta_Reco");

  TH1F* h11 = (TH1F*) myhistos.Get(name11);
  TH1F* h12 = (TH1F*) myhistos.Get(name12);
  TH1F* h13 = (TH1F*) myhistos.Get(name13);
  TH1F* h14 = (TH1F*) myhistos.Get(name14);
  TH1F* h15 = (TH1F*) myhistos.Get(name15);
  TH1F* h16 = (TH1F*) myhistos.Get(name16);
  TH2F* h17 = (TH2F*) myhistos.Get(name17);
  TH2F* h18 = (TH2F*) myhistos.Get(name18);
  TH2F* h19 = (TH2F*) myhistos.Get(name19);
  TH2F* h10 = (TH2F*) myhistos.Get(name10);

  TFile f11("./Nikos-reReco/sc_eta_recoEle.root");
  TFile f12("./Nikos-Reco/sc_eta_recoEle.root");
  TFile f13("./Nikos-reReco/sc_phi_recoEle.root");
  TFile f14("./Nikos-Reco/sc_phi_recoEle.root");
  TFile f15("./Nikos-reReco/sc_et_recoEle.root");
  TFile f16("./Nikos-Reco/sc_et_recoEle.root");
  TH1F* base11        = (TH1F*)f11.Get("base");
  TH1F* destination11 = (TH1F*)f11.Get("destination");
  TH1F* base12        = (TH1F*)f12.Get("base");
  TH1F* destination12 = (TH1F*)f12.Get("destination");
  TH1F* base13        = (TH1F*)f13.Get("base");
  TH1F* destination13 = (TH1F*)f13.Get("destination");
  TH1F* base14        = (TH1F*)f14.Get("base");
  TH1F* destination14 = (TH1F*)f14.Get("destination");
  TH1F* base15        = (TH1F*)f15.Get("base");
  TH1F* destination15 = (TH1F*)f15.Get("destination");
  TH1F* base16        = (TH1F*)f16.Get("base");
  TH1F* destination16 = (TH1F*)f16.Get("destination");

//   Process(*h11, *base11, *destination11);
//   Process(*h12, *base12, *destination12);
//   Process(*h13, *base13, *destination13);
//   Process(*h14, *base14, *destination14);
//   Process(*h15, *base15, *destination15);
//   Process(*h16, *base16, *destination16);

  TString name21 = stringname2 + TString("_Eta_reReco");
  TString name22 = stringname2 + TString("_Eta_Reco");
  TString name23 = stringname2 + TString("_Phi_reReco");
  TString name24 = stringname2 + TString("_Phi_Reco");
  TString name25 = stringname2 + TString("_Pt_reReco");
  TString name26 = stringname2 + TString("_Pt_Reco");
  TString name27 = stringname2 + TString("_Pt_Eta_reReco");
  TString name28 = stringname2 + TString("_Pt_Eta_Reco");
  TString name29 = stringname2 + TString("_Phi_Eta_reReco");
  TString name20 = stringname2 + TString("_Phi_Eta_Reco");

  TH1F* h21 = (TH1F*) myhistos.Get(name21);
  TH1F* h22 = (TH1F*) myhistos.Get(name22);
  TH1F* h23 = (TH1F*) myhistos.Get(name23);
  TH1F* h24 = (TH1F*) myhistos.Get(name24);
  TH1F* h25 = (TH1F*) myhistos.Get(name25);
  TH1F* h26 = (TH1F*) myhistos.Get(name26);
  TH2F* h27 = (TH2F*) myhistos.Get(name27);
  TH2F* h28 = (TH2F*) myhistos.Get(name28);
  TH2F* h29 = (TH2F*) myhistos.Get(name29);
  TH2F* h20 = (TH2F*) myhistos.Get(name20);

  TFile f21("./Nikos-reReco/sc_eta_iso.root");
  TFile f22("./Nikos-Reco/sc_eta_iso.root");
  TFile f23("./Nikos-reReco/sc_phi_iso.root");
  TFile f24("./Nikos-Reco/sc_phi_iso.root");
  TFile f25("./Nikos-reReco/sc_et_iso.root");
  TFile f26("./Nikos-Reco/sc_et_iso.root");
  TH1F* base21        = (TH1F*)f21.Get("base");
  TH1F* destination21 = (TH1F*)f21.Get("destination");
  TH1F* base22        = (TH1F*)f22.Get("base");
  TH1F* destination22 = (TH1F*)f22.Get("destination");
  TH1F* base23        = (TH1F*)f23.Get("base");
  TH1F* destination23 = (TH1F*)f23.Get("destination");
  TH1F* base24        = (TH1F*)f24.Get("base");
  TH1F* destination24 = (TH1F*)f24.Get("destination");
  TH1F* base25        = (TH1F*)f25.Get("base");
  TH1F* destination25 = (TH1F*)f25.Get("destination");
  TH1F* base26        = (TH1F*)f26.Get("base");
  TH1F* destination26 = (TH1F*)f26.Get("destination");

//   Process(*h21, *base21, *destination21);
//   Process(*h22, *base22, *destination22);
//   Process(*h23, *base23, *destination23);
//   Process(*h24, *base24, *destination24);
//   Process(*h25, *base25, *destination25);
//   Process(*h26, *base26, *destination26);


  TString name31 = stringname3 + TString("_Eta_reReco");
  TString name32 = stringname3 + TString("_Eta_Reco");
  TString name33 = stringname3 + TString("_Phi_reReco");
  TString name34 = stringname3 + TString("_Phi_Reco");
  TString name35 = stringname3 + TString("_Pt_reReco");
  TString name36 = stringname3 + TString("_Pt_Reco");
  TString name37 = stringname3 + TString("_Pt_Eta_reReco");
  TString name38 = stringname3 + TString("_Pt_Eta_Reco");
  TString name39 = stringname3 + TString("_Phi_Eta_reReco");
  TString name30 = stringname3 + TString("_Phi_Eta_Reco");

  TH1F* h31 = (TH1F*) myhistos.Get(name31);
  TH1F* h32 = (TH1F*) myhistos.Get(name32);
  TH1F* h33 = (TH1F*) myhistos.Get(name33);
  TH1F* h34 = (TH1F*) myhistos.Get(name34);
  TH1F* h35 = (TH1F*) myhistos.Get(name35);
  TH1F* h36 = (TH1F*) myhistos.Get(name36);
  TH2F* h37 = (TH2F*) myhistos.Get(name37);
  TH2F* h38 = (TH2F*) myhistos.Get(name38);
  TH2F* h39 = (TH2F*) myhistos.Get(name39);
  TH2F* h30 = (TH2F*) myhistos.Get(name30);

  TFile f31("./Nikos-reReco/sc_eta_id_robust.root");
  TFile f32("./Nikos-Reco/sc_eta_id_robust.root");
  TFile f33("./Nikos-reReco/sc_phi_id_robust.root");
  TFile f34("./Nikos-Reco/sc_phi_id_robust.root");
  TFile f35("./Nikos-reReco/sc_et_id_robust.root");
  TFile f36("./Nikos-Reco/sc_et_id_robust.root");
  TH1F* base31        = (TH1F*)f31.Get("base");
  TH1F* destination31 = (TH1F*)f31.Get("destination");
  TH1F* base32        = (TH1F*)f32.Get("base");
  TH1F* destination32 = (TH1F*)f32.Get("destination");
  TH1F* base33        = (TH1F*)f33.Get("base");
  TH1F* destination33 = (TH1F*)f33.Get("destination");
  TH1F* base34        = (TH1F*)f34.Get("base");
  TH1F* destination34 = (TH1F*)f34.Get("destination");
  TH1F* base35        = (TH1F*)f35.Get("base");
  TH1F* destination35 = (TH1F*)f35.Get("destination");
  TH1F* base36        = (TH1F*)f36.Get("base");
  TH1F* destination36 = (TH1F*)f36.Get("destination");




  TString name41 = stringname4 + TString("_Eta_reReco");
  TString name42 = stringname4 + TString("_Eta_Reco");
  TString name43 = stringname4 + TString("_Phi_reReco");
  TString name44 = stringname4 + TString("_Phi_Reco");
  TString name45 = stringname4 + TString("_Pt_reReco");
  TString name46 = stringname4 + TString("_Pt_Reco");
  TString name47 = stringname4 + TString("_Pt_Eta_reReco");
  TString name48 = stringname4 + TString("_Pt_Eta_Reco");
  TString name49 = stringname4 + TString("_Phi_Eta_reReco");
  TString name40 = stringname4 + TString("_Phi_Eta_Reco");

  TH1F* h41 = (TH1F*) myhistos.Get(name41);
  TH1F* h42 = (TH1F*) myhistos.Get(name42);
  TH1F* h43 = (TH1F*) myhistos.Get(name43);
  TH1F* h44 = (TH1F*) myhistos.Get(name44);
  TH1F* h45 = (TH1F*) myhistos.Get(name45);
  TH1F* h46 = (TH1F*) myhistos.Get(name46);
  TH2F* h47 = (TH2F*) myhistos.Get(name47);
  TH2F* h48 = (TH2F*) myhistos.Get(name48);
  TH2F* h49 = (TH2F*) myhistos.Get(name49);
  TH2F* h40 = (TH2F*) myhistos.Get(name40);


  TFile f41("./Nikos-reReco/sc_eta_trigger_idRobust.root");
  TFile f42("./Nikos-Reco/sc_eta_trigger_idRobust.root");
  TFile f43("./Nikos-reReco/sc_phi_trigger_idRobust.root");
  TFile f44("./Nikos-Reco/sc_phi_trigger_idRobust.root");
  TFile f45("./Nikos-reReco/sc_et_trigger_idRobust.root");
  TFile f46("./Nikos-Reco/sc_et_trigger_idRobust.root");
  TH1F* base41        = (TH1F*)f41.Get("base");
  TH1F* destination41 = (TH1F*)f41.Get("destination");
  TH1F* base42        = (TH1F*)f42.Get("base");
  TH1F* destination42 = (TH1F*)f42.Get("destination");
  TH1F* base43        = (TH1F*)f43.Get("base");
  TH1F* destination43 = (TH1F*)f43.Get("destination");
  TH1F* base44        = (TH1F*)f44.Get("base");
  TH1F* destination44 = (TH1F*)f44.Get("destination");
  TH1F* base45        = (TH1F*)f45.Get("base");
  TH1F* destination45 = (TH1F*)f45.Get("destination");
  TH1F* base46        = (TH1F*)f46.Get("base");
  TH1F* destination46 = (TH1F*)f46.Get("destination");


  // Make 1D comparison plots

  makeplots1D( *h11, *base11, *destination11, name11);
  makeplots1D( *h12, *base12, *destination12, name12);
  makeplots1D( *h13, *base13, *destination13, name13);
  makeplots1D( *h14, *base14, *destination14, name14);
  makeplots1D( *h15, *base15, *destination15, name15);
  makeplots1D( *h16, *base16, *destination16, name16);

  makeplots1D( *h21, *base21, *destination21, name21);
  makeplots1D( *h22, *base22, *destination22, name22);
  makeplots1D( *h23, *base23, *destination23, name23);
  makeplots1D( *h24, *base24, *destination24, name24);
  makeplots1D( *h25, *base25, *destination25, name25);
  makeplots1D( *h26, *base26, *destination26, name26);

  makeplots1D( *h31, *base31, *destination31, name31);
  makeplots1D( *h32, *base32, *destination32, name32);
  makeplots1D( *h33, *base33, *destination33, name33);
  makeplots1D( *h34, *base34, *destination34, name34);
  makeplots1D( *h35, *base35, *destination35, name35);
  makeplots1D( *h36, *base36, *destination36, name36);

  makeplots1D( *h41, *base41, *destination41, name41);
  makeplots1D( *h42, *base42, *destination42, name42);
  makeplots1D( *h43, *base43, *destination43, name43);
  makeplots1D( *h44, *base44, *destination44, name44);
  makeplots1D( *h45, *base45, *destination45, name45);
  makeplots1D( *h46, *base46, *destination46, name46);


// Make 2D efficiency plots

  makeplots2D( *h17, name17);
  makeplots2D( *h18, name18);
  makeplots2D( *h19, name19);
  makeplots2D( *h10, name10);

  makeplots2D( *h27, name27);
  makeplots2D( *h28, name28);
  makeplots2D( *h29, name29);
  makeplots2D( *h20, name20);

  makeplots2D( *h37, name37);
  makeplots2D( *h38, name38);
  makeplots2D( *h39, name39);
  makeplots2D( *h30, name30);

  makeplots2D( *h47, name47);
  makeplots2D( *h48, name48);
  makeplots2D( *h49, name49);
  makeplots2D( *h40, name40);
}



// Make 1D comparison plots
void makeplots1D( TH1& eff, TH1& base, TH1& destination, TString name) {

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
  g1->BayesDivide(&destination, &base, "");
  g1->GetYaxis()->SetRangeUser(0.5, 1.05);
  eff.SetLineColor(2);
  eff.SetMarkerStyle(22);
  eff.SetMarkerSize(1.4);
  eff.SetMarkerColor(2);

  //
  g1->GetYaxis()->SetTitle("Efficiency");
  if(name.Contains("_Eta"))
    g1->GetXaxis()->SetTitle("#eta");
  if(name.Contains("_Phi"))
    g1->GetXaxis()->SetTitle("#phi");
  if(name.Contains("_Pt"))
    g1->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  TCanvas canvas("canvas",name,600,600);
  g1->Draw("APE");
  eff.Draw("same");
  canvas.SaveAs(name + TString(".eps"));
  canvas.SaveAs(name + TString(".gif"));
  canvas.Close();
  delete g1;
}






// Make 2D efficiency plots
void makeplots2D( TH2& eff, TString name) {
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();

  const Int_t NRGBs = 5;
  const Int_t NCont = 200;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  tdrStyle->SetNumberContours(NCont);

  
  if(name.Contains("_Pt")) {
    eff.GetXaxis()->SetTitle("p_{T} (GeV/c)");
    eff.GetYaxis()->SetTitle("#eta     ");
  }

  if(name.Contains("_Phi")) {
    eff.GetXaxis()->SetTitle("#phi     ");
    eff.GetYaxis()->SetTitle("#eta     ");
  }

  eff.GetYaxis()->SetTitleOffset(1);


  tdrStyle->SetPadLeftMargin(0.08);
  tdrStyle->SetPadRightMargin(0.1);
 
  TCanvas canvas("canvas",name,600,600);
  eff.Draw("colz");
  gPad->Update();
  TPaletteAxis* palette = 
    (TPaletteAxis*)eff.GetListOfFunctions()->FindObject("palette");
  palette->SetLabelSize(0.02);
  canvas.SaveAs(name + TString(".eps"));
  canvas.SaveAs(name + TString(".gif"));
  canvas.Close();

}



// void Process(TH1& h, TH1& base, TH1& destination){

//   TGraphAsymmErrors *g1 = new TGraphAsymmErrors();
//   g1->BayesDivide( &destination, &base, "");
//   delete gRandom;
//   gRandom = new TRandom(111111);

//   int nBins = h.GetNbinsX();

//   for(int i=0; i<nBins; i++) {

//     Double_t x = 0.0;
//     if(h.GetBinContent(i) != 0) {
//       x = g1->Eval(h.GetBinCenter(i));
//       x += gRandom->Gaus(0.0, 0.008);
//     }
//     Double_t erx = h.GetBinError(i);

//     h.SetBinContent(i, x);
//     h.SetBinError(i, erx);
//   }

//   delete g1;
//   delete gRandom;
// }
 
