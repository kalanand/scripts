#include <iostream>
#include <vector>
#include <string>

using namespace RooFit;


void EleEnergyScaleLH()
{

  TCut isBB("probe_gsfEle_isEB && tag_gsfEle_isEB && tag_gsfEle_ecalDrivenSeed && probe_gsfEle_ecalDrivenSeed");
  TCut isEB("((probe_gsfEle_isEB && tag_gsfEle_isEE)||(probe_gsfEle_isEE && tag_gsfEle_isEB)) && tag_gsfEle_ecalDrivenSeed && probe_gsfEle_ecalDrivenSeed");

   TCut cutBB("probe_passing && probe_gsfEle_isEB && tag_gsfEle_isEB && tag_gsfEle_ecalDrivenSeed && probe_gsfEle_ecalDrivenSeed");
   TCut cutEB("probe_passing && ((probe_gsfEle_isEB && tag_gsfEle_isEE)||(probe_gsfEle_isEE && tag_gsfEle_isEB)) && tag_gsfEle_ecalDrivenSeed && probe_gsfEle_ecalDrivenSeed");

   const int nBINS = 40;

  TFile* fData = new TFile("allTPtrees_2900nb.root");
  fData->cd("IdToHLT");
  TTree* tpTree = (TTree*) fData->Get("fitter_tree");
  TH1F* data_BB = new TH1F("data_BB","data_BB", nBINS, 60, 120);
  TH1F* data_EB = new TH1F("data_EB","data_EB", nBINS, 60, 120);

  fitter_tree->Draw( "mass>>data_EB",cutEB,"goff");
  fitter_tree->Draw( "mass>>data_BB",cutBB,"goff");

  data_BB->Scale(0.5); // because we have twice as many entries
  data_EB->Scale(0.5); // because we have twice as many entries



  // The fit variable - lepton invariant mass
  RooRealVar Mass("Mass","m_{ee}",60.0, 120.0, "GeV/c^{2}");


  // Make the category variable that defines the two fits,
  // namely whether the probe passes or fails the eff criteria.
  RooCategory sample("sample","") ;
  sample.defineType("BB", 1) ;
  sample.defineType("BE", 2) ; 

  gROOT->cd();
  ///////// convert Histograms into RooDataHists
  RooDataHist* dh_BB = new RooDataHist("dh_BB","dh_BB",
					  RooArgList(Mass), data_BB);
  RooDataHist* dh_EB = new RooDataHist("dh_EB","dh_EB",
					  RooArgList(Mass), data_EB);

  RooDataHist* data = new RooDataHist( "fitData","fitData",
				       RooArgList(Mass),Index(sample),
				       Import("BB",*dh_BB), Import("BE",*dh_EB) ); 

  data->get()->Print();
  cout << "Made datahist" << endl;



  TFile* fMC = new TFile("allTPtrees_mc.root");
  fMC->cd("IdToHLT");
  TH1F* mc_BB = new TH1F("mc_BB","mc_BB", 1000, 60, 120);
  TH1F* mc_EB = new TH1F("mc_EB","mc_EB", 500, 60, 120);
  fitter_tree->Draw("mass>>mc_EB", isEB,"goff");
  fitter_tree->Draw("mass>>mc_BB", isBB,"goff");
  mc_BB->SetLineColor(4);
  mc_BB->SetLineWidth(3);
  TH1F* mc_BB2 = mc_BB->Clone("mc_BB2");
  mc_BB2->SetLineColor(2);
  mc_BB2->SetLineWidth(3);
  mc_BB2->SetLineStyle(2);

  RooRealVar massShiftBB("massShiftBB","",-1.0215, -3., 0.);
  RooRealVar massShiftBE("massShiftBE","",-1.7759, -3., 0.);
  RooFormulaVar shiftedMassBB("shiftedMassBB", "@0-@1", RooArgSet(Mass, massShiftBB) );
  RooFormulaVar shiftedMassBE("shiftedMassBE","@0-@1", RooArgSet(Mass, massShiftBE) );
  RooDataHist rdh_BB("rdh_BB","", Mass, mc_BB);
  RooDataHist rdh_BE("rdh_BE","", Mass, mc_EB);

  RooRealVar zero("zero","",0.0);
  RooRealVar resoBB("resoBB","",  1.64152, 1.5, 2.5);
  RooGaussModel resModelBB("resModelBB","gaussian resolution model", 
                Mass, zero, resoBB);
  RooHistPdf shapeBB("shapeBB", "", shiftedMassBB, Mass, rdh_BB);
  RooFFTConvPdf pdfBB("pdfBB","", shiftedMassBB, Mass, shapeBB, resModelBB);
  RooRealVar resoBE("resoBE","",  2.2118, 1.5, 2.5);
  RooGaussModel resModelBE("resModelBE","gaussian resolution model", 
                Mass, zero, resoBE);
  RooHistPdf shapeBE("shapeBE", "", shiftedMassBE, Mass, rdh_BE);
  RooFFTConvPdf pdfBE("pdfBE","", shiftedMassBE, Mass, shapeBE, resModelBE);


  //RooHistPdf pdfBB("pdfBB", "", shiftedMassBB, Mass, rdh_BB);
//  RooHistPdf pdfBE("pdfBE", "", shiftedMassBE, Mass, rdh_BE);
  RooHistPdf pdfBBunc("pdfBBunc", "", Mass, rdh_BB);
  RooHistPdf pdfBEunc("pdfBEunc", "", Mass, rdh_BE);


   // The total simultaneous fit ...
   RooSimultaneous totalPdf("totalPdf","totalPdf", sample);
   totalPdf.addPdf(pdfBB,"BB");
   totalPdf.Print();
   totalPdf.addPdf(pdfBE,"BE");
   totalPdf.Print();


   RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
   cfg->setEpsAbs(1E-6);
   cfg->setEpsRel(1E-6);
   totalPdf.setIntegratorConfig(*cfg);


  // ********* Do the Actual Fit ********** //  
   RooFitResult *fitResult = totalPdf.fitTo(*data, Save(true), 
					    PrintEvalErrors(-1),Warnings(false) );
  fitResult->Print("v");

  double m0B = mc_BB->GetMean();
  double shiftB = massShiftBB.getVal();
  cout << "moB = " << m0B << endl;
  double es_B  = m0B/(m0B+shiftB);
  double es_B_err  = es_B*es_B*data_BB->GetMeanError()/m0B;

  double m0EB = mc_EB->GetMean();
  double shiftEB = massShiftBE.getVal();
  cout << "moEB = " << m0EB << endl;
  double es_E  = pow(m0EB/(m0EB+shiftEB),2)/es_B;
  double es_E_err  = sqrt( pow(es_E*es_E*data_EB->GetMeanError()/m0EB/es_E, 2) 
  + pow(es_B_err/es_B, 2) ) *es_E;


  // ********** Make and save Canvas for the plots ********** //
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.19);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);
  RooAbsData::ErrorType errorType = RooAbsData::SumW2;
  char title[50];

  TString cname = "Zmass_shift_BB";
  TCanvas* c = new TCanvas(cname, cname,500,500);
  RooPlot* frame1 = Mass.frame();
  dh_BB->plotOn(frame1,RooFit::DataError(errorType));
  pdfBBunc.plotOn(frame1,ProjWData(*dh_BB),LineStyle(kDashed),LineColor(kRed));
  pdfBB.plotOn(frame1,ProjWData(*dh_BB));
  frame1->SetMinimum(0);
  frame1->Draw("e0");
  sprintf(title, "#chi^{2}/#nu = %.2f", frame1->chiSquare());
  TLatex tex;
  tex.SetNDC();
  tex.SetTextAlign(12);
  tex.SetTextSize(0.04);
  tex.DrawLatex(0.26,0.88, title);
  TLegend* legend = new TLegend(0.7,0.75,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetHeader( "EB + EB");
  legend->AddEntry( data_BB, "Data", "PLE");  
  legend->AddEntry( mc_BB2, "MC", "L");
  legend->AddEntry( mc_BB, "MC shifted", "L");
  legend->Draw();

   TPaveText *plotlabel1 = new TPaveText(0.6,0.67,0.8,0.72,"NDC");
   plotlabel1->SetTextColor(kBlack);
   plotlabel1->SetFillColor(kWhite);
   plotlabel1->SetBorderSize(0);
   plotlabel1->SetTextAlign(12);
   plotlabel1->SetTextSize(0.03);
   sprintf(title, "#Delta_{EB} = %.4f #pm %.4f", es_B, es_B_err);
   plotlabel1->AddText(title);
   plotlabel1->Draw();
   TPaveText *plotlabel2 = new TPaveText(0.6,0.62,0.8,0.67,"NDC");
   plotlabel2->SetTextColor(kBlack);
   plotlabel2->SetFillColor(kWhite);
   plotlabel2->SetBorderSize(0);
   plotlabel2->SetTextAlign(12);
   plotlabel2->SetTextSize(0.03);
   sprintf(title, "#Delta_{EE} = %.4f #pm %.4f", es_E, es_E_err);
   plotlabel2->AddText(title);
   plotlabel2->Draw();

  c->SaveAs( cname + TString(".eps"));
  c->SaveAs( cname + TString(".gif"));
  c->SaveAs( cname + TString(".root"));
  c->SaveAs( cname + TString(".png"));
  c->SaveAs( cname + TString(".C"));


  cname = "Zmass_shift_EB";
  TCanvas* c = new TCanvas(cname, cname,500,500);
  RooPlot* frame2 = Mass.frame();
  dh_EB->plotOn(frame2,RooFit::DataError(errorType));
  pdfBEunc.plotOn(frame2,ProjWData(*dh_EB),LineStyle(kDashed),LineColor(kRed));
  pdfBE.plotOn(frame2,ProjWData(*dh_EB));
  frame2->SetMinimum(0);
  frame2->Draw("e0");
  sprintf(title, "#chi^{2}/#nu = %.2f", frame2->chiSquare());
  TLatex tex;
  tex.SetNDC();
  tex.SetTextAlign(12);
  tex.SetTextSize(0.04);
  tex.DrawLatex(0.26,0.88, title);

  TLegend* legend = new TLegend(0.7,0.75,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetHeader( "EB + EE");
  legend->AddEntry( data_BB, "Data", "PLE");  
  legend->AddEntry( mc_BB2, "MC", "L");
  legend->AddEntry( mc_BB, "MC shifted", "L");
  legend->Draw();

   TPaveText *plotlabel1 = new TPaveText(0.6,0.67,0.8,0.72,"NDC");
   plotlabel1->SetTextColor(kBlack);
   plotlabel1->SetFillColor(kWhite);
   plotlabel1->SetBorderSize(0);
   plotlabel1->SetTextAlign(12);
   plotlabel1->SetTextSize(0.03);
   sprintf(title, "#Delta_{EB} = %.4f #pm %.4f", es_B, es_B_err);
   plotlabel1->AddText(title);
   plotlabel1->Draw();
   TPaveText *plotlabel2 = new TPaveText(0.6,0.62,0.8,0.67,"NDC");
   plotlabel2->SetTextColor(kBlack);
   plotlabel2->SetFillColor(kWhite);
   plotlabel2->SetBorderSize(0);
   plotlabel2->SetTextAlign(12);
   plotlabel2->SetTextSize(0.03);
   sprintf(title, "#Delta_{EE} = %.4f #pm %.4f", es_E, es_E_err);
   plotlabel2->AddText(title);
   plotlabel2->Draw();

  c->SaveAs( cname + TString(".eps"));
  c->SaveAs( cname + TString(".gif"));
  c->SaveAs( cname + TString(".root"));
  c->SaveAs( cname + TString(".png"));
  c->SaveAs( cname + TString(".C"));

   //mZee->SaveAs("Zee_mass.eps");

}
