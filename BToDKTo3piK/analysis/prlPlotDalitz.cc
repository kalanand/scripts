// $Id: prlPlotDalitz.cc,v 1.1 2006/09/13 22:39:34 fwinkl Exp $
// Make PRL Dalitz plot

void prlPlotDalitz()
{
  const Int_t dalitzBINS = 30;

  // Use this times many toy events than data events for projections
  const Int_t scale = 100;

  data = read(dataTree);
  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/resutls/fitData.par");

  useYieldFitVars();
  pdfOnResDK.useDalitz();

  RooRealVar mym12(*m12);
  RooRealVar mym13(*m13);
  RooRealVar mym23(*m23);
  mym12.setBins(dalitzBINS);
  mym13.setBins(dalitzBINS);
  mym23.setBins(dalitzBINS);
  mym12.SetTitle("m^{2}(#pi^{+}#pi^{0})");
  mym13.SetTitle("m^{2}(#pi^{-}#pi^{0})");
  mym23.SetTitle("m^{2}(#pi^{+}#pi^{-})");


  RooArgList vars(mym12,mym13,mym23);

  Int_t Ntoy = scale*data->numEntries();
  cout << "Generating "<<Ntoy<<" toy events for projections..."<<endl;

  RooDataSet* genData = pdfOnResDK.generate(Ntoy);
  genData->addColumn(*s23);


  gROOT->SetStyle("BABAR");
  gStyle->SetPalette(1);
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetNdivisions(506,"XY");
  gStyle->SetOptStat(0);
  gStyle->SetHistMinimumZero();

  for (int i=0; i<vars.getSize(); i++) {
    RooRealVar& var = (RooRealVar&)vars[i];

    TCanvas* can = new TCanvas(var.GetName(),var.GetTitle(),400,400);

    TH1* hData = var.createHistogram("hData");
    TH1* hPdf = var.createHistogram("hPdf");

    hData->Sumw2();
    hPdf->Sumw2();
    data->fillHistogram(hData,RooArgList(var));
    genData->fillHistogram(hPdf,RooArgList(var));

    hPdf->Scale(1.0/scale);
  
    hData->SetTitle("");
    hData->SetMarkerStyle(8);
    hData->SetMarkerSize(0.6);

    hPdf->SetMarkerStyle(22);
    hPdf->SetMarkerColor(kBlue);
    hPdf->SetMarkerSize(0.6);
    hPdf->SetLineColor(kBlue);
   
    hData->Draw();
    hPdf->Draw("hist same");
    TString filename = "prlPlotDalitz-"+TString(var.GetName());

    can->SaveAs(filename+".eps");
    can->SaveAs(filename+".root");
  }
}
