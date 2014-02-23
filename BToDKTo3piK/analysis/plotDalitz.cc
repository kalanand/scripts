// $Id: plotDalitz.cc,v 1.2 2006/09/13 22:39:33 fwinkl Exp $
// Script to plot the three Dalitz projections

void plotDalitz(RooRealVar* var)   // m12, m13 or m23
{
  if (var==0) return;

  const Int_t BINS = 40;

  RooBinning bins;
  bins.setMin(0);
  bins.setMax(3);
  bins.addUniform(BINS,0,3);   // uniform binning


  // Add one bin for veto window
  if (var==m23) {
    bins.addBoundary(dalitzCfg->m23VetoMin());
    bins.addBoundary(dalitzCfg->m23VetoMax());
    setupPdfHolder(BdkPdfDKDalitz::POLAR, m12, s13_23);
  }

  // Set binning
  var->setBinning(bins);

  // Read final fit result
  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/results/fitData.par");
  
  // Plot data
  data = read(dataTree);
  //  data->addColumn(*s23);

  RooPlot* p = var->frame();
  data->plotOn(p);

  // Dataset with categories
  RooDataSet catData("catData", "catData", RooArgSet(*Hdtrkchge));
  Hdtrkchge->setIndex(-1);
  catData.add(*Hdtrkchge);
  Hdtrkchge->setIndex(1);
  catData.add(*Hdtrkchge);

  // Plot PDF
  RooAbsPdf* pdf = pdfOnResDK.getPdf();

  plot1dIntCfg.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  plot1dIntCfg.setEpsAbs(1e-30);
  plot1dIntCfg.setEpsRel(1e-2);
  pdf->setIntegratorConfig(plot1dIntCfg);

  RooDataSet* projData = pdfOnResDK.generate();
  m13->setBins(20);
  RooDataSet projDataB("projDataB","",projData,RooArgSet(*m13,*Hdtrkchge));
  pdf->plotOn(p,ProjWData(projDataB));

  
  TCanvas* can = new TCanvas("can",var->GetTitle(),600,600);
  can->SetTopMargin(0.05);
  can->SetLeftMargin(0.12);
  can->SetRightMargin(0.05);
  gStyle->SetTitleOffset(1.4,"Y");

  p->SetTitle("");

  // Make Y axis label (RooFit labels get screwed up by veto window bin)
  Double_t gev = 3.0 / BINS;
  TString ylabel = "Events / (";
  ylabel.Form("Events / (%g GeV^{2}/c^{4})",gev);
  p->GetYaxis()->SetTitle(ylabel);
  p->Draw();

  TString filename = "plotDalitz-";
  filename += var->GetName();
  can->SaveAs(filename+".eps");
  can->SaveAs(filename+".root");
}
