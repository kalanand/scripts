// $Id: plotDalitzAmp.cc,v 1.2 2008/01/03 18:11:23 fwinkl Exp $
// Simply plot the Dalitz plane and all three projections from toy MC
// Run in batch mode with:
//   root -b -q -l setup.cc setupAuxPdfs.cc plotDalitzAmp.cc
//   root -b -q -l 'setup.cc(BdkPdfDKDalitz::CART) setupAuxPdfs.cc(BdkPdfDKDalitz::CART) plotDalitzAmp.cc'

void plotDalitzAmp(BdkPdfDalitzBase& pdf, Int_t nEvents = 1000,
                   TString dalitzVars = "m13:m12", TString filename = "")
{
  if (filename=="") filename = pdf.GetName();

  BdkAbsDDalitzAmp* amp(0);
  // Set Dbar amplitude to zero if it is a 'DK' pdf
  if (TString(pdf.ClassName())=="BdkPdfDKDalitz") {
    ((BdkPdfDKDalitz&)pdf).setCPparams(0,0,0);
    amp = ((BdkPdfDKDalitz&)pdf).pdfType()->dalitzAmp();
  }
  else if (TString(pdf.ClassName())=="BdkPdfDDalitz") {
    amp = ((BdkPdfDDalitz&)pdf).pdfType()->dalitzAmp();
  }

  const Int_t BINS = 50;

  gStyle->SetOptStat(0);

  if (amp) {
    cout << "Fit fractions:" << endl;
    amp->fitFractions().Print("v");
    amp->fitFractions().writeToFile(filename+"-fitFrac.txt");
  }

  cout << "Generating " << nEvents << " events from " << pdf.GetName() << endl;
  RooDataSet* data = pdf.generate(nEvents);

  // This variable is used in the s23 RooFormulaVar and needs to be set
  // to the total squared sum of the mother and all daughter masses
  BdkDalitzBase* basePdf = (BdkDalitzBase*)pdf.getPdf();
  mtotal->setVal(basePdf->M()**2 +
                 basePdf->m1()**2 +
                 basePdf->m2()**2 +
                 basePdf->m3()**2);
                 
  data->addColumn(*s23);
  
  TCanvas* c = new TCanvas("c",pdf.GetTitle(),1000,1000);
  c->Divide(2,2);

  c->cd(1);  
  TGraph* gDalitz = ((BdkDalitzBase*)pdf.getPdf())->drawBoundary(200);
  gDalitz->SetLineWidth(2);

  data->tree().Draw(dalitzVars);
  //  gDalitz->Draw("c same");

  c->cd(2);
  TH1* hm12 = new TH1D("hm12",TString("m12 ")+pdf.GetTitle(),BINS,0,3);
  data->tree().Draw("m12>>hm12","","pe");  

  c->cd(3); 
  TH1* hm13 = new TH1D("hm13",TString("m13 ")+pdf.GetTitle(),BINS,0,3);
  data->tree().Draw("m13>>hm13","","pe");

  c->cd(4);
  TH1* hm23 = new TH1D("hm23",TString("m23 ")+pdf.GetTitle(),BINS,0,3);
  data->tree().Draw("m23>>hm23","","pe");

  c->SaveAs(filename+".png");
  c->SaveAs(filename+".eps");
}

void plotDalitzAmp(Int_t nEvents = 50000)
{
  /*
  // amp *= m0*width to correct for different BW definitions
  amp = kskpAmp;
  for (int i=0; i<amp->nComps(); i++) {
    if (amp->nameRes(i)!="a0-" && amp->nameRes(i)!="Nonres")
      amp->ampRes(i)->setVal(amp->ampRes(i)->getVal()*amp->massRes(i)->getVal()*amp->gammaRes(i)->getVal());
  }
  amp = kskpBarAmp;
  for (int i=0; i<amp->nComps(); i++) {
    if (amp->nameRes(i)!="a0-" && amp->nameRes(i)!="Nonres")
      amp->ampRes(i)->setVal(amp->ampRes(i)->getVal()*amp->massRes(i)->getVal()*amp->gammaRes(i)->getVal());
  }
  */

  plotDalitzAmp(*kskpDPdf, nEvents, "m23:m12", "kskpDPdf");
  plotDalitzAmp(*kskpDbarPdf, nEvents, "m23:m12", "kskpDbarPdf");

  /*  
  plotDalitzAmp(dalitzHolderN.sigGoodD0Type(), nEvents);
  plotDalitzAmp(*ksppPdf, nEvents);
  plotDalitzAmp(*kkpPdf, nEvents);
  plotDalitzAmp(*kskkPdf, nEvents, "m12:m23");
  
  plotDalitzAmp(*kskpDPdf, nEvents, "m23:m12", "kskpDPdf");
  plotDalitzAmp(*kskpDbarPdf, nEvents, "m23:m12", "kskpDbarPdf");
  */
}

