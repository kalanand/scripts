// Generate shaped signal from flat signal MC
void fitShapedFlat() {

  bool doNorm = false;
  Bool_t firstCallToSetEff = kTRUE;

  RooRealVar tru_S23("tru_S23","tru_S23",0,3); //This is pi+ pi0  rho+
  RooRealVar tru_S31("tru_S31","tru_S31",0,3); //This is pi- pi0  rho-
  RooDataSet* dataFlat = RooDataSet::read("sigFlatTree.dat",
					  RooArgSet(tru_S31,tru_S23),"Q");
 

  // easy pointer to signal PDF:
  BdkPdfDDalitz pdf("pdf","pdf",tru_S23,tru_S31,1); 
  if(doNorm==true) {               
    pdf.parameters(); 
    BdkDDalitzAmp::normalizeAll(); 
    pdf.parameters().writeToFile("pars.txt");  
  }
  else { pdf.parameters().readFromFile("pars.txt"); }

             

  ofstream outFile("fitShapedFlat.dat");
  
//   remove eff func before shaping data:
  setEff(pdf, kFALSE, firstCallToSetEff);

  // shape the data:
  BdkDDalitz* thePdf =  (BdkDDalitz*)(pdf.getPdf());


  double max = 0;  
  double avg = 0;

  int numEnt = dataFlat->numEntries();

  RooArgSet * pdfset = thePdf->getParameters(RooArgSet());

  // find max:
  for (int e = 0; e < numEnt; ++e){
     RooArgSet * event = dataFlat->get(e);
     (*pdfset) = (*event);
    double val = thePdf->getVal();
    avg += val;
    if (val > max) max = val;
  }

  avg /= numEnt;
  cout << "average PDF value = " << avg <<". max PDF value = " << max 
       << endl;


  // init a RDS with the same vars as flat:
  RooDataSet * dSet = new RooDataSet(TString(dataFlat->GetName()) 
				     + ".shaped",
				       TString(dataFlat->GetTitle()) 
				     + ".shaped",
				       RooArgSet(tru_S23,tru_S31));

  // accept/reject:
  for (int e = 0; e < numEnt; ++e){
    RooArgSet * event = dataFlat->get(e);
    (*pdfset) = (*event);
    double val = thePdf->getVal(event);
    double rand = max * RooRandom::uniform();
    if (val > rand) dSet->add(*event);
  }




  // reset eff func before fitting:
  setEff(pdf, kFALSE, firstCallToSetEff);


  RooDataSet * gSet = pdf.generate(dSet->numEntries());
	
  TH2 * hFlat = dataFlat->createHistogram(tru_S23, tru_S31);
  TH2 * hDSet = dSet->createHistogram(tru_S23, tru_S31);
  TH2 * hGSet = gSet->createHistogram(tru_S23, tru_S31);

  char * canName = "Flat Sig MC";	
  hFlat->SetTitle(canName);
  hDSet->SetTitle("Shaped by sig PDF");
  hGSet->SetTitle("Generated from sig PDF, eff=1");
	
  TCanvas * can = new TCanvas(canName, canName, 1185, 435);
  gStyle->SetOptStat(0);
  can->Divide(3,1);
  can->cd(1);
  hFlat->Draw("colz");
  can->cd(2);
  hDSet->Draw("colz");
  can->cd(3);
  hGSet->Draw("colz");  
  can->SaveAs("fitShapedFlat.gif");
  can->SaveAs("fitShapedFlat.eps");
  can->Close();
  delete can;
  // Plot chi^2
  plotPull(*hGSet, *hDSet);



  // write the shaped dataset:
  TString fileName("sigShapedFlat.root");
  cout << "Writing shaped data to " << fileName << endl;	
  TFile file(fileName, "recreate");
  dSet->tree().Write();
  file.Close();


  RooArgSet allPars = pdf.parameters();
  RooRealVar* var1 = allPars.find("pdf.pdf.dalitzAmp.Nonres_amp");
  RooRealVar* var2 = allPars.find("pdf.pdf.dalitzAmp.Nonres_phase");
  RooRealVar* var3 = allPars.find("pdf.pdf.dalitzAmp.Rho+_amp");
  RooRealVar* var4 = allPars.find("pdf.pdf.dalitzAmp.Rho+_phase");
  RooRealVar* var5 = allPars.find("pdf.pdf.dalitzAmp.Rho-_amp");
  RooRealVar* var6 = allPars.find("pdf.pdf.dalitzAmp.Rho-_phase");
  RooRealVar* var7 = allPars.find("pdf.pdf.dalitzAmp.Rho0_amp");
  RooRealVar* var8 = allPars.find("pdf.pdf.dalitzAmp.Rho0_phase");
  var1->setConstant(false);
  var2->setConstant(false);
  var5->setConstant(false);
  var6->setConstant(false);
  var7->setConstant(false);
  var8->setConstant(false);
  var3->setVal(1.0);
  var3->setConstant(true);
  var4->setVal(0.0);
  var4->setConstant(true);

  RooFitResult *fitRes = thePdf->fitTo(*dSet,"mrt");
}










// apply or remove the efficiency function, changing the normalization
// coefficients appropriately:
void setEff(BdkPdfDDalitz & wrapper, Bool_t applyEff, Bool_t & firstCall) {
  BdkDDalitz * pdf = (BdkDDalitz*)(wrapper.getPdf());
  
  if (applyEff) {
    pdf->setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());    
    wrapper.parameters().readFromFile("../BToDKTo3piK/params/dalitzNorm.par");
    cout << "---- setting efficiency function" << endl;
  }
  else {
    // If the PDF has an efficiency function, remove it:
    pdf->setEfficiencyFunc(0);
    
    // On the first call, normalize. Subsequently read from file:
    if (kTRUE == firstCall) {
      pdf->dalitzAmp()->calNorm();
      wrapper.parameters().writeToFile("fitShapedFlat-noEffNorm.par");
      firstCall = kFALSE;
    }
    else {
      wrapper.parameters().readFromFile("fitShapedFlat-noEffNorm.par");
    }
    cout << "---- removing efficiency function, params=" << endl;
  }
}








void plotPull(TH2& hGSet, TH2& hDSet) {

  int BINS1 = hGSet.GetNbinsX();
  int BINS2 = hGSet.GetNbinsY();
  TAxis* xaxis = hGSet.GetXaxis();
  TAxis* yaxis = hGSet.GetYaxis();
  Double_t MAX1 = xaxis->GetBinUpEdge(xaxis->GetLast());
  Double_t MAX2 = yaxis->GetBinUpEdge(yaxis->GetLast());
  Double_t MIN1 = xaxis->GetBinLowEdge(xaxis->GetFirst());
  Double_t MIN2 = yaxis->GetBinLowEdge(yaxis->GetFirst());
  
  TH2D hchi2("hchi2","",BINS1,MIN1,MAX1,BINS2,MIN2,MAX2);  
  TH1D pullhist("pullhist","",25,-4,4); 
  int ndof = hchi2.GetNbinsX() * hchi2.GetNbinsY();
  ndof--;   // histos are normalized to each other
  TCanvas bCanvas("bCanvas", "", 800, 600);

  Double_t chi2Total = 0;    // total chi2
  for (int x=1; x <= hchi2.GetNbinsX(); x++) {
    for (int y=1; y <= hchi2.GetNbinsY(); y++) {
      
      Double_t bin1 =  hDSet.GetBinContent(x,y);
      Double_t bin2 = hGSet.GetBinContent(x,y); 
      
      if (bin1==0 || bin2==0) ndof--;    // no data -> one less dof
      else {
	Double_t sqrtChi2 = (bin1-bin2);
	sqrtChi2 /= sqrt(hDSet.GetBinError(x,y)**2 
			 + hGSet.GetBinError(x,y)**2);        
	hchi2.SetBinContent(x,y,sqrtChi2);
	chi2Total += sqrtChi2**2;
	pullhist.Fill(sqrtChi2);
      }
    }
  }
  hchi2.SetMaximum(4);
  hchi2.SetMinimum(-4);
  hchi2.SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  pullhist.SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  pullhist.Fit("gaus");
  gStyle->SetOptStat(0);
  bCanvas.Clear();
  bCanvas.Divide(2,2);
  bCanvas.cd(1);
  hDSet.Draw("colz");
  bCanvas.cd(2);
  hGSet.Draw("colz");
  bCanvas.cd(3);
  hchi2.Draw("colz");
  bCanvas.cd(4);
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  pullhist.Draw();
  bCanvas.Update();
  bCanvas.SaveAs("fitShapedFlat_2.gif");
  bCanvas.SaveAs("fitShapedFlat_2.eps");
  bCanvas.Close();
}
