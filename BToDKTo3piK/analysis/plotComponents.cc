void plotComponents() {

  BdkPdfDDalitz & pdf = dalitzHolderN.DpiGoodD0Type();
  
  // remove efficiency function and m23 veto
  pdf.setEfficiencyFunc(0);
  dalitzCfg->setM23VetoMass(0,0);

  BdkDDalitzAmp * amp = pdf.pdfType()->dalitzAmp();

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.2,"Y");

  TGraph* gDalitz = pdf.pdfType()->drawBoundary(200);
  gDalitz->SetLineWidth(2);

  const int nComps = amp->nComps();
  
  for (int i = 0; i < nComps; ++i){
    TString section;
    section += i;

    // Set only the amp of this resonance to 1, all others to 0:
    for (int j = 0; j < nComps; ++j){
      amp->ampRes(j)->setVal(0);
      pdf.fixAll();
    }
    amp->ampRes(i)->setVal(1);
    amp->ampRes(i)->setConstant(kFALSE);
    pdf.parametersFree().Print("V");

    int nEvents = 10000;
    if (amp->nameRes(i) == "Omega") {
      nEvents = 1000;  // fewer events for narrow omega.
    }
    
    data = pdf.generate(nEvents);
    TCanvas can(section, section, 500, 500);
    TH2 * hist = m12->createHistogram(amp->nameRes(i),*m13);
    hist->SetTitle(amp->nameRes(i));
    hist->GetXaxis()->SetTitle(m12->GetTitle());
    hist->GetYaxis()->SetTitle(m13->GetTitle());

    hist->Draw();  // axis
    gDalitz->Draw("c same");   // boundary
    data->tree().Draw("m13:m12","","same");  // data
    //data->tree().Project(hist->GetName(),"m13:m12");  // data
    //hist->Draw("col");
    //gDalitz->Draw("c same");   // boundary

    TString canName = "toy-";
    canName += amp->nameRes(i);
    canName += ".eps";
    can.SaveAs(canName);
  }
}


  
