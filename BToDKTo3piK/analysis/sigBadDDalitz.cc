// Study the Dalitz plot of sigBadD:

void sigBadDDalitz() {
  mass12mc = new RooRealVar("d0pppmcmass", "M(#pi^{+}#pi^{0})", 0, 2, "GeV/c^{2}");
  mass13mc = new RooRealVar("d0ppmmcmass", "M(#pi^{-}#pi^{0})", 0, 2, "GeV/c^{2}");

  // RooFormularVar's to calculate squared Dalitz masses
  s12mc = new RooFormulaVar("m12mc","M(#pi^{+}#pi^{0})^{2}","@0*@0",*mass12mc);
  s13mc = new RooFormulaVar("m13mc","M(#pi^{-}#pi^{0})^{2}","@0*@0",*mass13mc);

  // Can be used to access squared Dalitz masses in a RooDataset
  m12mc = new RooRealVar("m12mc","M(#pi^{+}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");
  m13mc = new RooRealVar("m13mc","M(#pi^{-}#pi^{0})^{2}", 0, 3, "GeV^{2}/c^{4}");


  allVars.add(*mass12mc);
  allVars.add(*mass13mc);

  readCut = cutSigReg + cutDKBadD;
  data = read(sigTree);

  readCut = cutSigReg;
  RooDataSet * dataG = (RooDataSet*)read(sigTree);

  ((RooDataSet*)data)->addColumn(*s12mc);
  ((RooDataSet*)data)->addColumn(*s13mc);

  ((RooDataSet*)dataG)->addColumn(*s12mc);
  ((RooDataSet*)dataG)->addColumn(*s13mc);


  TCanvas * can = new TCanvas("can", "can", 1200, 1200);
  can->Divide(2,2);
  can->cd(1);
  ((RooDataSet*)data)->tree()->Draw("m12mc-m12:m12");
  can->cd(2);
  ((RooDataSet*)dataG)->tree()->Draw("m12mc-m12:m12");
  can->cd(3);
  ((RooDataSet*)data)->tree()->Draw("m12mc-m12", "abs(m12mc-m12)<0.5");
  can->cd(4);
  ((RooDataSet*)dataG)->tree()->Draw("m12mc-m12", "abs(m12mc-m12)<0.5");

  /*

  TH2F * his1 = new TH2F("hist1", "hist1", 100, 0, 3, 100, 0, 3);
  ((RooDataSet*)dataG)->tree()->Project("hist1", "m12:m13", "abs(m12mc-m12)>1");
  hist1->SetMarkerStyle(6);
  hist1->SetMarkerSize(3);
  hist1->Draw();

  TH2F * his2 = new TH2F("hist2", "hist2", 100, 0, 3, 100, 0, 3);
  ((RooDataSet*)dataG)->tree()->Project("hist2", "m12mc:m13mc", "abs(m12mc-m12)>1");
  hist2->SetMarkerStyle(6);
  hist2->SetMarkerSize(3);
  hist2->SetMarkerColor(4);
  hist2->Draw("same");
  */
}
