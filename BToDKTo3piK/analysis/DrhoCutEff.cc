// Calculates the Drho cut efficiency for a given tree

double DrhoCutEff(TTree * tree) {

  TCut cutNoVeto = cutNN && cutKsVeto && cutDeltaE && cutmES && cutMD && cutBadD;

  TCut cutWithVeto = cutNoVeto && cutDtoKpi;

  cout << "Cut without Drho veto: " << endl;
  cutNoVeto.Print();
  cout << endl << "Cut with Drho veto: " << endl;
  cutWithVeto.Print();
  cout << endl;

  readCut = cutNoVeto;
  RooDataSet * dataNoVeto = read(tree);
  RooPlot * plNoVeto = mixneumass->frame();
  dataNoVeto->plotOn(plNoVeto);

  readCut = cutWithVeto;
  RooDataSet * dataWithVeto = read(tree);
  RooPlot * plWithVeto = mixneumass->frame();
  dataWithVeto->plotOn(plWithVeto);
  
  TCanvas * can = new TCanvas("can", "can", 1000, 500);
  can->Divide(2,1);
  can->cd(1);
  plNoVeto->Draw();
  can->cd(2);
  plWithVeto->Draw();


  double eff = (double)(dataWithVeto->numEntries()) / dataNoVeto->numEntries();
  cout << "eff = " << dataWithVeto->numEntries() << " / "
       << dataNoVeto->numEntries() << " = " << eff << endl;

  return eff;
}
  
