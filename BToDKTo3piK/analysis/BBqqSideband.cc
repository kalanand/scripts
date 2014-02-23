// makes a 1D N1 fit of the mES sideband

void BBqqSideband(TCut cut, TString name) {  
  readCut = cut;
//  readCut = "Deltae<-0.045&&Deltae>-0.07";
//   readCut = "Deltae>0.09";
//   readCut += cutBadD;
//  readCut = "mes<5.26";

  // Read DK
//  wght = 0.00655;
  RooDataSet * mc0 = read(chainDK);
  RooDataSet * mc  = randRep(mc0,"randAdd<0.00655");
 //  wght = 0.228;
  RooDataSet * mcBch0 = read(chainBch);
  RooDataSet * mcBch = randRep(mcBch0,"randAdd<0.228"); 

//  wght = 0.227;
  RooDataSet * mcB00 = read(chainB0);
  RooDataSet * mcB0 = randRep(mcB00,"randAdd<0.227");

//  wght = 0.629;
  RooDataSet * mcCc0 = read(chainCc);
  RooDataSet * mcCc = randRep(mcCc0,"randAdd<0.629"); 

//  wght = 0.645;
  RooDataSet * mcUds0 = read(chainUds);
  RooDataSet * mcUds = randRep(mcUds0,"randAdd<0.645");

//    cout<<" comb b "<<mc->numEntries()<<endl;                                                                                  
  cout<<" signal "<<mc->numEntries()<<endl;                                                                                  
  cout<<"bch "<<mcBch->numEntries()<<endl; 
  cout<<"b0 "<<mcB0->numEntries()<<endl; 
  cout<<"cc "<<mcCc->numEntries()<<endl; 
  cout<<"uds "<<mcUds->numEntries()<<endl; 
  mc->append(*mcBch);
  mc->append(*mcB0);
  mc->append(*mcCc);
  mc->append(*mcUds);

  cout<<mc->numEntries()<<endl; 

  RooDataSet * offData = read(chainOffData);
  data = read(chainOnData);
 
 if(TString(name)=="mes-side") { 
   mc->plotOn(mesFrame);
   mesFrame->getAttMarker()->SetMarkerColor(kGreen);
   data->plotOn(mesFrame);
   data->plotOn(mesFrame);
   mesFrame->Draw();
 }
  if(TString(name)!="mes-side") {
    mc->plotOn(DeltaeFrame);
    DeltaeFrame->getAttMarker()->SetMarkerColor(kGreen);
    data->plotOn(DeltaeFrame);
    data->plotOn(DeltaeFrame);
    DeltaeFrame->Draw();
  }

  fitOption = "mr";
  BdkPdfSum mcPdf("mcPdf", "mcPdf", *nnHolder.BBBadD0(), *nnHolder.qqBadD0());
  fit(mcPdf, *mc);
  BdkPdfSum dataPdf("dataPdf", "dataPdf", *nnHolder.BBBadD0(), *nnHolder.qqBadD0());
  fit(dataPdf, *data);

  mc->plotOn(nnFrame);
  mcPdf.getPdf()->plotOn(nnFrame);
  cout<<"chisquare "<<nnFrame->chiSquare()<<endl;
  RooArgSet argSet(*(nnHolder.qqBadD0()->getPdf()));
  mcPdf.getPdf()->plotCompOn(nnFrame, argSet);
  nnFrame->getAttLine()->SetLineColor(kRed);
  
  RooPlot * nnFrame2 = nnout->frame();
  data->plotOn(nnFrame2);
  dataPdf.getPdf()->plotOn(nnFrame2);
  dataPdf.getPdf()->plotCompOn(nnFrame2, argSet);
  nnFrame2->getAttLine()->SetLineColor(kRed);

  TCanvas * can1 = new TCanvas("can2", "canmc", 400, 400);
  nnFrame->Draw();
  can1->Print(TString(name)+"-mc-n1.eps");
  TCanvas * can2 = new TCanvas("can3", "candata", 400, 400);
  nnFrame2->Draw();
  can2->Print(TString(name)+"-data-n1.eps");
}
