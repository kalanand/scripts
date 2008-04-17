// studies the effect of a cut on mixneumass to remove the DpiX
// peaking background

void mixNeuMassStudy(TCut & typeCut = cutDPiX){
  mes->setMin(5.2);
  
  TCut basic("(-0.07<Deltae&&Deltae<0.06&&1.83<d0mass&&d0mass<1.895&&nnout>0.1&&bknnout>0.25)");

  RooDataSet * mcBB = (RooDataSet *)dataSetFromTree(bbTree, basic);
  RooPlot * mesP1 = mes->frame();
  mesP1->SetTitle("BB MC");
  mcBB->plotOn(mesP1);
  
  TCut cut2(basic);
  cut2 += "(mixneumass<1.84||mixneumass>1.89)";
  RooDataSet * mcBB2 = (RooDataSet *)dataSetFromTree(bbTree, cut2);
  RooPlot * mesP2 = mes->frame();
  mesP2->SetTitle("BB MC, m(Kpi) cut");
  mcBB2->plotOn(mesP2);
  
  
  TCut cut3(basic);
  cut3 += typeCut;
  RooDataSet * mcDpiX1 = (RooDataSet *)dataSetFromTree(bbTree, cut3);
  RooPlot * mesP3 = mes->frame();
  mesP3->SetTitle("DpiX MC");
  mcDpiX1->plotOn(mesP3);
  
  TCut cut4(cut3);
  cut4 += "(mixneumass<1.84||mixneumass>1.89)";
  RooDataSet * mcDpiX2 = (RooDataSet *)dataSetFromTree(bbTree, cut4);
  RooPlot * mesP4 = mes->frame();
  mesP4->SetTitle("DpiX MC, m(Kpi) cut");
  mcDpiX2->plotOn(mesP4);
  
  
  TCanvas * can3 = new TCanvas("can3", "can3", 1000, 500);
  can3->Divide(2,1);

  RooPlot * mesPBB = mes->frame(30);
  mesPBB->SetTitle("BB MC");
  mcBB->plotOn(mesPBB);
  mesPBB->getAttMarker()->SetMarkerColor(kBlue);
  mesPBB->getAttLine()->SetLineColor(kBlue);
  mcBB2->plotOn(mesPBB);
  mesPBB->getAttMarker()->SetMarkerColor(kRed);
  mesPBB->getAttLine()->SetLineColor(kRed);
  can3->cd(1);
  mesPBB->Draw();

  RooPlot * mesPDpiX = mes->frame(30);
  mesPDpiX->SetTitle("DpiX MC");
  mcDpiX1->plotOn(mesPDpiX, "L");
  mesPDpiX->getAttMarker()->SetMarkerColor(kBlue);
  mesPDpiX->getAttLine()->SetLineColor(kBlue);
  mcDpiX2->plotOn(mesPDpiX, "L");
  mesPDpiX->getAttMarker()->SetMarkerColor(kRed);
  mesPDpiX->getAttLine()->SetLineColor(kRed);
  can3->cd(2);
  mesPDpiX->Draw();
  
  
  TCanvas * can = new TCanvas("can", "can", 1000, 1000);
  can->Divide(2,2);
  can->cd(1);
  mesP1->Draw();
  can->cd(2);
  mesP2->Draw();
  can->cd(3);
  mesP3->Draw();
  can->cd(4);
  mesP4->Draw();
  
  
  TCut truth = basic;
  truth += "excltruth==16";
  RooDataSet * mcSig = (RooDataSet *)dataSetFromTree(sigTree, truth);
  
  TCut truth2 = truth;
  truth2 += "(mixneumass<1.8395||mixneumass>1.8895)";
  RooDataSet * mcSig2 = (RooDataSet*)dataSetFromTree(sigTree, truth2);
  
  cout << "Signal efficiency of mixneumass cut = " << mcSig2->numEntries() 
       << " / " << mcSig->numEntries() << endl;
  
  TCanvas * can2 = new TCanvas("can2", "can2", 1000, 1000);
  can2->Divide(2,2);

  can2->cd(1);
  TH2 * bbDalitzAll = new TH2D("bbDalitzAll","Dalitz all",100,0,3,100,0,3);
  RooDataSet * mcDpiX1Minus = chargeCut(mcDpiX1, -1);
  mcDpiX1Minus->tree().Draw("d0pppupmass**2:d0ppmupmass**2>>bbDalitzAll");


  can2->cd(2);
  TList * mcDpiX1List = mcDpiX1->split(*Hdtrkchge);
  TH2 * bbDalitzCut = new TH2D("bbDalitzCut","Dalitz wcut",100,0,3,100,0,3);
  RooDataSet * mcDpiX2Minus = chargeCut(mcDpiX2, -1);
  mcDpiX2Minus->tree().Draw("d0pppupmass**2:d0ppmupmass**2>>bbDalitzCut");


  TH2 * mdeall = mcDpiX1->createHistogram(*Deltae, *mes);
  mdeall->SetTitle("DpiX MC");
  TH2 * mdecut = mcDpiX2->createHistogram(*Deltae, *mes);
  mdecut->SetTitle("DpiX MC, m(Kpi) cut");


  can2->cd(3);
  mdeall->Draw();
  can2->cd(4);
  mdecut->Draw();
  
}

