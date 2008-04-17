
TH1* diffDataMC(RooRealVar* var, TCut cut)
{
  RooRealVar* treeVar = var;
  if (var==qprime) treeVar = nnout;
  else if (var==dprime) treeVar = bknnout;

  // Create weighted MC histogram
  TH1* hTemp = treeVar->createHistogram("hTemp");
  weightedMCHisto(hTemp,treeVar->GetName(),cut,false);
  
  // Add transformed NN variables
  RooDataHist histMC("histMC","histMC",RooArgList(*treeVar),hTemp);
  histMC.addColumns(RooArgList(*dprimeF,*qprimeF));

  readCut = cut;
  data = read(dataTree);
  
  //  TH1* hData = hMC->Clone("hData");
  //  dataTree->Project(hData->GetName(),var->GetName(),cut);

  RooDataHist histData("histData","histData",RooArgSet(*var),*data);

  TH1* hMC = histMC.createHistogram("hMC",*var);
  TH1* hData = histData.createHistogram("hData",*var);

  // Normalize MC to data
  hMC->Sumw2();
  hMC->Scale(hData->Integral()/hMC->Integral());

  TH1* hDiff = hData->Clone("hDiff");
  hDiff->Sumw2();


  hDiff->Add(hMC,-1);

  return hDiff;
}


TH2D* diffDataMC2D(TCut cut)
{
  TString var = "d0ppmupmass**2:d0pppupmass**2"; 
  TH2D* hMC = (TH2D*)m12->createHistogram("hMC",YVar(*m13));
  TH2D* hData = (TH2D*)m12->createHistogram("hData",YVar(*m13));
  
  weightedMCHisto(hMC,var,cut,false);
  //  hMC->Draw();

  
  dataTree->Project(hData->GetName(),var,cut);

  hData->Sumw2();
  hData->Scale(1/hData->Integral());
  hMC->Sumw2();
  hMC->Scale(1/hMC->Integral());

  TH2D* hDiff = (TH2D*)hData->Clone("hDiff");

  hDiff->Add(hMC,-1);

  return hDiff;
}
