void TwoDFit(RooRealVar *var1, RooRealVar *var2,  RooAbsData * dat) {
//check whether var1 is the same as var2

  if(var1 == var2 ) {
    cout<<" not good, should be two different input RooRealVar variables";
    return ;
  }

  RooProdPdf * pdfVarProd[BdkEvtTypes::NTYPES];

  RooExtendPdf * extPdfVar[BdkEvtTypes::NTYPES];

  int varFlag1 = -1;
  int varFlag2 = -1;
  RooPlot * frame1 = 0;
  RooPlot * frame2 = 0;
  if (var1 == mes) {
    varFlag1 = BdkPdfProdAll::MES;
    frame1 = mesFrame;
  }
  else if (var1 == Deltae) {
    varFlag1 = BdkPdfProdAll::DELTAE;
    frame1 = DeltaeFrame;
  }
  else if (var1 == d0mass) {
    varFlag1 = BdkPdfProdAll::MD;
    frame1 = d0massFrame;
  }

   if (var2 == mes) {
    varFlag2 = BdkPdfProdAll::MES;
    frame2 = mesFrame;
  }
  else if (var2 == Deltae) {
    varFlag2 = BdkPdfProdAll::DELTAE;
    frame2 = DeltaeFrame;
  }
  else if (var2 == d0mass) {
    varFlag2 = BdkPdfProdAll::MD;
    frame2 = d0massFrame;
  }

 RooExtendPdf * extPdfVar[BdkEvtTypes::NTYPES];
 RooAbsReal * normVar[BdkEvtTypes::NTYPES];
 RooArgSet extPdfList;
 RooArgSet evtNums;

 RooAddPdf * addPdf=0; 

 BdkPdfProdAll * ProdType[BdkEvtTypes::NTYPES];

 for (int t = 0; t < BdkEvtTypes::NTYPES; ++t) {
    normVar[t] = pdfOnResDK.numEvt(t);
    evtNums.add(*(pdfOnResDK.numEvt(t)));
    ProdType[t] = (BdkPdfProdAll *) pdfOnResDK.pdf(t);
    pdfVarProd[t] = new RooProdPdf( TString(ProdType[t]->GetName())+TString(var1->GetName())
                                        +TString(var2->GetName()) + ".2dProdPdf",
				    TString(ProdType[t]->GetTitle()) +TString(var1->GetTitle())
                                        +TString(var2->GetTitle())+ " 2dProdPdf",
				    *(ProdType[t]->getVarPdf(varFlag1)->getPdf()),
				    *(ProdType[t]->getVarPdf(varFlag2)->getPdf()) );

    extPdfVar[t] = new RooExtendPdf(TString(pdfVarProd[t]->GetName()) + ".2dextPdf",
				    TString(pdfVarProd[t]->GetTitle()) + " 2dextPdf",	
                                    *(pdfVarProd[t]),
				    *normVar[t] );
    extPdfList.add(*extPdfVar[t]);
  }

  addPdf = new RooAddPdf( TString(pdfOnResDK.GetName()) + TString(var1->GetName()) 
				+ TString(var2->GetName())  + ".sumPdf",
                          TString(pdfOnResDK.GetTitle())+ TString(var1->GetTitle())
			        + TString(var2->GetTitle()) + " sumPdf",
                          extPdfList, evtNums);

  fitOption = "mer";
  RooFitResult * result = 0; 
  if(doFit) {
       result =  addPdf->fitTo(*dat,fitOption);
  }
  if(0 != result ) {	
	result->Print("V");
  }	

  dat->plotOn(frame1);
  addPdf->plotOn(frame1);
  dat->plotOn(frame2);
  addPdf->plotOn(frame2);
  frame1->getAttLine()->SetLineColor(kRed);
  frame2->getAttLine()->SetLineColor(kRed);
  cout<<" chisqurare of "<< var1->GetName() << " plot: " << frame1->chiSquare() <<endl;
  cout<<" chisqurare of "<< var2->GetName() << " plot: " << frame2->chiSquare() <<endl;
  RooArgSet CompSet(*(extPdfVar[9]),*(extPdfVar[8]));
  addPdf->plotCompOn(frame1, CompSet);
  addPdf->plotCompOn(frame2, CompSet);
  frame1->getAttLine()->SetLineColor(kBlue);
  frame2->getAttLine()->SetLineColor(kBlue);
  CompSet.add(*(extPdfVar[7]));
  CompSet.add(*(extPdfVar[6]));
  CompSet.add(*(extPdfVar[5]));
  CompSet.add(*(extPdfVar[4]));
  CompSet.add(*(extPdfVar[3]));
  CompSet.add(*(extPdfVar[2]));
  addPdf->plotCompOn(frame1, CompSet);
  addPdf->plotCompOn(frame2, CompSet);
  frame1->getAttLine()->SetLineColor(kGreen);
  frame2->getAttLine()->SetLineColor(kGreen);
  CompSet.add(*(extPdfVar[0]));
  addPdf->plotCompOn(frame1, CompSet);
  addPdf->plotCompOn(frame2, CompSet);
  frame1->getAttLine()->SetLineColor(kCyan);
  frame2->getAttLine()->SetLineColor(kCyan);


  TCanvas * Cann1 = new TCanvas(TString(var1->GetName()), TString(var1->getTitle()), 600, 600);
  frame1->Draw();
  TCanvas * Cann2 = new TCanvas(TString(var2->GetName()), TString(var2->getTitle()), 600, 600);
  frame2->Draw();
} 

