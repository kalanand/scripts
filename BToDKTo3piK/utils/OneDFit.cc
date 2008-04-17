void OneDFit(RooRealVar *var, RooAbsData * dat=0 ) {

  RooExtendPdf * extPdfVar[BdkEvtTypes::NTYPES];

  int varFlag = -1;
  RooPlot * frame = 0;
  if (var == mes) {
    varFlag = BdkPdfProdAll::MES;
    frame = mes->frame();
  }
  else if (var == Deltae) {
    varFlag = BdkPdfProdAll::DELTAE;
    frame = Deltae->frame();
  }
  else if (var == d0mass) {
    varFlag = BdkPdfProdAll::MD;
    frame = d0mass->frame();
  }

 RooExtendPdf * extPdfVar[BdkEvtTypes::NTYPES];
 RooAbsReal * normVar[BdkEvtTypes::NTYPES];
 RooArgSet extPdfList;
 RooArgSet evtNums;

 RooAddPdf * addPdf=0; 
 int totEvts = 0;

 BdkPdfProdAll * ProdType[BdkEvtTypes::NTYPES];

 for (int t = 0; t < BdkEvtTypes::NTYPES; ++t) {
    normVar[t] = pdfOnResDK.numEvt(t);
    evtNums.add(*(pdfOnResDK.numEvt(t)));
    ProdType[t] = (BdkPdfProdAll *) pdfOnResDK.pdf(t);
    extPdfVar[t] = new RooExtendPdf(TString(ProdType[t]->getVarPdf(varFlag)->GetName()) + ".extPdf",
                                    TString(ProdType[t]->getVarPdf(varFlag)->GetTitle()) + " extPdf",
                                    *(ProdType[t]->getVarPdf(varFlag)->getPdf()),
				    *normVar[t] );
    extPdfList.add(*extPdfVar[t]);
   
    totEvts += (int) normVar[t]->getVal();
   
  }

  addPdf = new RooAddPdf( TString(pdfOnResDK.GetName()) + TString(var->GetName()) + ".sumPdf",
                          TString(pdfOnResDK.GetTitle())+ TString(var->GetTitle()) + " sumPdf",
                          extPdfList, evtNums);

  cout<<" total evts " << totEvts << endl; 
  if(0 == dat) {
	dat = (RooAbsData*) addPdf->generate(RooArgSet(*var), totEvts);
        cout<<" problematic here "<<endl;  
  }

  
  fitOption = "mer";
  RooFitResult * result = 0; 
  if(doFit) {
        result =  addPdf->fitTo(*dat,fitOption);
  }
  if(0 != result ) {	
	result->Print("V");
  }	

  dat->plotOn(frame);
  addPdf->plotOn(frame);
  frame->getAttLine()->SetLineColor(kRed);
  cout<<" chisqare of "<< var->GetName() << " plot: "<<frame->chiSquare()<<endl;
  RooArgSet CompSet(*(extPdfVar[9]),*(extPdfVar[8]));
  addPdf->plotCompOn(frame, CompSet);
  frame->getAttLine()->SetLineColor(kBlue);
  CompSet.add(*(extPdfVar[7]));
  CompSet.add(*(extPdfVar[6]));
  CompSet.add(*(extPdfVar[5]));
  CompSet.add(*(extPdfVar[4]));
  CompSet.add(*(extPdfVar[3]));
  CompSet.add(*(extPdfVar[2]));
  addPdf->plotCompOn(frame, CompSet);
  frame->getAttLine()->SetLineColor(kGreen); 
  CompSet.add(*(extPdfVar[0]));
  addPdf->plotCompOn(frame, CompSet);
  frame->getAttLine()->SetLineColor(kCyan); 
  TCanvas * Cann = new TCanvas(TString(var->GetName()), TString(var->GetTitle()), 600, 600);
  frame->Draw();
  
} 

