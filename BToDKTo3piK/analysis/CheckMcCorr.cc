#include <fstream>

void CheckMcCorr(const char* resultFile, int todo=0,  
	       const char * parInputFile = "analysis/defaultParInputFile.cc",
	       const char * floatFile = "analysis/defaultFloatFile.cc") {
  if( todo > 6 ) {
   cout<<" please input a number not bigger than 6 "<< endl;
   return;
  }  
  int ToDo = 6 - todo; 
  // Execute commands before generating events:
  cout << "reading parInputFile " << parInputFile << endl;
  gROOT->Macro(parInputFile);

  // Fix everything, then execute the file that determines what's floating:
  pdfOnResDK.fixAll();
  cout << "executing floatFile " << floatFile << endl;
  gROOT->Macro(floatFile);

  // Store initial (true) values:
  RooArgSet paramsFree = pdfOnResDK.parametersFree();
  int initValues[300];

  TIterator *iter= paramsFree.createIterator();
  RooRealVar *arg ;
  int i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    initValues[i] = arg->getVal();
  }
  
// Read DK signal
// weight 0.00655; 
  cout<<"===============  signal sample  ============="<<endl; 
  readCut = cutGoodD;
  RooDataSet * data0 = read(chainDK);
  (RooDataSet *)data = randRep(data0, "randAdd<=0.00655");
//  (RooDataSet *)data = randRep(data0, "randAdd<=1.0");
  int sigdk = data->numEntries(); 
  RooTable * table = data->table(*Hdtrkchge);
  table->Print("V");
  double plusEnt = 0; 
  double minsEnt = 0; 
  for(i=0; i<sigdk; ++i) {
   // Get the original event:
    RooArgSet * kchargeset = data->get(i);
    RooCategory * kk = (RooCategory *) kchargeset->find(Hdtrkchge->GetName()); 
    if(1 == kk->getIndex()) {
	plusEnt +=1;
    }
    if(-1 == kk->getIndex()) {
	minsEnt +=1;
    }	
  }

  delete data0;
  readCut = cutBadD;
  RooDataSet * DKBad0 = read(chainDK);
  RooDataSet * dataBadDK = randRep(DKBad0, "randAdd<=0.00655");
//  RooDataSet * dataBadDK = randRep(DKBad0, "randAdd<=1.0");
  int numBadDk = dataBadDK->numEntries();
  delete DKBad0;
  ((RooDataSet *) data)->append(*dataBadDK);

// Dpi sample: Good one,  DE-mes correlation, BAD one, mD-DE, N2-DE correlation. 
// Bch weight  0.228;
  cout<<"=================   D0 pi  ============="<<endl;
  int numDpi=0;
  int numBadDpi=0;
  RooDataSet * dataDpi = 0;
  if( ToDo < 6 ) {
      readCut = cutGoodD;
      RooDataSet * dataGoodDpi0 = read(chainDpi);
//      RooDataSet * repGoodDpi = replace(dataGoodDpi0, Deltae,pdfOnResDK.DpiGoodD0()); 
      dataDpi = randRep(dataGoodDpi0,"randAdd>=0.228&&randAdd<=.456");
      numDpi = dataDpi->numEntries();
      delete dataGoodDpi0;
 
      readCut = cutBadD;
      RooDataSet * dataBadDpi0 = read(chainDpi);
      RooDataSet * dataBadDpi = randRep(dataBadDpi0,"randAdd>=0.228&&randAdd<=.456");
      numBadDpi = dataBadDpi->numEntries();
      delete dataBadDpi0;
      dataDpi->append(*dataBadDpi);
      ((RooDataSet *) data)->append(*dataDpi);
  }

//Dpipi sample
//Bch weight 0.228
  int numDpipi = 0;
  RooDataSet * dataDpipi = 0;
  if( ToDo < 5 ) {
     cout<<"=================  D0 pi pi  ============="<<endl;
     readCut ="";
     RooDataSet * dataBchDpipi0 = read(chainBchDpipi);
     dataDpipi = randRep(dataBchDpipi0,"randAdd<=0.228");
     delete dataBchDpipi0;

//B0 weight 0.227   
     readCut = "";
     RooDataSet * dataB0Dpipi0 = read(chainB0Dpipi);
     RooDataSet * dataB0Dpipi = randRep(dataB0Dpipi0,"randAdd<=0.227");
     dataDpipi->append(*dataB0Dpipi);
     numDpipi = dataDpipi->numEntries();
     delete dataB0Dpipi0;
     ((RooDataSet *) data)->append(*dataDpipi);
  }
   

//charmless sample
//bch weight 0.228

  int numChrmless = 0;
  RooDataSet * dataChrmless = 0; 
  if( ToDo < 4 ) { 
     cout<<"=================  Charmless  ============="<<endl;
     readCut ="";
     RooDataSet * dataChrmless0 = read(chainCharmless);
     dataChrmless = randRep(dataChrmless0,"randAdd<=0.228");
     numChrmless = dataChrmless->numEntries();
     delete dataChrmless0;
     ((RooDataSet *) data)->append(*dataChrmless);
  }

//BBcomb sample
//BBcomb bad  mes-N2 correlation
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  RooDataSet * dataBBcomb = 0; 
  if( ToDo < 3 && 1 != ToDo ) {
     cout<<"=================  BB combinatoric   ============="<<endl;
     //bch weight 0.228
     readCut = cutBadD;  
     RooDataSet * dataBchBadComb0 = read(chainBchComb);
//     dataBBcomb = randRep(dataBchBadComb0,"0.2<=randAdd&&randAdd<=0.428");
//     delete dataBchBadComb0;
     //b0 weight 0.227
     readCut = cutBadD;  
     RooDataSet * dataB0BadComb0 = read(chainB0Comb);
//     RooDataSet * dataBadB0Comb = randRep(dataB0BadComb0, "0.2<=randAdd&&randAdd<=0.427");
//     dataBBcomb->append(*dataBadB0Comb); 
     dataBchBadComb0->append(*dataB0BadComb0);
     RooDataSet * repBB = replace(dataBchBadComb0, bknnout, pdfOnResDK.BBBadD0());
     dataBBcomb = randRep(dataBchBadComb0,"0.2<=randAdd&&randAdd<=0.428"); 
     numBBcomBad = dataBBcomb->numEntries();
     delete dataB0BadComb0;

     //BBcomb good
     readCut = cutGoodD;
     RooDataSet * dataBchGoodComb0 = read(chainBchComb);
     RooDataSet * dataBchGoodComb = 0;
     RooDataSet * repBchGood = 0;
     if( 0!=dataBchGoodComb0->numEntries() ){
        dataBchGoodComb = randRep(dataBchGoodComb0,"0.2<=randAdd&&randAdd<=0.428");
        numBBcomGood += dataBchGoodComb->numEntries();
        dataBBcomb->append(*dataBchGoodComb);
     }
     delete dataBchGoodComb0;

     //b0 weight 0.227
     readCut = cutGoodD;
     RooDataSet * dataB0GoodComb0 = read(chainB0Comb);
     RooDataSet * dataB0GoodComb = 0;
     if( 0!=dataB0GoodComb0->numEntries() ){
         dataB0GoodComb=randRep(dataB0GoodComb0,"0.2<=randAdd&&randAdd<=0.427");
         numBBCombGood+=dataB0GoodComb->numEntries();
         dataBBcomb->append(*dataB0GoodComb);
     }
     delete dataB0GoodComb0;
     ((RooDataSet *) data)->append(*dataBBcomb);
  }

//Continuum Sample
//continuum bad
// cc  weight = 0.629;
  int numQQbad = 0;
  int numQQgood = 0;
  RooDataSet * dataCont = 0; 
  if( ToDo < 2 ) {
      cout<<"=================   comtinuum   ============="<<endl;
      readCut = cutBadD;
      RooDataSet * dataBadCc0 = read(chainCc);
      dataCont = randRep(dataBadCc0, "randAdd<=0.629");
// uds weight = 0.645;
      RooDataSet * dataBadUds0 = read(chainUds);
      RooDataSet * dataBadUds = randRep(dataBadUds0, "randAdd<=0.645");
      dataCont->append(*dataBadUds);
      numQQbad = dataCont->numEntries();

//continuum good, mD-DE correlation
      readCut = cutGoodD;
      RooDataSet * dataGoodCc0 = read(chainCc);
      RooDataSet * dataGoodCc = randRep(dataGoodCc0, "randAdd<=0.629");
      numQQgood = dataGoodCc->numEntries();
      dataCont->append(*dataGoodCc);
      delete dataGoodCc0;
      ((RooDataSet *) data)->append(*dataCont);
  }

  cout<<"================= end reading and loading samples =============="<<endl;

  cout<<" D0K sample: Good "<< sigdk << " Bad "<<numBadDk<<endl;
  cout<<" Dpi sample: Good "<< numDpi <<" Bad "<<numBadDpi<<endl;
  cout<<" Dpp sample: " << numDpipi <<endl;
  cout<<" charmless : " << numChrmless <<endl;
  cout<<" Comb sample: Good "<< numBBcomGood <<" Bad "<< numBBcomBad <<endl;
  cout<<" Cont sample: Good "<< numQQgood <<" Bad "<< numQQbad <<endl;
  
  if (0 == data || data->numEntries() == 0) {
    cout << "emptry data set" << endl;
    return;
  }

  cout<<"==============set up the initial fit value =========="<<endl;

  RooRealVar * bbb = pdfOnResDK.find("pdfOnResDK.BBBadD0NumEvts");
  RooRealVar * bbr = pdfOnResDK.find("pdfOnResDK.BBGoodD0Frac");
  RooRealVar * dpp = pdfOnResDK.find("pdfOnResDK.DPiPiNumEvts");
  RooRealVar * dpg = pdfOnResDK.find("pdfOnResDK.DpiGoodD0NumEvts");
  RooRealVar * dpr = pdfOnResDK.find("pdfOnResDK.DpiBadD0Frac");
  RooRealVar * chm = pdfOnResDK.find("pdfOnResDK.charmlessNumEvts");
  RooRealVar * qqb = pdfOnResDK.find("pdfOnResDK.qqBadD0NumEvts");
  RooRealVar * qqr = pdfOnResDK.find("pdfOnResDK.qqGoodD0Frac");
  RooRealVar * sig = pdfOnResDK.find("pdfOnResDK.sigGoodD0NumEvts");
  RooRealVar * sigr = pdfOnResDK.find("pdfOnResDK.sigBadD0Frac");
  RooRealVar * asym = pdfOnResDK.find("pdfOnResDK.asymSigGoodD");
  sig->setVal(sigdk);
  sigr->setVal(double(numBadDk)/double(sigdk));
  asym->setVal( (plusEnt-minsEnt)/(plusEnt+minsEnt) );
  dpg->setVal(numDpi);
  if(ToDo>5) { 
        dpg->setConstant(kTRUE);
	dpr->setVal(0);
        dpr->setConstant(kTRUE);
  }
  dpp->setVal(numDpipi);
//  if(ToDo>4) { dpp->setConstant(kTRUE); }
  dpp->setConstant(kTRUE); 
  chm->setVal(numChrmless);
  if(ToDo>3) { chm->setConstant(kTRUE); }
  bbb->setVal(numBBcomBad);
  if(ToDo>2 && 1 != ToDo ) {
	bbb->setConstant(kTRUE);
  	bbr->setVal(0);
    	bbr->setConstant(kTRUE);
  } else if( 1 != ToDo) {
	bbr->setVal(double(numBBcomGood)/double(numBBcomBad));
  }
   
  qqb->setVal(numQQbad);
  if(ToDo>1) {
	qqb->setConstant(kTRUE);
 	qqr->setVal(0);
 	qqr->setConstant(kTRUE);
  } 
  cout<<"========== end of input the intial values =========="<<endl;

  // Fit:
  cout<<"Fitting the variables"<<endl;
  paramsFree.Print("V");
  RooFitResult * result = fit(pdfOnResDK, *data);

  // Report the results:
  ofstream ofs(resultFile);
  Int_t status = -100;
  if(0 != result ) status = result->status(); 
  ofs << "fitStatus " << status << " " << paramsFree.getSize() ;

  iter= paramsFree.createIterator();
  i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    double finalValue = arg->getVal();
    double error = arg->getError();
    ofs << " " << initValues[i] << " " << finalValue << " " << error << " ";
  }
  ofs << endl;
  
  if (doPlot) {
    // Plot: 
    RooArgSet chargeSet(*Hdtrkchge);
    RooDataSet redData("redData", "reduced Data for plotting", chargeSet);
    Hdtrkchge->setIndex(-1);
    redData.add(chargeSet);
    Hdtrkchge->setIndex(1);
    redData.add(chargeSet);
 
    data->plotOn(mesFrame);
    pdfOnResDK.getPdf()->plotOn(mesFrame,"L", 1.0, RooAbsPdf::Relative,&redData);
    mesFrame->getAttLine()->SetLineColor(kRed);

    RooArgSet args(*(pdfOnResDK.qqBadD0()->getPdf()));
    pdfOnResDK.getPdf()->plotCompOn(mesFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    mesFrame->getAttLine()->SetLineColor(kRed);

    data->plotOn(DeltaeFrame);
    pdfOnResDK.getPdf()->plotOn(DeltaeFrame,"L", 1.0, RooAbsPdf::Relative,&redData);
    DeltaeFrame->getAttLine()->SetLineColor(kRed);

    data->plotOn(d0massFrame);
    pdfOnResDK.getPdf()->plotOn(d0massFrame,"L", 1.0, RooAbsPdf::Relative,&redData);
    d0massFrame->getAttLine()->SetLineColor(kRed);
    
    data->plotOn(nnFrame);
    pdfOnResDK.getPdf()->plotOn(nnFrame,"L", 1.0, RooAbsPdf::Relative,&redData);
    nnFrame->getAttLine()->SetLineColor(kRed);
    
    data->plotOn(bknnFrame);
    pdfOnResDK.getPdf()->plotOn(bknnFrame,"L", 1.0, RooAbsPdf::Relative,&redData);
    bknnFrame->getAttLine()->SetLineColor(kRed);
    

    RooArgSet args(*(pdfOnResDK.qqBadD0()->getPdf()));
    pdfOnResDK.getPdf()->plotCompOn(mesFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(DeltaeFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(d0massFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(nnFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(bknnFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    mesFrame->getAttLine()->SetLineColor(kBlue);
    DeltaeFrame->getAttLine()->SetLineColor(kBlue);
    d0massFrame->getAttLine()->SetLineColor(kBlue);
    nnFrame->getAttLine()->SetLineColor(kBlue);
    bknnFrame->getAttLine()->SetLineColor(kBlue);
    args.add(*(pdfOnResDK.BBBadD0()->getPdf()));

    args.add(*(pdfOnResDK.DPiPi()->getPdf()));
    pdfOnResDK.getPdf()->plotCompOn(mesFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(DeltaeFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(d0massFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(nnFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    pdfOnResDK.getPdf()->plotCompOn(bknnFrame, args, "L", 1.0, RooAbsPdf::Relative, &redData);
    mesFrame->getAttLine()->SetLineColor(kGreen);
    DeltaeFrame->getAttLine()->SetLineColor(kGreen);
    d0massFrame->getAttLine()->SetLineColor(kGreen);
    nnFrame->getAttLine()->SetLineColor(kGreen);
    bknnFrame->getAttLine()->SetLineColor(kGreen);
    
    TCanvas *c1 = new TCanvas("c1","test bch",800,600);
    c1->Divide(3,2);
    
    c1->cd(1);
    mesFrame->Draw();
    c1->cd(2);
    DeltaeFrame->Draw();
    c1->cd(3);
    d0massFrame->Draw();
    c1->cd(4);
    nnFrame->Draw();
    c1->cd(5);
    bknnFrame->Draw();  
    c1->Print("FullMc.eps");
  }

  cout << "exiting fitFullMc" << endl;
}



