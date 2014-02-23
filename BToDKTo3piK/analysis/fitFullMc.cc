#include <fstream>

void fitFullMc(const char* resultFile, 
	       const char * parInputFile = "analysis/defaultParInputFile.cc",
	       const char * floatFile = "analysis/defaultFloatFile.cc") {
  
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
  int numBadDk = dataBadDK->numEntries();
  delete DKBad0;

// Dpi sample  
// Bch weight  0.228;
  cout<<"=================   D0 pi  ============="<<endl;
  readCut = cutGoodD;
  RooDataSet * dataGoodDpi0 = read(chainDpi);
  RooDataSet * dataGoodDpi = randRep(dataGoodDpi0,"randAdd>=0.228&&randAdd<=.456");
  int numDpi = dataGoodDpi->numEntries();
  delete dataGoodDpi0;
 
  readCut = cutBadD;
  RooDataSet * dataBadDpi0 = read(chainDpi);
  RooDataSet * dataBadDpi = randRep(dataBadDpi0,"randAdd>=0.228&&randAdd<=.456");
  int numBadDpi = dataBadDpi->numEntries();
  delete dataBadDpi0;

//Dpipi sample
//Bch weight 0.228
  cout<<"=================  D0 pi pi  ============="<<endl;
  readCut ="";
  RooDataSet * dataBchDpipi0 = read(chainBchDpipi);
  RooDataSet * dataBchDpipi = randRep(dataBchDpipi0,"randAdd<=0.028");
  delete dataBchDpipi0;

//B0 weight 0.227   
  readCut = "";
  RooDataSet * dataB0Dpipi0 = read(chainB0Dpipi);
  RooDataSet * dataB0Dpipi = randRep(dataB0Dpipi0,"randAdd<=0.227");
  int numDpipi = dataB0Dpipi->numEntries() + dataBchDpipi->numEntries();
  delete dataB0Dpipi0;
  
//charmless sample
//bch weight 0.228
  cout<<"=================  Charmless  ============="<<endl;
  readCut ="";
  RooDataSet * dataChrmless0 = read(chainCharmless);
  RooDataSet * dataChrmless = randRep(dataChrmless0,"randAdd<=0.228");
  int numChrmless = dataChrmless->numEntries();
  delete dataChrmless0;

//BBcomb sample
//BBcomb bad
//bch weight 0.228
  cout<<"=================  BB combinatoric   ============="<<endl;
  readCut = cutBadD;  
  RooDataSet * dataBchBadComb0 = read(chainBchComb);
//  RooDataSet * dataBchBadComb = randRep(dataBchBadComb0,"randAdd<=0.228");
//  delete dataBchBadComb0;
//b0 weight 0.227
  RooDataSet * dataB0BadComb0 = read(chainB0Comb);
//  RooDataSet * dataB0BadComb = randRep(dataB0BadComb0,"randAdd<=0.227");
  dataBchBadComb0->append(*dataB0BadComb0);
  RooDataSet * dataBBcomb0 = replace(dataBchBadComb0, bknnout, pdfOnResDK.BBBadD0());
  RooDataSet * dataBBcomb = randRep(dataBBcomb0, "randAdd<=0.228"); 
//  int numBBcomBad = dataBchBadComb->numEntries() + dataB0BadComb->numEntries();
  int numBBcomBad = dataBBcomb->numEntries();
  delete dataBchBadComb0;
  delete dataB0BadComb0;

//BBcomb good
  readCut = cutGoodD;
  RooDataSet * dataBchGoodComb0 = read(chainBchComb);
  RooDataSet * dataBchGoodComb = 0;
  int numBBcomGood = 0;
  if( 0!=dataBchGoodComb0->numEntries() ){
     dataBchGoodComb = randRep(dataBchGoodComb0,"randAdd<=0.228");
     numBBcomGood += dataBchGoodComb->numEntries();       
  }
  delete dataBchGoodComb0;

//b0 weight 0.227
  RooDataSet * dataB0GoodComb0 = read(chainB0Comb);
  RooDataSet * dataB0GoodComb = 0;
  if( 0!=dataB0GoodComb0->numEntries() ){
    dataB0GoodComb=randRep(dataB0GoodComb0,"randAdd<=0.227");
    numBBCombGood+=dataB0GoodComb->numEntries();
  }
  delete dataB0GoodComb0;

//Continuum Sample
//continuum bad
// cc  weight = 0.629;
  cout<<"=================   comtinuum   ============="<<endl;
  readCut = cutBadD;
  RooDataSet * dataBadCc0 = read(chainCc);
  RooDataSet * dataBadCc = randRep(dataBadCc0, "randAdd<=0.629");
// uds weight = 0.645;
  RooDataSet * dataBadUds0 = read(chainUds);
  RooDataSet * dataBadUds = randRep(dataBadUds0, "randAdd<=0.645");
  int numQQbad = dataBadCc->numEntries()+dataBadUds->numEntries();

//continuum good
  readCut = cutGoodD;
  RooDataSet * dataGoodCc0 = read(chainCc);
  RooDataSet * dataGoodCc = randRep(dataGoodCc0, "randAdd<=0.629");
  int numQQgood = dataGoodCc->numEntries();

  cout<<"================= end reading samples =============="<<endl;

  cout<<"============ start to load data for fitting ========="<<endl;
  ((RooDataSet *) data)->append(*dataBadDK);
  ((RooDataSet *) data)->append(*dataGoodDpi);
  ((RooDataSet *) data)->append(*dataBadDpi);
  ((RooDataSet *) data)->append(*dataB0Dpipi);
  ((RooDataSet *) data)->append(*dataBchDpipi);
  ((RooDataSet *) data)->append(*dataChrmless);
  if(0!= dataBchGoodComb) {
     ((RooDataSet *) data)->append(*dataBchGoodComb);
  } 
  if(0 != dataB0GoodComb) { 
	 ((RooDataSet *) data)->append(*dataB0GoodComb);
  }	
  ((RooDataSet *) data)->append(*dataBBcomb);
  ((RooDataSet *) data)->append(*dataBadCc);
  ((RooDataSet *) data)->append(*dataBadUds);
  ((RooDataSet *) data)->append(*dataGoodCc);
  cout<<"============== end of loading data ================="<<endl;

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
  RooRealVar * dpp = pdfOnResDK.find("pdfOnResDK.DPiPiNumEvts");
  RooRealVar * dpg = pdfOnResDK.find("pdfOnResDK.DpiGoodD0NumEvts");
  RooRealVar * chm = pdfOnResDK.find("pdfOnResDK.charmlessNumEvts");
  RooRealVar * qqb = pdfOnResDK.find("pdfOnResDK.qqBadD0NumEvts");
  RooRealVar * sig = pdfOnResDK.find("pdfOnResDK.sigGoodD0NumEvts");
  RooRealVar * asym = pdfOnResDK.find("pdfOnResDK.asymSigGoodD");
  bbb->setVal(numBBcomBad);
  dpp->setVal(numDpipi);
  dpg->setVal(numDpi);
  chm->setVal(numChrmless);
  qqb->setVal(numQQbad);
  sig->setVal(sigdk);
  asym->setVal( (plusEnt-minsEnt)/(plusEnt+minsEnt) );
 
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



