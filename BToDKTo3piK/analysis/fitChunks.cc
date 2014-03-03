// $Id: fitChunks.cc,v 1.3 2006/04/13 02:21:26 fwinkl Exp $
// Create and fit weighted MC sample

void fitChunks(const char* resultFile = 0, 
               int todo = ALL_TYPES, int bin = 0,
	       Replacer replacer = Replacer(),
	       const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
	       const char * parInputFile = "analysis/defaultParInputFile.cc",
	       const char * floatFile = "analysis/defaultFloatFile.cc") {

  RooArgSet* paramsOnResDK = new RooArgSet(pdfOnResDK.parameters());

  //reload the params
  readOnResDKPar();

  // Execute commands before generating events:
  if (parInputFile) {
    cout << "reading parInputFile " << parInputFile << endl;
    gROOT->Macro(parInputFile);
  }

  // Fix everything and float some:
  if (floatFile) {
    pdfOnResDK.fixAll();
    cout << "executing floatFile " << floatFile << endl;
    gROOT->Macro(floatFile);
  }

  delete data; // clean up before starting
  data = 0;

  // Apply the signal cut to all trees (to save time later)
  TCut cut = cutSigReg;  
  cout << endl <<" ============== Precutting the trees ==========="<<endl;  
  cout << "Cut: "<<cut.GetTitle()<<endl;

  gROOT->cd();
  TTree* _sigTree = sigTree->CopyTree(cut);
  TTree* _bpTree = bpTree->CopyTree(cut);
  TTree* _b0Tree = b0Tree->CopyTree(cut);
  TTree* _udsTree = udsTree->CopyTree(cut);
  TTree* _ccTree = ccTree->CopyTree(cut);
  TTree* _dpiTree = dpiTree->CopyTree(cut);


  // Dpi:
  int numDpiGood = 0;
  int numDpiBad = 0;
  if(todo & DPI_G_BIT) {
    cout<< endl<<"=================   DPi good ============="<<endl;
    numDpiGood += addChunk(data, _dpiTree, cutDPiGoodD, DPI_WEIGHT, bin);
  }
  if(todo & DPI_B_BIT) {
    cout<< endl<<"=================   DPi bad ============="<<endl;
    numDpiBad  += addChunk(data, _dpiTree, cutDPiBadD,  DPI_WEIGHT, bin);
  }

  // DPiX:
  int numDPiX = 0;
  if( todo & DPIX_BIT) {
    cout<<"=================  DPiX  ============="<<endl;
    numDPiX += addChunk(data, _bpTree, cutDPiX, BP_WEIGHT, bin);
    numDPiX += addChunk(data, _b0Tree, cutDPiX, B0_WEIGHT, bin);    
  }
   
  // DKX:
  int numDKX = 0;
  if( todo & DKX_BIT) {
    cout<<"=================  DKX  ============="<<endl;
    numDKX += addChunk(data, _bpTree, cutDKX, BP_WEIGHT, bin);
    numDKX += addChunk(data, _b0Tree, cutDKX, B0_WEIGHT, bin);    
  }

  // BBcomb:
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  if( todo & BB_G_BIT) {
    cout<< endl<<"=================  BB good  ============="<<endl;
    numBBcomGood += addChunk(data, _bpTree, cutBBGoodD, BP_WEIGHT, bin);
    numBBcomGood += addChunk(data, _b0Tree, cutBBGoodD, B0_WEIGHT, bin);
  }
  if( todo & BB_B_BIT) {
    cout<< endl<<"=================  BB bad   ============="<<endl;
    RooDataSet * bbData = 0;
    numBBcomBad += addChunk(bbData, _bpTree, cutBBBadD, BP_WEIGHT, bin);
    numBBcomBad += addChunk(bbData, _b0Tree, cutBBBadD, B0_WEIGHT, bin);

    if (replacer._typeBit & BB_B_BIT) {
      //      bbData = replace(bbData, replacer._var, pdfOnResDK.BBBadD0());  
    }

    if (0 != bbData) {
      if (0 != data) {
	((RooDataSet*)data)->append(*bbData);
      }
      else {
	data = bbData;
      }
    }
  }
  
  // Continuum:
  int numQQBad = 0;
  int numQQGood = 0;
  if( todo & QQ_G_BIT) {
    cout<< endl<<"=================   continuum good  ============="<<endl;
    numQQGood += addChunk(data, _udsTree, cutqqGoodD, UDS_WEIGHT, bin);
    numQQGood += addChunk(data, _ccTree, cutqqGoodD, CC_WEIGHT, bin);
  }
  
  // put continuum bad in a different data set, so we can measure its asymmetry:
  RooDataSet * qqData = 0;
  if( todo & QQ_B_BIT) {
    cout<< endl<<"=================   continuum bad  ============="<<endl;    
    numQQBad += addChunk(qqData, _udsTree, cutqqBadD, UDS_WEIGHT, bin);
    numQQBad += addChunk(qqData, _ccTree, cutqqBadD, CC_WEIGHT, bin);
  }
  
  double qqBAsym = 0;

  if (0 != qqData) {
    qqBAsym = trueAsym(qqData);
    // then add qqData to data:
    if (0 != data && 0 != qqData) {
      ((RooDataSet*)data)->append(*qqData);
    }
    else {
      data = qqData;
    }
  }

  Double_t gAsym = 0;
  if (data) {
    RooTable * table = data->table(*Hdtrkchge);
    table->Print("V");
    gAsym = trueAsym(data);
  }
  cout << "Global asymmetry = "<<gAsym<<endl;

  // Read DK signal good & bad:
  int numSigGood = 0;
  Double_t sigAsym = 0;
  if (todo & SIG_G_BIT) {
    cout<< endl << "===============  signal good sample  ============="<<endl; 
    RooDataSet* sigData = 0;
    /*
    numSigGood = addChunk(sigData, _sigTree, cutDKGoodD, SIG_WEIGHT, bin);    
    if (replacer._typeBit & SIG_G_BIT) {
      //      data = replace(sigData, replacer._var, pdfOnResDK.sigGoodD0());  
    }
    */

    // Generate signal events from flat signal MC
    TTree* _sigFlatTree = sigFlatTree->CopyTree(cut+cutDKGoodD);
    RooDataSet* sigFlatBpData = read(_sigFlatTree,0,"Hdtrkchge>0",allVars);
    RooDataSet* sigFlatBmData = read(_sigFlatTree,0,"Hdtrkchge<0",allVars);

    RooDataSet* sigBpData = generateFromFlat(*sigFlatBpData,
                              (BdkDKDalitz*)dalitzHolderP.sigGoodD0()->getPdf());
    RooDataSet* sigBmData = generateFromFlat(*sigFlatBmData,
                              (BdkDKDalitz*)dalitzHolderN.sigGoodD0()->getPdf());
    delete _sigFlatTree;

    Double_t effBp = (Double_t)sigBpData->numEntries()/sigFlatBpData->numEntries();
    Double_t effBm = (Double_t)sigBmData->numEntries()/sigFlatBmData->numEntries();

    cout << "Reweighting efficencies: B+ = "<<effBp<<endl
         << "                         B- = "<<effBm<<endl;

    numSigGood = addChunk(sigData, &sigBpData->tree(), "", SIG_WEIGHT/effBp, bin);
    numSigGood += addChunk(sigData, &sigBmData->tree(), "", SIG_WEIGHT/effBm, bin);

    if(0 != sigData) {  
      RooTable * table = sigData->table(*Hdtrkchge);
      table->Print("V");
      sigAsym = trueAsym(sigData)-gAsym;
      cout << "Signal Asymmetry = "<<sigAsym<<endl;

      if (0 != data) ((RooDataSet*)data)->append(*sigData);
      else data = sigData;        
    }
  }

  int numSigBad = 0; 
  if (todo & SIG_B_BIT) {
    cout<< endl << "===============  signal bad sample  ============="<<endl; 
    numSigBad = addChunk(data, _sigTree, cutDKBadD, SIG_WEIGHT, bin);
  }




  cout<<"================ end reading and loading samples ============="<<endl;

  cout<<" D0K sample: Good "<< numSigGood << " Bad "<< numSigBad << endl;
  cout<<" Dpi sample: Good "<< numDpiGood <<" Bad "<< numDpiBad << endl;
  cout<<" DPiX sample: " << numDPiX <<endl;
  cout<<" DKX sample: " << numDKX <<endl;
  cout<<" BB sample: Good "<< numBBcomGood <<" Bad "<< numBBcomBad <<endl;
  cout<<" qq sample: Good "<< numQQGood <<" Bad "<< numQQBad <<endl;
   
  if (0 == data || data->numEntries() == 0) {
    cout << "empty data set" << endl;
    return;
  }

  
  cout<<"============== setup the initial fit values =========="<<endl;

  Int_t totBBNumEvts = numDKX + numDPiX + numBBcomBad;

  setVar(paramsOnResDK, "pdfOnResDK.qqGoodD0Frac",
	 numQQBad ? (double)numQQGood/numQQBad : 0);

  setVar(paramsOnResDK, "pdfOnResDK.BBGoodD0Frac", 
	 numBBcomBad ? (double)numBBcomGood/numBBcomBad : 0);

  setVar(paramsOnResDK, "pdfOnResDK.DpiBadD0Frac", 
         numDpiGood ? (double)numDpiBad/numDpiGood : 0);

  setVar(paramsOnResDK, "pdfOnResDK.DKXFrac", 
         numDPiX ? (double)numDKX/numDPiX : 0);

  setVar(paramsOnResDK, "pdfOnResDK.DPiXFrac", 
         totBBNumEvts ? (double)numDPiX/totBBNumEvts : 0);

  setVar(paramsOnResDK, "pdfOnResDK.sigBadD0Frac", 
         numSigGood ? (double)numSigBad/numSigGood : 0);


  setVarFix0(paramsOnResDK, "pdfOnResDK.DpiGoodD0NumEvts", numDpiGood);
 
  setVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD",sigAsym);
  if (0 == numSigGood){
    fixVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD");
  }

  setVarFix0(paramsOnResDK, "pdfOnResDK.qqBadD0NumEvts", numQQBad);
  setVarFix0(paramsOnResDK, "pdfOnResDK.sigGoodD0NumEvts", numSigGood); 
  setVarFix0(paramsOnResDK, "pdfOnResDK.totBBNumEvts", totBBNumEvts);



  setVar(paramsOnResDK, "pdfOnResDK.asymQqBadD", qqBAsym - gAsym);
  if (0 == numQQBad){
    fixVar(paramsOnResDK, "pdfOnResDK.asymQqBadD");
  }

  setVar(paramsOnResDK, "pdfOnResDK.globalAsym", gAsym);	 
  
 
  cout<<"========== end of setting up the intial values =========="<<endl;

  // fix more vars:
  if (moreFixFile) {
    cout << "executing moreFixFile " << moreFixFile << endl;
    gROOT->Macro(moreFixFile);
  }

  // Store initial (true) values:
  RooArgSet* initValues = (RooArgSet*)pdfOnResDK.parametersFree().snapshot(false);

  // Fit:
  cout<<"Fitting the variables"<<endl;
  pdfOnResDK.parametersFree().Print("v");
  RooFitResult* result = 0;
  result = fit(pdfOnResDK, *data);

  // Report the results:
  ostringstream os;

  Int_t status = -100;
  Double_t minNll = 100;
  if(0 != result ) {
    status = result->status();
    minNll = result->minNll();
  }
  os << "fitStatus " << status << " " << minNll
     << " " << initValues->getSize()<<endl;
                                                                                
  TIterator* iter = pdfOnResDK.parametersFree().createIterator();
  RooRealVar* arg;
  while(arg = (RooRealVar*)iter->Next()) {
    double initValue = ((RooRealVar*)initValues->find(arg->GetName()))->getVal();
    double finalValue = arg->getVal();
    double error = arg->getError();
    os << arg->GetName() << " " << initValue << " "
       << finalValue << " " << error <<endl;
  }

  cout << os.str();

  if (resultFile) {
    ofstream ofs(resultFile);
    ofs << os.str();
    ofs.close();
  }
  
//  stdPlot(data);

  delete paramsOnResDK;
  delete _sigTree;
  delete _bpTree;
  delete _b0Tree;
  delete _ccTree;
  delete _udsTree;
  delete _dpiTree;

  cout << "exiting fitChunks, todo = " << printBits(todo) << endl;
}
