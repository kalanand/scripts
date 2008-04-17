// For variable var (=mes, Deltae, or d0mass), goes over all the 
// chains and fits them to get the parameters:

RooPlot * getRanParFrame = 0;

void getRanParAll(RooRealVar * var, 
		  int todo = 0,
		  Bool_t writeToTFile = kFALSE) {
  if("" == fitOption || "mer" == fitOption ) {
      fitOption = "rm";
   }

  TString fileName = "getRanPar-";
  fileName += var->GetName();
  fileName += ".root";
 
  // open TFile. User's responsibility to close it (or it closes when
  // exiting root):
  if (kTRUE == writeToTFile) {
    delete tFile;
    tFile = new TFile(fileName, "RECREATE");
  }

  // Determine which PDF to use:
  int varFlag = -1;
  RooPlot * frame = 0;
  if (var == mes) {
    varFlag = BdkPdfProdAll::MES;
    frame = mesFrame;
  }
  else if (var == Deltae) {
    varFlag = BdkPdfProdAll::DELTAE;
    frame = DeltaeFrame;
  }
  else if (var == d0mass) {
    varFlag = BdkPdfProdAll::MD;
    frame = d0massFrame;
  }
  else {
    cout << "getRanParAll(): cannot work with variable " << var->GetName()
	 << ". Exiting." << endl;
    return;
  }


  Deltae->setFitMax(0.06);
  mes->setFitMin(5.272);
  d0mass->setFitMin(1.83);
  d0mass->setFitMax(1.895);

  if (0 == todo || 1 == todo) {
    // DK good D:
    readCut = cutGoodD;
    data = read(chainRanDK);
    getRanPar(data, pdfOnResDK.sigGoodD0Prod().getVarPdf(varFlag), 
	      "DK-good-D", var);
  }

  if (0 == todo || 2 == todo) {
    // DK bad D:
    readCut = cutBadD;
    data = read(chainRanDK);
    getRanPar(data, pdfOnResDK.sigBadD0Prod().getVarPdf(varFlag), 
	      "DK-bad-D", var);
  }

  if (0 == todo || 3 == todo) {
    // Dpi good D:
    readCut = cutGoodD; 
    data = read(chainRanBchDpi);
    getRanPar(data, pdfOnResDK.DpiGoodD0Prod().getVarPdf(varFlag), 
	      "Dpi-good-D", var);
  }

  if (0 == todo || 4 == todo) {
    // Dpi bad D:
    readCut = cutBadD;
    data = read(chainRanBchDpi);
    getRanPar(data, pdfOnResDK.DpiBadD0Prod().getVarPdf(varFlag), 
	      "Dpi-bad-D", var);
  }

 
  if (0 == todo || 5 == todo) {
    // DKother:
    TCut DkCut1 = "((B1decmode>9&&B1decmode<19)||(109<B1decmode&&B1decmode<149))||((B2decmode<19&&B2decmode>9)||(109<B2decmode&&B2decmode<149))";
    RooAbsData * data0 = 0;
    int numDKother = addChunk(data0, chainRanBchComb,DkCut1+cutBadD,1);
    numDKother += addChunk(data0, chainRanB0Comb,DkCut1+cutBadD, 1);
    cout<< " DK(*) Other "<< numDKother <<" entries" <<endl; 
    getRanPar(data0, pdfOnResDK.charmlessProd().getVarPdf(varFlag), 
	      "charmless", var);
  }

  if (0 == todo || 6 == todo) {
    //dpipi
    TCut DpiOthCut = "((B1decmode>99&&B1decmode<109)||B1decmode>149)||((B2decmode<109&&B2decmode>99)||B2decmode>149)"; 
    RooAbsData * data5 = 0;
    int numDpiOther = addChunk(data5, chainRanBchComb,DpiOthCut+cutBadD, 1);
    numDpiOther += addChunk(data5, chainRanB0Comb,DpiOthCut+cutBadD, 1);
    getRanPar(data5, pdfOnResDK.DPiPiProd().getVarPdf(varFlag),
             "BB-Dpipi",var);
    cout<<" D(*) pi has "<<numDpiOther<<" entries. "<<endl;	
    delete data5;
  }

  if (0 == todo || 7 == todo) {
    // BB bad D:
    RooAbsData * data1 = 0;
    TCut Bbo = "B1decmode<=0&&B2decmode<=0";
    int numBBbad = addChunk(data1, chainRanBchComb, Bbo+cutBadD, 1);
    numBBbad += addChunk(data1, chainRanB0Comb,Bbo+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " BB comb bad has "<<numBBbad<<" entries " <<endl;
    getRanPar(data1, pdfOnResDK.BBBadD0Prod().getVarPdf(varFlag), 
	      "BB-comb-bad-D", var);
    delete data1; 
  }

  if (0 == todo || 8 == todo) {
    // BB good D:
    RooAbsData * data2 = 0;
    int numBBgood = addChunk(data2, chainRanBchComb, cutGoodD, 1);
    numBBgood += addChunk(data2, chainRanB0Comb, cutGoodD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " BB comb good has "<<numBBgood<<" entries " <<endl;
    getRanPar(data2, pdfOnResDK.BBGoodD0Prod().getVarPdf(varFlag), 
	      "BB-comb-good-D", var);
    delete data2;  
  }


  if (0 == todo || 9 == todo) {
    // uds + cc bad
    RooAbsData * data3 = 0;
    int numQQbad = addChunk(data3, chainRanCc, cutBadD, 1);
    numQQbad += addChunk(data3, chainRanUds, cutBadD, (WEIGHT_UDS/WEIGHT_CC)); 
    cout << " continuum bad D has " << numQQbad << " entries " <<endl; 
    getRanPar(data3, pdfOnResDK.qqBadD0Prod().getVarPdf(varFlag), 
	      "Cont-Bad-D", var);
    delete data3;
  }

  if (0 == todo || 10 == todo) {
    // cc bad D:
    RooAbsData * data4 = 0;
    int numQQgood = addChunk(data4, chainRanCc,cutGoodD, 1);
    numQQgood += addChunk(data4, chainRanUds, cutGoodD, (WEIGHT_UDS/WEIGHT_CC));
    cout << " continuum good D has " << numQQgood << " entries " <<endl; 
    getRanPar(data4, pdfOnResDK.qqGoodD0Prod().getVarPdf(varFlag), 
	      "Cont-Good-D", var);
    delete data4;  
  }


  if (0 == todo || 11 == todo) {
    // off resonance:
    readCut = "";
    data = read(chainRanOffData);
    getRanPar(data, pdfOnResDK.qqBadD0Prod().getVarPdf(varFlag), 
	      "off-res", var);
  }


    
  if (kTRUE == writeToTFile) {
    if(0 != tFile) tFile.Close();
  }
}


//-------------------------------------------------------------
// Does common tasks for all 3 variables:
void getRanPar(RooAbsData* dat, 
	       BdkPdfAbsBase * pdf, 
	       const char * title,
	       RooRealVar * var) {

  cout << "============================================================"
       << endl
       << "=== Fitting \"" << title << "\"" << endl
       << "=== " << var->GetName() << " distribution" << endl
       << "============================================================"
       << endl;

  if (0 == dat || dat->numEntries() == 0) {
    cout << "empty data set" << endl;
    return;
  }

  // adaptive setting of # of bins:
  int nbins = 70;
  if (dat->numEntries() < 1500) {
    nbins = 30;
  }

  dat->SetTitle(title);

  // Fit: 
  fitResult = fit(*pdf, *dat);
  if (0 != fitResult) {
    fitResult->Print();
    if (0 != tFile) {
      fitResult->Write();
    }
  }

  // Plot:  

  TCanvas * can = new TCanvas("can", "can", 400, 400);
  RooPlot * getRanParFrame = 0;
  getRanParFrame = var->frame(nbins);
  getRanParFrame->SetName(title);
  getRanParFrame->SetTitle(title);
  dat->plotOn(getRanParFrame);
  pdf->getPdf()->plotOn(getRanParFrame);
  getRanParFrame->Draw();

  // Print out some info:
  cout << "chi^2 of the plot = " << getRanParFrame->chiSquare() << endl;

  // Save into files:
  TString baseFile = var->GetName();
  baseFile += ".";
  baseFile += title;

  TString parFile = baseFile + ".par";
  pdf->parameters().writeToFile(parFile);

  if (0 != fitResult) {
    ofstream parF(parFile, ios_base::app);
    printFitResult(fitResult, parF);
  }

  TString epsFile = baseFile + ".eps";
  can->SaveAs(epsFile);

  if (0 != tFile) {
    getRanParFrame->Write();  
  }
}

