// $Id: getParams.cc,v 1.8 2006/06/19 17:58:16 fwinkl Exp $
// For variable var (=mes, Deltae, or d0mass), goes over all the 
// chains and fits them to get the parameters:
//
// To use this script for the mES sideband do:
// change cutSigReg -> cutSBmES
// setupChains(false)
// setupPlots()
// mes.setRange(5.2,5.3)


RooPlot * getParamsFrame = 0;

void getParamsAll(RooRealVar * var, 
		  int todo = 0,
		  Bool_t writeToTFile = kFALSE) {
  fitOption = "rm";
  optOption = "c2";

  pdfOnResDK.useNnComb();
  pdfOnResDK.useMd();
  pdfOnResDK.useMes();

  readOnResDKPar(false);
  fixxy();

  TString fileName = "getParams-";
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
  /*
  if (var == mes) {
    varFlag = BdkPdfProdAll::MES;
    frame = mesFrame;
  }
  */
  if (var == Deltae) {
    varFlag = BdkPdfProdAll::DELTAE;
    frame = DeltaeFrame;
  }
  else if (var == d0mass) {
    varFlag = BdkPdfProdAll::MD;
    frame = d0massFrame;
  }
  else if (var == qprime) {
    varFlag = BdkPdfProdAll::NNCONT;
    frame = qprimeFrame;
  }
  else if (var == dprime) {
    varFlag = BdkPdfProdAll::NNCOMB;
    frame = dprimeFrame;
  }
  else {
    cout << "getParamsAll(): cannot work with variable " << var->GetName()
	 << ". Exiting." << endl;
    return;
  }

  if (todo & SIG_G_BIT) {
    // DK good D:
    getParams(sigTree, pdfOnResDK.sigGoodD0N().getVarPdf(varFlag), 
	      cutDKGoodD, var, frame);
  }

  if (todo & SIG_B_BIT) {
    // DK bad D:
    getParams(sigTree, pdfOnResDK.sigBadD0N().getVarPdf(varFlag), 
	      cutDKBadD, var, frame);
  }

  if (todo & DPI_G_BIT) {
    // Dpi good D:
    getParams(dpiTree, pdfOnResDK.DpiGoodD0N().getVarPdf(varFlag), 
	      cutDPiGoodD, var, frame);
  }

  if (todo & DPI_B_BIT) {
    // Dpi bad D:
    getParams(dpiTree, pdfOnResDK.DpiBadD0N().getVarPdf(varFlag), 
	      cutDPiBadD, var, frame);
  }

  if (todo & DKX_BIT) {
    // DKX
    getParams(bbTree, pdfOnResDK.DKXN().getVarPdf(varFlag), 
	      cutDKX, var, frame);
  }

  if (todo & DPIX_BIT) {
    // DPiX
    getParams(bbTree, pdfOnResDK.DPiXN().getVarPdf(varFlag), 
	      cutDPiX, var, frame);
  }

  if (todo & BB_B_BIT) {
    // BB bad D:
    getParams(bbTree, pdfOnResDK.BBBadD0N().getVarPdf(varFlag), 
	      cutBBBadD, var, frame);
  }

  if (todo & BB_G_BIT) {
    // BB good D:
    getParams(bbTree, pdfOnResDK.BBGoodD0N().getVarPdf(varFlag), 
	      cutBBGoodD, var, frame);
  }

  if (todo & QQ_B_BIT) {
    // qq bad D
    getParams(qqTree, pdfOnResDK.qqBadD0N().getVarPdf(varFlag), 
	      cutqqBadD, var, frame);
  }

  if (todo & QQ_G_BIT) {
    // qq good D
    getParams(qqTree, pdfOnResDK.qqGoodD0N().getVarPdf(varFlag), 
	      cutqqGoodD, var, frame);
  }

  /*
  if (0 == todo || 11 == todo) {
    // off resonance:
    readCut = "";
    getParams(chainOffData, pdfOnResDK.qqBadD0Prod().getVarPdf(varFlag), 
	      "off-res", var, frame);
  }
  */

  if (kTRUE == writeToTFile) {
    if(0 != tFile) tFile->Close();
  }

}


//-------------------------------------------------------------
// Does common tasks for all variables
// The signal region cut is added to cut
void getParams(TTree* tree, 
	       BdkPdfAbsBase* pdf,
               TCut cut,
	       RooRealVar* var, RooPlot* frame) {

  const char* title = cut.GetName();
  cout << "============================================================"
       << endl
       << "=== Fitting \"" << title << "\"" << endl
       << "=== " << var->GetName() << " distribution" << endl
       << "============================================================"
       << endl;

  if (tree->GetEntries() == 0) {
    cout << "empty tree" << endl;
    return;
  }

  readCut = cut + cutSigReg;
  //  readCut = cut + cutSBmES;
  // Read data:
  data = read(tree);
  if (0 == data || data->numEntries() == 0) {
    cout << "empty data set" << endl;
    return;
  }

  data->SetTitle(title);

  // Fit: 
  fitResult = fit(*pdf, *data);
  if (0 != fitResult) {
    //    fitResult->Print();
    if (0 != tFile) {
      tFile->cd();
      fitResult->SetName(TString("fr_")+title);
      fitResult->Write();
    }
  }

  // Plot:  

  TGaxis::SetMaxDigits(3);
  gStyle->SetTitleYOffset(1.4);
  TCanvas * can = new TCanvas("can", "can", 0, 0, 500, 400);
  can->SetLeftMargin(0.12);

  getParamsFrame = (RooPlot*)frame->Clone(title);
  getParamsFrame->SetTitle(title);
  data->plotOn(getParamsFrame);
  pdf->getPdf()->plotOn(getParamsFrame);
  getParamsFrame->GetYaxis()->SetTitleOffset(1.2);
  getParamsFrame->Draw();

  // Print out some info:
  cout << "chi^2 of the plot = " << getParamsFrame->chiSquare() << endl;

  // Save into files:
  TString baseFile = var->GetName();
  baseFile += ".";
  baseFile += title;

  TString parFile = baseFile + ".par";
  pdf->parameters().writeToFile(parFile);
  cout << "Writing fit result to "<<parFile<<endl;

  // Append covariance matrix to par file
  ofstream of;
  of.open(parFile, ios_base::out | ios_base::app);
  printCovMatrix(fitResult,of);
  of.close();

  TString epsFile = baseFile + ".eps";
  can->SaveAs(epsFile);

  TString rootFile = baseFile + ".root";
  can->SaveAs(rootFile);

  if (0 != tFile) {
    getParamsFrame->Write();  
  }
}
