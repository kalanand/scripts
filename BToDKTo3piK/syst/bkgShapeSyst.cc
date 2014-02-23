// $Id: bkgShapeSyst.cc,v 1.6 2006/07/11 21:44:03 fwinkl Exp $
//
// Systematics for data-MC agreement of 1D and Dalitz fit variables
// Fit data in mES sideband with sum of all event types
// Write out .par file that can be used with systAnyParFile()

RooAbsCollection* origSBParams;  // nominal sideband region parameters
RooAbsCollection* origParams;    // nominal signal region paramters

void bkgShapeSyst(Bool_t doPlot = kFALSE) {

  // use random selection ntuples with sidebands
  setupChains(false);

  // Initialize and save original parameters to calculate difference later
  useBothFitVars();
  readOnResDKPar();
  origParams = pdfOnResDK.parameters().snapshot(false);

  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/SB-all.par");
  origSBParams = pdfOnResDK.parameters().snapshot(false);

  // Don't read mes because the range of the RooRealVar will kill the sideband
  RooArgSet vars(*Deltae,*bknnout,*nnout,*mass12,*mass13,*Hdtrkchge);

  RooDataSet* d = read(dataTree,0,cutSBmES,&vars);
  data = makeDalitzDataSet(*d);    // flip B+ Dalitz variables
  delete d;
  
  doEvtType(BdkEvtTypes::DPiX, "DPiX", doPlot);
  doEvtType(BdkEvtTypes::DKX, "DKX", doPlot);
  doEvtType(BdkEvtTypes::BB_BAD_D, "BBbad", doPlot);
  doEvtType(BdkEvtTypes::QQ_BAD_D, "QQbad", doPlot);
}


// do everything for one event type
void doEvtType(BdkEvtTypes::Type floatEvtType, TString name, Bool_t doPlot = kFALSE)
{
  /*
  fitVarAll(*data, floatEvtType, doPlot, name);
  RooArgSet set1(pdfOnResDK.parametersFree());
  fixAll(set1); 
  set1.writeToFile(TString("bkgShapeSyst-1D-")+name+".par");
  */

  fitVar(*data, BdkPdfProdAll::DELTAE, floatEvtType, doPlot, name);
  fitVar(*data, BdkPdfProdAll::NNCONT, floatEvtType, doPlot, name);
  fitVar(*data, BdkPdfProdAll::NNCOMB, floatEvtType, doPlot, name);
  //fitVar(*data, BdkPdfProdAll::DALITZ, floatEvtType, doPlot, name);
}


// fit var in data, only float floatEvtType
void fitVar(RooDataSet& data,
            BdkPdfProdAll::Var var, 
            BdkEvtTypes::Type floatEvtType,
            Bool_t doPlot = kFALSE,
            TString name = "")

{
  // Disable all vars
  pdfOnResDK.useDalitz(false);
  pdfOnResDK.useNnCont(false);
  pdfOnResDK.useNnComb(false);
  pdfOnResDK.useDE(false);
  pdfOnResDK.useMes(false);
  pdfOnResDK.useMd(false);
  // Enable selected var
  pdfOnResDK.useVar(var);
  pdfOnResDK.getPdf();

  readOnResDKPar(false);
  // read fit fractions in SB
  //  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/numEvts_mesSB.par");
  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/SB-all.par");


  // Fix all parameters except the ones from floatEvtType
  for (int i=0; i<BdkEvtTypes::NTYPES; i++) {
    if (i != floatEvtType) pdfOnResDK.prodN(i)->fixAll();
    // The +/- PDFs are not shared for the Dalitz. Fix the positive PDF.
    if (var==BdkPdfProdAll::DALITZ) pdfOnResDK.prodP(i)->fixAll();
  }

  BdkPdfAbsBase* bdkPdf = pdfOnResDK.prodN(floatEvtType)->getVarPdf(var);
  
  // for NN variables fix the bifurcated Gaussian if fraction is small
  if (var==BdkPdfProdAll::NNCONT || var==BdkPdfProdAll::NNCOMB) {
    if (((BdkPdfBifurGaussGauss*)bdkPdf)->fracGauss()->getVal()>0.8)
      ((BdkPdfBifurGaussGauss*)bdkPdf)->bifur().fixAll();
  }

  

  // print floating parameters
  bdkPdf->parametersFree().Print("v");

  if (bdkPdf->parametersFree().getSize()==0) {
    cout << "No free parameters in "<<bdkPdf->GetName()<<". Skipping."<<endl;
    return;
  }

  // get the associated RooRealVar
  RooRealVar* r = (RooRealVar*)(RooArgList(bdkPdf->dependents())).at(0);
  
  // only use the negative part of the PDF
  RooAbsPdf* pdf = pdfOnResDK.addPdfNeg();

  RooPlot* p;
  if (doPlot) {
    // plot original PDF
    p = r->frame(Range("plot"));
    data.plotOn(p);
    pdf->plotOn(p,LineStyle(kDashed));
  }
   
  // fit
  fitOption = "mr";
  optOption = "c3";
  fit(*pdf, data);

  TString filename = "bkgShapeSyst-"+TString(r->GetName())+"-"+name;

  // Calculate difference to original SB parameters
  RooArgSet set(pdfOnResDK.parametersFree());
  RooArgSet* diff = (RooArgSet*)origSBParams->selectCommon(set)->snapshot(false);
  sub(*diff,set);
  cout << "Difference to original sideband parameters:"<<endl;
  diff->Print("v");

  // Apply difference to signal region parameters
  RooArgSet* newParams = (RooArgSet*)origParams->selectCommon(set)->snapshot(false);
  sub(*newParams, *diff);

  fixAll(*newParams);
  newParams->writeToFile(filename+".par");
  delete diff;
  delete newParams;

  if (doPlot) {
    // plot new PDF
    pdf->plotOn(p);

    TCanvas *can = new TCanvas("can","",400,400);
    p->SetTitle("");
    p->Draw();
    can->SaveAs(filename+".eps");
  }

}



// only float floatEvtType
void fitVarAll(RooDataSet& data,
               BdkEvtTypes::Type floatEvtType,
               Bool_t doPlot = kFALSE,
               TString name = "")

{
  // Disable all vars
  pdfOnResDK.useDalitz(false);
  pdfOnResDK.useNnCont();
  pdfOnResDK.useNnComb();
  pdfOnResDK.useDE();
  pdfOnResDK.getPdf();

  readOnResDKPar(false);
  // read fit fractions in SB
  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/numEvts_mesSB.par");

  // Fix all parameters except the ones from floatEvtType
  for (int i=0; i<BdkEvtTypes::NTYPES; i++) {
    if (i != floatEvtType) pdfOnResDK.prodN(i)->fixAll();
  }

  BdkPdfAbsBase* bdkPdf = pdfOnResDK.prodN(floatEvtType);

  if (((BdkPdfBifurGaussGauss*)pdfOnResDK.prodN(floatEvtType)->getNnContPdf())->fracGauss()->getVal()>0.7)
    ((BdkPdfBifurGaussGauss*)pdfOnResDK.prodN(floatEvtType)->getNnContPdf())->bifur().fixAll();

  if (((BdkPdfBifurGaussGauss*)pdfOnResDK.prodN(floatEvtType)->getNnCombPdf())->fracGauss()->getVal()>0.7)
    ((BdkPdfBifurGaussGauss*)pdfOnResDK.prodN(floatEvtType)->getNnCombPdf())->bifur().fixAll();


  // print floating parameters
  bdkPdf->parametersFree().Print("v");

  if (bdkPdf->parametersFree().getSize()==0) {
    cout << "No free parameters in "<<bdkPdf->GetName()<<". Skipping."<<endl;
    return;
  }

  // only use the negative part of the PDF
  RooAbsPdf* pdf = pdfOnResDK.addPdfNeg();

  // fit
  fitOption = "mr";
  optOption = "c4";
  fit(*pdf, data);


}
