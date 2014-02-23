// $Id: rndToy.cc,v 1.1 2007/05/11 12:31:25 fwinkl Exp $
// Full toy MC study with random values for rB, gamma and delta

RooDataSet* proto;
RooAbsPdf* pdf;

// BdkBatchMCStudy doesn't know how to deal with pointers
RooCategory kcharge(*Hdtrkchge,"*Hdtrkchge");

// setup the PDF
void setupPdf() {

  proto = 0;   // don't need it for pdfOnResDK toys

  pdf = pdfOnResDK.getPdf();
  readOnResDKPar();

  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/results/fitData.par");

  // No systematic errors
  pdfOnResDK.setNLLYieldsSystBit(0);

  // Set correct names for BdkBatchMCStudy
  pdf->SetName("pdfOnResDK");
  m12->SetName("*m12");
  m13->SetName("*m13");
  Deltae->SetName("*Deltae");
  qprime->SetName("*qprime");
  dprime->SetName("*dprime");
}
