// plots mes distribution of BBcomb MC and overlays it with the pdf shape from the
// MC and from the data fit:

void compareBBParams() 
{
  data = read(chainBBComb);  
  BdkPdfAbsBase & pdf = pdfOnResDK.BBBadD0Prod();

  data->plotOn(mesFrame);

  // plot MC parameters:
  pdf.getPdf()->plotOn(mesFrame);
  
  // get and plot data parameters:
  pdf.parameters().readFromFile("analysis/defaultParInputFileData.par");
  pdf.getPdf()->plotOn(mesFrame);
  mesFrame->getAttLine()->SetLineColor(kRed);
  mesFrame->getAttLine()->SetLineStyle(2);

  TCanvas * can = new TCanvas("can", "can", 800, 500);
  mesFrame->Draw();
  can->Print("compareBBParams.eps");
}
