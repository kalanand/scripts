void getBinD0Mass(char * filename) 
{
RooRealVar * d0mass = new RooRealVar("d0mass", "d0 mass", 1.8645, 1.805, 1.924); 
TFile * f = new TFile(filename);
//TFile * f = new TFile("d0pi_3pi_data.root");
TH1F * histo = (TH1F*) f->Get("h500");  
TCanvas * can = new TCanvas("newcan", "new Can", 600, 600);
histo->Draw();
RooDataHist bdata("bdata", "d0 mass binned", RooArgList(*d0mass), histo);
BdkPdfCBPolyn * d0histPdf = new BdkPdfCBPolyn("d0histPdf", "d0 mass pdf", *d0mass, 0, 2);
d0histPdf->parameters()->readFromFile("D0pi3pi-D0mass.par");
d0histPdf->parameters()->Print("V");
RooFitResult * fitResult = 0;
fitResult = d0histPdf->getPdf()->fitTo(bdata, "mr");
if(0 != fitResult ) {  fitResult->Print("V"); }
RooPlot * d0massFrame = d0mass->frame();
bdata->plotOn(d0massFrame);
d0histPdf->getPdf()->plotOn(d0massFrame);
cout<<"plot chisquare "<<d0massFrame->chiSquare()<<endl;
RooArgSet args(*(d0histPdf->polyn().getPdf()));
d0histPdf->getPdf()->plotOn(d0massFrame, Components(args), DrawOption("L"), LineColor(kRed));
d0massFrame->Draw();
}
  


