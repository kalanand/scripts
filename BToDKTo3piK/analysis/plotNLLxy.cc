// $Id: plotNLLxy.cc,v 1.1 2006/05/29 22:56:48 fwinkl Exp $
// Plot NLL for different values of x and y

void plotNLLxy() {

  TCanvas* c = new TCanvas("c","plotNLLxy",1200,600);
  c->Divide(5,2);
 
  BdkPdfDKDalitz& pdf = dalitzHolderN.sigGoodD0Type();

  pdf.x()->setRange(-4,4);
  pdf.y()->setRange(-4,4);

  int i = 1;
  for (double x = -0.2; x<=0.2; x+=0.05, i++) {

    pdf.x()->setVal(x);
    pdf.y()->setVal(x);

    data = pdf.generate(1000);
  
    c->cd(i);
    RooPlot* px = pdf.getPdf()->plotNLLOn(pdf.x()->frame(),data);
    RooPlot* py = pdf.getPdf()->plotNLLOn(pdf.y()->frame(),data);
    
    px->getAttLine()->SetLineColor(kBlue);
    py->getAttLine()->SetLineColor(kRed);

    px->Draw();
    py->Draw("same");
 
    delete data;
    c->Update();
  }
  c->SaveAs("plotNLLxy.eps");
}
