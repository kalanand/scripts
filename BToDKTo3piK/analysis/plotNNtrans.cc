// $Id: plotNNtrans.cc,v 1.2 2006/08/04 22:26:34 fwinkl Exp $
// Plot the NN transformation

void plotNNtrans(TTree* tree, TCut cut) 
{
  TCanvas* c = new TCanvas("plotNNtrans","NN transformation",800,400);  
  c->Divide(2,1);
  
  readCut = cut;
  RooDataSet* d = read(tree);
  
  c->cd(1);
  RooPlot* p1 = nnout->frame(0.1,1.0,20);
  d->plotOn(p1);
  p1->SetTitle("");
  p1->Draw();

  c->cd(2);
  RooPlot* p2 = qprime->frame(-5,5,20);
  d->plotOn(p2);
  p2->SetTitle("");
  p2->Draw();

  c->SaveAs("plotNNtrans.eps");
  delete d;
}
