// $Id: compareDist.cc,v 1.1 2006/05/16 22:13:07 fwinkl Exp $
// Compare two datasets or one dataset with two different cuts

void compareDistAll()
{

  TCanvas * can = new TCanvas("can", "compareDist", 1200, 600);
  can->Divide(4,2,0.004,0.004);
  TGaxis::SetMaxDigits(3);
  can->SetLeftMargin(0.12);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetOptStat(0);

  ostringstream os;

  int i = 1;
  can->cd(i);

  gROOT->cd();
  TTree* cutTree = dataTree->CopyTree(cutSigReg);
  TCut cut1 = cutRun14;
  TCut cut2 = cutRun5;

  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, Deltae, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, mes, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, d0mass, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, nnout, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, bknnout, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, mass12, os);
  can->cd(i++); compare(cutTree, cutTree, cut1, cut2, mass13, os);


  TString file = "compareDist";

  cout << os.str();
  cout << "KS probabilities saved to "<<file<<".txt"<<endl;
  ofstream of(file+".txt");
  of << os.str();
  of.close();

  can->SaveAs(file+".eps");
  can->SaveAs(file+".root");

  delete cutTree;
}



// Compare data in tree1 and tree2 applying cut1 and cut2
RooPlot* compare(TTree* tree1, TTree* tree2,
                 TCut cut1, TCut cut2,
                 RooRealVar* var, ostringstream& os)
{
  const char* title = ""; //cut.GetName();

  RooPlot* plot = 0;

  plot = var->frame();
  //    plot->SetName(title);
  plot->SetTitle(title);
  
  TH1* h1 = var->createHistogram("h1");
  TH1* h2 = var->createHistogram("h2");
  h1->Sumw2();
  h2->Sumw2();
  
  tree1->Project(h1->GetName(),var->GetName(),cut1);
  tree2->Project(h2->GetName(),var->GetName(),cut2);
  
  h1->Scale(100/h1->Integral());
  h2->Scale(100/h2->Integral());
  h1->SetMarkerStyle(8);
  h2->SetMarkerStyle(22);
  h2->SetMarkerColor(kBlue);
  h1->SetMarkerSize(0.6);
  h2->SetMarkerSize(0.6);
  
  Double_t ks = h1->KolmogorovTest(h2);
  TString s;
  s.Form("%-10s%10.4f",var->GetName(),ks);
  os << s <<endl;
  
  plot->addTH1(h1,"pe");
  plot->addTH1(h2,"pe");

  plot->Draw();
  return plot;
}
