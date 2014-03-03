// $Id: getNNhists.cc,v 1.2 2006/05/19 17:34:33 fwinkl Exp $
// Save the histograms for the NnCont PDF to a file


// Get all histograms for var
void getNNhistsAll(RooRealVar* var) 
{
  
  TH1D* h;
  TList list;

  list.Add(getHist(sigTree, cutDKGoodD, var));
  list.Add(getHist(sigTree, cutDKBadD, var));
  list.Add(getHist(dpiTree, cutDPiBadD, var));
  list.Add(getHist(dpiTree, cutDPiGoodD, var));
  list.Add(getHist(bbTree, cutDKX, var));
  list.Add(getHist(bbTree, cutDPiX, var));
  list.Add(getHist(bbTree, cutBBBadD, var));
  list.Add(getHist(bbTree, cutBBGoodD, var));
  list.Add(getHist(qqTree, cutqqGoodD, var));
  list.Add(getHist(qqTree, cutqqBadD, var));
    
  list.Print();

  TCanvas *can = new TCanvas("can",var->GetName(),1000,500);
  can->Divide(5,2);
  gStyle->SetOptStat(0);
  
  TFile f(TString(var->GetName())+".root","RECREATE");
  for (int i=0; i<list.GetEntries(); i++) {    
    ((TH1D*)list.At(i))->Write();
    can->cd(i+1);
    ((TH1D*)list.At(i))->Draw();
  }
  cout << "Histograms saved to "<<f.GetName()<<endl;
  f.Close();  
  can->SaveAs(TString(var->GetName())+".eps");
}


//-------------------------------------------------------------
// Does common tasks for all variables
// The signal region cut is added to cut
TH1D* getHist(TTree *tree, TCut cut, RooRealVar *var)
{

  const Int_t BINS = var->getBins();
  const char* title = cut.GetName();
  
  if (tree->GetEntries() == 0) {
    cout << "empty tree" << endl;
    return 0;
  }

  cout << "Getting histogram for "<<title<<endl;

  TH1D* h = new TH1D(title,title,BINS,var->getMin(),var->getMax());
  h->GetXaxis()->SetTitle(var->GetTitle());
  tree->Project(h->GetName(),var->GetName(),cut+cutSigReg,"goff");

  return h;
}
  
