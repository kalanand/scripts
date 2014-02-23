// $Id: compareR16R18.cc,v 1.5 2006/04/25 17:51:47 fwinkl Exp $
// Script to compare R16 and R18 data

// if this is true, compareAll() uses the flat signal MC
Bool_t useFlatSignalMC = kFALSE;

void compareR16R18()
{
  // Make clones so that changes don't affect the original vars
  RooRealVar mynnout(*nnout);
  RooRealVar mybknnout(*bknnout);

  mynnout.setRange(0,1); 
  mynnout.setBins(20);
  mybknnout.setRange(0,1);
  mybknnout.setBins(20);

  compareAll(cutSigReg, Deltae);
  compareAll(cutSigReg, mes);
  compareAll(cutSigReg, d0mass);
  compareAll(cutDeltaE+cutmES+cutMD+cutDtoKpi+cutKsVeto+TCut("bknnout>0.25"), &mynnout);
  compareAll(cutDeltaE+cutmES+cutMD+cutDtoKpi+cutKsVeto+TCut("nnout>0.1"), &mybknnout);
  compareAll(cutSigReg, m12);   // Dalitz plot
}


// Compare all event types for var
void compareAll(TCut cut, RooRealVar* var)
{

  TCanvas * can = new TCanvas("can", "compareR16R18", 1200, 600);
  can->Divide(4,2,0.004,0.004);
  TGaxis::SetMaxDigits(3);
  can->SetLeftMargin(0.12);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetOptStat(0);

  ostringstream os;
  os << "KS-"<<var->GetName()<<":"<<endl;

  int i = 1;
  can->cd(i);
  if (useFlatSignalMC)
    compare(sigFlatTree16, sigFlatTree18, cutDKGoodD+cut, var, os)->Draw();
  else
    compare(sigTree16, sigTree18, cutDKGoodD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  if (useFlatSignalMC)
    compare(sigFlatTree16, sigFlatTree18, cutDKBadD+cut, var, os)->Draw();
  else
    compare(sigTree16, sigTree18, cutDKBadD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(dpiTree16, dpiTree18, cutDPiGoodD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(dpiTree16, dpiTree18, cutDPiBadD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(bbTree16, bbTree18, cutDKX+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(bbTree16, bbTree18, cutDPiX+cut, var, os)->Draw();
  //  can->cd(++i); can->Update();
  //  compare(bbTree16, bbTree18, cutBBGoodD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(bbTree16, bbTree18, cutBBBadD+cut, var, os)->Draw();
  //  can->cd(++i); can->Update();
  //  compare(qqTree16, qqTree18, cutqqGoodD+cut, var, os)->Draw();
  can->cd(++i); can->Update();
  compare(qqTree16, qqTree18, cutqqBadD+cut, var, os)->Draw();

  TString file = "compareR16R18-";
  file += var->GetName();

  cout << os.str();
  cout << "KS probabilities saved to "<<file<<".txt"<<endl;
  ofstream of(file+".txt");
  of << os.str();
  of.close();

  can->SaveAs(file+".eps");
  can->SaveAs(file+".root");
}


// Compare event type defined by cut for var
RooPlot* compare(TTree* tree16, TTree* tree18,
                 TCut cut, RooRealVar* var, ostringstream& os)
{
  const char* title = cut.GetName();

  RooPlot* plot = 0;

  // Dalitz
  if (var==m12 || var==m13) {
    plot = compareDalitz(tree16,tree18,cut,os);
  }
  // other 1D variables
  else {
    plot = var->frame();
    plot->SetName(title);
    plot->SetTitle(title);

    TH1* h16 = var->createHistogram("r16");
    TH1* h18 = var->createHistogram("r18");
    h16->Sumw2();
    h18->Sumw2();

    tree16->Project(h16->GetName(),var->GetName(),cut);
    tree18->Project(h18->GetName(),var->GetName(),cut);

    h16->Scale(100/h16->Integral());
    h18->Scale(100/h18->Integral());
    h16->SetMarkerStyle(8);
    h18->SetMarkerStyle(22);
    h18->SetMarkerColor(kBlue);
    h16->SetMarkerSize(0.6);
    h18->SetMarkerSize(0.6);
    
    Double_t ks = h16->KolmogorovTest(h18);
    TString s;
    s.Form("%-10s%10.4f",title,ks);
    os << s <<endl;

    plot->addTH1(h16,"pe");
    plot->addTH1(h18,"pe");
  }

  return plot;
}


// Compare two Dalitz plots
RooPlot* compareDalitz(TTree* tree16, TTree* tree18,
                       TCut cut, ostringstream& os)
{
  const Int_t BINS = 30;
  const char* title = cut.GetName();

  RooDataSet* data16 = read(tree16,0,cut,allVars);
  RooDataSet* data18 = read(tree18,0,cut,allVars);

  TH2F* h16 = data16.createHistogram(*m12,*m13,BINS,BINS);
  TH2F* h18 = data18.createHistogram(*m12,*m13,BINS,BINS);

  TH2F *hchi2 = new TH2F("hchi2",title,BINS,0,3,BINS,0,3);  
  TVector3 chi2test = chi2test2d(h16,h18,hchi2);

  RooPlot* p = new RooPlot(*m12,*m13);
  p->addObject(hchi2,"colz");

  Double_t ks = h16->KolmogorovTest(h18);
  TString s;
  s.Form("%-10s%10.4f",title,ks);
  os << s <<endl;

  delete data16;
  delete data18;
  delete h16;
  delete h18;
  return p;
}
