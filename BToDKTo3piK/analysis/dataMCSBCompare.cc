// $Id: dataMCSBCompare.cc,v 1.14 2006/05/20 02:50:54 fwinkl Exp $
// Sideband studies. Compares data to weighted BBbar and qqbar MC

int MINEVENTS = 30;   // minimum number of events per bin

// Does the comparison for all sidebands and collects the results
void MCSBCompareDalitzAll(TCut addCut = cutBasic) 
{
  TCanvas *mainCavas = new TCanvas("mainCanvas","Dalitz plots",600,1000);
  mainCanvas->Divide(3,4);

  const int nnCuts_n = 4;
  const int minEvents_n = 2;

  TCut nnCuts[nnCuts_n] = {"","nnout<0.25","nnout>0.25","nnout>0.6"};
  int minEvents[minEvents_n] = {15,30};
  
  TVector3 r1[nnCuts_n][minEvents_n];
  TVector3 r2[nnCuts_n][minEvents_n];
  TVector3 r3[nnCuts_n][minEvents_n];
  TVector3 r4[nnCuts_n][minEvents_n];
  
  for (int i = 0; i<nnCuts_n; i++) {
    for (int j = 0; j<minEvents_n; j++) {
      MINEVENTS = minEvents[j];

      r1[i][j] = dataMCSBCompareDalitz(cutSBUpperDE && addCut && nnCuts[i],
                                       mainCanvas,1);
      r2[i][j] = dataMCSBCompareDalitz(cutSBLowerDE && addCut && nnCuts[i],
                                       mainCanvas,4);  

      r3[i][j] = dataMCSBCompareDalitz(cutSBmES && addCut && nnCuts[i],
                                       mainCanvas,7);  
      r4[i][j] = dataMCSBCompareDalitz(cutSBMD && addCut && nnCuts[i],
                                       mainCanvas,10);
    }
  }
  mainCanvas->Draw();
  /*
  mainCanvas->SaveAs("dataMCSBCompare.root");
  mainCanvas->SaveAs("dataMCSBCompare.eps");
  mainCanvas->SaveAs("dataMCSBCompare.png");
  */

  for (int i = 0; i<nnCuts_n; i++) {
    for (int j = 0; j<minEvents_n; j++) { 
      cout << nnCuts[i].GetTitle() << " minEvents = "<<minEvents[j]<<endl;
      printf ("%15s %10s %10s\n","sideband","chi2/ndof","p-value(%)");
      printf ("%15s %10.3f %10.2f\n","upper DE",r1[i][j][0]/r1[i][j][1],r1[i][j][2]*100);
      printf ("%15s %10.3f %10.2f\n","lower DE",r2[i][j][0]/r2[i][j][1],r2[i][j][2]*100);
      printf ("%15s %10.3f %10.2f\n","mES",r3[i][j][0]/r3[i][j][1],r3[i][j][2]*100);
      printf ("%15s %10.3f %10.2f\n","mD",r4[i][j][0]/r4[i][j][1],r4[i][j][2]*100);
      cout << endl;
    }
  }

}


// Does a comparison of data and MC Dalitz plots using "cut".
// Results are drawn either in a new canvas or in
// pad number "pad" of canvas "can"
//
// Returns a 3-vector with (total chi2, ndof, p-value)
TVector3 dataMCSBCompareDalitz(TCut cut = 0,
                               TCanvas *can = 0, Int_t pad = 1)
{
  gROOT->cd();
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  gStyle->SetOptStat(1001110);

  // Create new canavs if not supplied
  if (!can) {
    can = new TCanvas("can2","can2",900,300);
    can->Divide(3,1);
  }

  const int MAXBINS = 30;     // Maximum number of bins for Dalitz plots 

  // Data histogram  
  // determine the number of bins to have reasonable statistics
  // since we have more MC than data, the data histogram defines the limit
  TTree *tempTree = dataTree->CopyTree(cut);
  Long64_t N = tempTree->GetEntries();  // number of events in total after cuts

  if (N==0) {
    cout << "No entries in "<<dataTree->GetName()<<" after cut "<<
         << cut->GetTitle() << endl;
    return;
  }

  int BINS = (int)sqrt(2*N/MINEVENTS);
  if (BINS>MAXBINS) BINS = MAXBINS;
  
  cout << "Using "<<BINS<<"^2 bins for histograms "
       << "("<<(double)N/(BINS**2/2.0)<<" events per bin on average)"<<endl;

  TH2D *hdata = new TH2D("hdata","",BINS,0,3,BINS,0,3);
  tempTree->Draw("d0pppupmass**2:d0ppmupmass**2>>hdata","","goff");
  hdata->Sumw2();    // make sure the errors are scaled correctly
  hdata->Scale(1000/hdata->Integral());

  delete tempTree;

  TH2D *hMC = new TH2D("hMC","",BINS,0,3,BINS,0,3);
  hMC->Sumw2();
  TString vars("d0pppupmass**2:d0ppmupmass**2");
  weightedMCHisto(hMC,vars,cut);
  hMC->Scale(1000/hMC->Integral());

  //
  // Plot chi^2
  //
  TH2D *hchi2 = new TH2D("hchi2","",BINS,0,3,BINS,0,3);  
  TVector3 chi2test = chi2test2d(hdata,hMC,hchi2);

  can->cd(pad);
  hdata->SetTitle(TString("Data (")+TString(cut.GetName())+");s12;s13");
  hdata->Draw("colz");

  can->cd(pad+1);
  hMC->SetTitle(TString("MC (")+TString(cut.GetName())+");s12;s13");
  hMC->SetMaximum(hdata->GetMaximum());
  hMC->Draw("colz");

 
  can->cd(pad+2);  
  //  hchi2->SetMinimum(-3);
  hchi2->SetMaximum(3);

  //  if (hchi2->GetMaximum() > 10) hchi2->SetMaximum(10);
  hchi2->SetTitle(TString("#sqrt{chi2} (")+TString(cut.GetName())+");s12;s13");
  hchi2->Draw("colz");

  //  Double_t chi2Prob = TMath::Prob(0.5*chi2Total,Int_t(0.5*ndof));
  cout << "Chi2          = " << chi2test[0] << endl;
  cout << "ndof          = " << chi2test[1] << endl;
  cout << "Chi2/ndof     = " << chi2test[0]/chi2test[1] << endl;
  cout << "Chi2 prob     = " << chi2test[2] << endl;

  return chi2test;
}  


// Compare all 1D variables
void dataMCSBCompareAll()
{
  const int nnCuts_n = 4;
  TCut nnCuts[nnCuts_n] = {"","nnout<0.25","nnout>0.25","nnout>0.6"};

  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.3,"Y");


  TString s;
  Double_t ks = -1;
  for (int i=0; i<nnCuts_n; i++) {
    TString result;
    TCanvas *can = new TCanvas("can","",1200,400);
    can->Divide(3,1,0.004,0.004);

    can->cd(1);
    ks = dataMCCompare(*Deltae, cutSBmES+cutBasic+nnCuts[i]);
    s.Form("%.3f\n",ks); result += s;
    can->cd(2);
    ks = dataMCCompare(*Deltae, cutSBLowerMD+cutBasic+nnCuts[i]);
    s.Form("%.3f\n",ks); result += s;
    can->cd(3);
    ks = dataMCCompare(*Deltae, cutSBUpperMD+cutBasic+nnCuts[i]);
    s.Form("%.3f\n",ks); result += s;

    cout << "-----------------------------------------------------"<<endl;
    cout << "Chi2 prob. "<<nnCuts[i].GetTitle()<<endl<<result;
    cout << "-----------------------------------------------------"<<endl;

    TString filename("dataMCSBCompare-Deltae_");
    filename += i;
    can->SaveAs(filename+".eps");
    can->SaveAs(filename+".root");
  }

  TString result;
  TCanvas *can = new TCanvas("can","",1200,800);
  can->Divide(3,2,0.004,0.004);

  can->cd(1);
  ks = dataMCCompare(*nnout, cutSBmES+cutBasic);
  s.Form("%.3f\n",ks); result += s;
  can->cd(2);
  ks = dataMCCompare(*nnout, cutSBLowerMD+cutBasic);
  s.Form("%.3f\n",ks); result += s;
  can->cd(3);
  ks = dataMCCompare(*nnout, cutSBUpperMD+cutBasic);
  s.Form("%.3f\n",ks); result += s;
  can->cd(4);
  ks = dataMCCompare(*nnout, cutSBLowerDE+cutBasic);
  s.Form("%.3f\n",ks); result += s;
  can->cd(5);
  ks = dataMCCompare(*nnout, cutSBUpperDE+cutBasic);
  s.Form("%.3f\n",ks); result += s;
    
  cout << "-----------------------------------------------------"<<endl;
  cout << "Chi2 prob. "<<endl<<result;
  cout << "-----------------------------------------------------"<<endl;
  
  TString filename("dataMCSBCompare-nnout");
  can->SaveAs(filename+".eps");
  can->SaveAs(filename+".root");
}


void dataMCCompareAll(Bool_t addSignal = kTRUE, Bool_t normalize = kTRUE)
{
  TCut cut("",cutDeltaE+cutmES+cutMD+cutNNq+cutDtoKpi+cutKsVeto+"bknnout<0.5");

  RooArgList vars(*Deltae,*mes,*nnout,*bknnout,*d0mass,*mass12);
  TArrayD pchi2(vars.getSize());

  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.3,"Y");  
  TCanvas* can = new TCanvas("can","",900,600);
  can->Divide(3,2);

  for (int i=0; i<vars.getSize(); i++) {
    can->cd(i+1);
    Double_t p = dataMCCompare(*(RooRealVar*)vars.at(i),cut,addSignal,normalize);
    pchi2.AddAt(p,i);
    can->Update();
  }

  cout << "Chi2 probabilities:"<<endl;
  for (int i=0; i<vars.getSize(); i++) {
    cout << setw(15) << vars.at(i)->GetName() << ": "
         << setw(10) << setprecision(4) << pchi2[i] << endl;
  }
}

// Plot and compare data and MC for var with cut and return chi2 prob.
Double_t dataMCCompare(RooRealVar& var, TCut cut, 
                       Bool_t addSignal = kFALSE, Bool_t normalize = kTRUE)
{
  TH1 *hdata = var.createHistogram("hdata");
  hdata->Sumw2();
  hdata->SetTitle(cut.GetName());
  hdata->SetMarkerStyle(22);
  hdata->SetMarkerSize(0.6);
  hdata->SetMinimum(0);
  dataTree->Project(hdata->GetName(),var.GetName(),cut);
  
  TH1 *hMC = var.createHistogram("hMC");
  weightedMCHisto(hMC,var.GetName(),cut,addSignal);
  hMC->SetMinimum(0);
  hMC->SetMarkerStyle(20);
  hMC->SetMarkerSize(0.6);
  hMC->SetMarkerColor(kBlue);
  hMC->SetLineColor(kBlue);

  Double_t pchi2 = -1;
  if (normalize) {
    hdata->DrawNormalized("pe",100);
    hMC->DrawNormalized("pe same",100);
    pchi2 = hMC->Chi2Test(hdata,"",1);   // 1 constraint
  }
  else {
    hdata->Draw("pe");
    hMC->Draw("pe same");
    pchi2 = hMC->Chi2Test(hdata,"",0);   // 0 constraints
  }

  return pchi2;
}


//-------------------------------------------------------------------
// Returns the fraction of qq events in data using an R2 fit:
double findQQFrac(const RooAbsData * data,
		  const RooAbsData * mcBB,
		  const RooAbsData * mcQQ) {

  // Set global fit option (used by fit macro)
  fitOption = "mr";

  // Start by setting up a Gaussian and a bifurcated Gaussian PDF for qq:
  BdkPdfBifurGauss bifurQq("bifurQq", "bifurQq", *R2);
  bifurQq.mean()->setVal(0.3);
  bifurQq.sigl()->setVal(0.1);
  bifurQq.sigr()->setVal(0.1);

  // and one for BB:
  BdkPdfBifurGauss bifurBb("bifurBb", "bifurBb", *R2);
  bifurBb.mean()->setVal(0.15);
  bifurBb.sigl()->setVal(0.1);
  bifurBb.sigr()->setVal(0.1);

  // And their sum:
  BdkPdfSum r2Pdf("r2Pdf", "r2Pdf", bifurBb, bifurQq);
  
  TCanvas * can = new TCanvas("qqFracCan", "qqFracCan", 1200, 300);
  can->Divide(3,1);

  // Fit the individual MC samples:
  fit(bifurQq, *mcQQ);
  RooPlot * plQQ = R2->frame();
  mcQQ->plotOn(plQQ);
  bifurQq.getPdf()->plotOn(plQQ);
  can->cd(1);
  plQQ->Draw();

  fit(bifurBb, *mcBB);
  RooPlot * plBB = R2->frame();
  mcBB->plotOn(plBB);
  bifurBb.getPdf()->plotOn(plBB);
  can->cd(2);
  plBB->Draw();

  // Fix all pars except the BB fraction and fit the data:
  r2Pdf.fixAll();
  r2Pdf.fixFractions(kFALSE);
  RooFitResult * result = fit(r2Pdf, *data);
  result->Print("V");

  RooPlot *plData = R2->frame();
  data->plotOn(plData);
  r2Pdf.getPdf()->plotOn(plData);

  RooArgSet setQq(*bifurQq.getPdf());
  r2Pdf.getPdf()->plotCompOn(plData, setQq);

  RooArgSet setBb(*bifurBb.getPdf());
  r2Pdf.getPdf()->plotCompOn(plData, setBb);

  can->cd(3);
  plData->Draw();

  return (1-r2Pdf.fraction(0)->getVal());
}
