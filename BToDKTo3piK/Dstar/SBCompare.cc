// adapted from Frank's code. Does a comparison of backgrounds 
// in data and MC sidebands. Also, does comparison of backgrounds 
// in signal region and in sidebands of Monte Carlo.
   
void SBCompare() {

  for(int i=0; i<4; i++) {
  TCut SBCuts("Dmass>1.93 && Dmass<1.99 && (hh_mass<0.489||hh_mass>0.508)");
  TVector3 r;
  if(i==0) r = dataMCSBCompare(SBCuts, "dataMCSBCompare");
  else if(i==1) r = SignalSBCompare(SBCuts,"SignalSBCompare-CombnRef",0);
  else if(i==2) r = SignalSBCompare(SBCuts,"SignalSBCompare-Comb",1);
  else if(i==3) r = SignalSBCompare(SBCuts,"SignalSBCompare-Ref",2);

  printf ("%40s %10s %10s\n","sideband","chi2/ndof","p-value(%)");
  printf ("%40s %10.3f %10.2f\n",SBCuts.GetTitle(),r[0]/r[1],r[2]*100);
  cout << endl;
  }
}








// Does a comparison of data and MC using "cut" and "qqFrac" as the
// qqbar fraction. Results are drawn either in a new canvas or in
// pad number "pad" of canvas "can"
//
// Returns a 3-vector with (total chi2, ndof, p-value)
TVector3 dataMCSBCompare(TCut cut=0, TString &plotName)
{

  TFile hello("hhPi0.root","READ");
//   TTree* ccTree = (TTree*)hello.Get("pipipi0_ccbar");
  TTree* dataTree = (TTree*)hello.Get("pipipi0_data");

  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  gStyle->SetOptStat(1001110);

  Double_t weight1= 179.018*13.0;
  Double_t weight2= 332.850*5.5;
  Double_t weight3= 215.286*5.5;
  Double_t weight4= 162.902*21.0;
  Double_t sumweight = weight1+weight2+weight3+weight4;
  weight1 /= sumweight;
  weight2 /= sumweight;
  weight3 /= sumweight;
  weight4 /= sumweight;

  TChain chain1("ntp1");
  chain1.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/ccbar*.root");
  chain1.SetWeight(weight1,"global");
  TChain chain2("ntp1");
  chain2.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bbbar*.root");
  chain2.SetWeight(weight2,"global");
  TChain chain3("ntp1");
  chain3.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bplusbminus*.root");
  chain3.SetWeight(weight3,"global");
  TChain chain4("ntp1");
  chain4.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/uds*.root");
  chain4.SetWeight(weight4,"global");

  TChain chain("chain");
  chain.Add( &chain1 );
  chain.Add( &chain2 );
  chain.Add( &chain3 );
  chain.Add( &chain4 );
  TTree* ccTree = (TTree*) &chain;
  TCut chaincut("abs(delta_M-0.1455)<0.0006 && Dmass>1.72 &&  Dmass<1.99");


  TCanvas* can = new TCanvas("can","Dalitz plots",800,600);
  can->Divide(2,2);
  Int_t pad = 1;

  const int MAXBINS = 50;     // Maximum number of bins for Dalitz plots 
  const int MINEVENTS = 15; // minimum number of events per bin

  // Data histogram  
  // determine the number of bins to have reasonable statistics
  // since we have more MC than data, the data histogram defines the limit
  gROOT->cd();
  TTree *tempTree = dataTree->CopyTree(cut );
  Long64_t N = tempTree->GetEntries();  // number of events in total after cuts
  int BINS = (int)sqrt(2*N/MINEVENTS);
  if (BINS>MAXBINS) BINS = MAXBINS;

  
  cout << "Using "<<BINS<<"^2 bins for histograms "
       << "("<<(double)N/(BINS**2/2.0)<<" events per bin on average)"<<endl;



  TH2D *hdata = new TH2D("hdata","",BINS,0,3,BINS,0,3);
  tempTree->Draw("fit_S23:fit_S31>>hdata","","goff");
  hdata->Sumw2();   
  hdata->Scale(1000/hdata->Integral());

  TH2D *hMC = new TH2D("hMC","",BINS,0,3,BINS,0,3);
  ccTree->Draw("fit_S23:fit_S31>>hMC",(cut)&&(chaincut)&&"Flag!=1 && Flag!=2"
	       ,"goff");
  hMC->Sumw2();  
  hMC->Scale(1000/hMC->Integral());

 


 //
  // Plot chi^2
  //
  TH2D *hchi2 = new TH2D("hchi2","",BINS,0,3,BINS,0,3);  
  TH1D* pullhist = new TH1D("pullhist","",25,-4,4); 
  int ndof = hchi2->GetNbinsX() * hchi2->GetNbinsY();
  ndof--;   // histos are normalized to each other

  Double_t chi2Total = 0;    // total chi2
  for (int x=1; x <= hchi2->GetNbinsX(); x++) {
    for (int y=1; y <= hchi2->GetNbinsY(); y++) {
      
      Double_t bin1 =  hdata->GetBinContent(x,y);
      Double_t bin2 = hMC->GetBinContent(x,y); 

      if (bin1==0 || bin2==0) ndof--;    // no data -> one less dof
      else {
        Double_t sqrtChi2 = (bin1-bin2);
        sqrtChi2 /= sqrt(hdata->GetBinError(x,y)**2 
			 + hMC->GetBinError(x,y)**2);        
        hchi2->SetBinContent(x,y,sqrtChi2);
        chi2Total += sqrtChi2**2;
	pullhist->Fill(sqrtChi2);
      }
    }
  }

  can->cd(pad);
  hdata->SetTitle("Data: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hdata->Draw("colz");

  can->cd(pad+1);
  hMC->SetTitle("MC: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hMC->SetMaximum(hdata->GetMaximum());
  hMC->Draw("colz");

 
  can->cd(pad+2);  
  hchi2->SetMaximum(4);
  hchi2->SetMinimum(-4);

  hchi2->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hchi2->Draw("colz");

  can->cd(pad+3);  
  pullhist->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  pullhist->Fit("gaus");
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  pullhist->Draw();

  can->Draw();
  can->SaveAs(plotName + TString(".eps"));
  can->SaveAs(plotName + TString(".gif"));
  can->Close();

  Double_t chi2Prob = TMath::Prob(0.5*chi2Total,Int_t(0.5*ndof));
  cout << "Chi2          = " << chi2Total << endl;
  cout << "ndof          = " << ndof << endl;
  cout << "Chi2/ndof     = " << chi2Total/ndof << endl;
  cout << "Chi2 prob     = " << chi2Prob << endl;


  delete hdata;
  delete hMC;
  delete hchi2;
  delete pullhist;
  delete can;
  return TVector3(chi2Total,ndof,chi2Prob);
}  












// Does a comparison of signal and sideband using ccbar MC.
// Results are drawn either in a new canvas or in
// pad number "pad" of canvas "can"
// Returns a 3-vector with (total chi2, ndof, p-value)
TVector3 SignalSBCompare(TCut cut=0, TString &plotName, Int_t bkgCode=0) {

  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  gStyle->SetOptStat(1001110);

  Double_t weight1= 179.018*13.0;
  Double_t weight2= 332.850*5.5;
  Double_t weight3= 215.286*5.5;
  Double_t weight4= 162.902*21.0;
  Double_t sumweight = weight1+weight2+weight3+weight4;
  weight1 /= sumweight;
  weight2 /= sumweight;
  weight3 /= sumweight;
  weight4 /= sumweight;

  TChain chain1("ntp1");
  chain1.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/ccbar*.root");
  chain1.SetWeight(weight1,"global");
  TChain chain2("ntp1");
  chain2.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bbbar*.root");
  chain2.SetWeight(weight2,"global");
  TChain chain3("ntp1");
  chain3.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/bplusbminus*.root");
  chain3.SetWeight(weight3,"global");
  TChain chain4("ntp1");
  chain4.Add("/nfs/farm/babar/AWG33/Charm/kalanand/PiPiPi0/uds*.root");
  chain4.SetWeight(weight4,"global");

  TChain chain("chain");
  chain.Add( &chain1 );
  chain.Add( &chain2 );
  chain.Add( &chain3 );
  chain.Add( &chain4 );
  TTree* ccTree = (TTree*) &chain;
  TCut chaincut("abs(delta_M-0.1455)<0.0006 && Dmass>1.72 &&  Dmass<1.99");

  TCanvas* can = new TCanvas("can","Dalitz plots",800,600);
  can->Divide(2,2);
  Int_t pad = 1;


  const int MAXBINS = 50;     // Maximum number of bins for Dalitz plots 
  const int MINEVENTS = 15; // minimum number of events per bin

  // Sideband histogram  
  // determine the number of bins to have reasonable statistics
  // since we have more signal events than sideband events, the sideband
  //  histogram defines the limit



  TString flagstring;
  /*
    bkgCode==0 : no signal events  => comb. background + KPiPi0 reflection
    bkgCode==1 : no signal events, ro reflection events => only comb. background
    bkgCode==2 : KPiPi0 reflection events only

  */
  if(bkgCode==0) flagstring="!(Flag==1||Flag==2||Flag==11||Flag==101)";
  else if(bkgCode==1) flagstring="!(Flag==1||Flag==2||Flag==11||Flag==101||Flag==10||Flag==20)";
  else if(bkgCode==2) flagstring="Flag==10||Flag==20";

  TCut flagCut(flagstring);




  gROOT->cd();
  TTree *tempTree = ccTree->CopyTree( (cut)&&(chaincut)&&(flagCut) );

  Long64_t N = tempTree->GetEntries();  // number of events in total after cuts
  int BINS = (int)sqrt(2*N/MINEVENTS);
  if (BINS>MAXBINS) BINS = MAXBINS;


  TH2D *hsig = new TH2D("hsig","",BINS,0,3,BINS,0,3);
  ccTree->Draw("fit_S23:fit_S31>>hsig",
	       (chaincut)&&(flagCut)&&"abs(Dmass-1.8605)<0.0167 && (hh_mass<0.489 || hh_mass>0.508)",
	       "goff");
  hsig->Sumw2();   
  hsig->Scale(1000/hsig->Integral());

  TH2D *hSB = new TH2D("hSB","",BINS,0,3,BINS,0,3);
  tempTree->Draw("fit_S23:fit_S31>>hSB",(cut)&&(chaincut)&&(flagCut) ,"goff");
  hSB->Sumw2();  
  hSB->Scale(1000/hSB->Integral());

 


 //
  // Plot chi^2
  //
  TH2D *hchi2 = new TH2D("hchi2","",BINS,0,3,BINS,0,3);  
  TH1D* pullhist = new TH1D("pullhist","",20,-4,4); 
  int ndof = hchi2->GetNbinsX() * hchi2->GetNbinsY();
  ndof--;   // histos are normalized to each other

  Double_t chi2Total = 0;    // total chi2
  for (int x=1; x <= hchi2->GetNbinsX(); x++) {
    for (int y=1; y <= hchi2->GetNbinsY(); y++) {
      
      Double_t bin1 =  hsig->GetBinContent(x,y);
      Double_t bin2 = hSB->GetBinContent(x,y); 

      if (bin1==0 || bin2==0) ndof--;    // no data -> one less dof
      else {
        Double_t sqrtChi2 = (bin1-bin2);
	sqrtChi2 /= sqrt(hsig->GetBinError(x,y)**2 + hSB->GetBinError(x,y)**2);
        hchi2->SetBinContent(x,y,sqrtChi2);
        chi2Total += sqrtChi2**2;
	pullhist->Fill(sqrtChi2);
      }
    }
  }

  can->cd(pad);
  hsig->SetTitle("Signal: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hsig->Draw("colz");

  can->cd(pad+1);
  hSB->SetTitle("Sideband: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hSB->SetMaximum(hsig->GetMaximum());
  hSB->Draw("colz");
 
  can->cd(pad+2);  
  hchi2->SetMaximum(4);
  hchi2->SetMinimum(-4);
  hchi2->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  hchi2->Draw("colz");

  can->cd(pad+3);  
  pullhist->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
  pullhist->Fit("gaus");
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.9);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  pullhist->Draw();

  can->Draw();
  can->SaveAs(plotName + TString(".eps"));
  can->SaveAs(plotName + TString(".gif"));
  can->Close();

  Double_t chi2Prob = TMath::Prob(0.5*chi2Total,Int_t(0.5*ndof));
  cout << "Chi2          = " << chi2Total << endl;
  cout << "ndof          = " << ndof << endl;
  cout << "Chi2/ndof     = " << chi2Total/ndof << endl;
  cout << "Chi2 prob     = " << chi2Prob << endl;

  delete hsig;
  delete hSB;
  delete hchi2;
  delete pullhist;
  delete can;
  return TVector3(chi2Total,ndof,chi2Prob);
}  








