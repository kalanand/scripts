// makes plots of all the PDF components 

void plotAllPdfs(int nEvents = 10000,
		 Bool_t doDalitz = kTRUE, 
		 Bool_t doQ = kFALSE, 
		 Bool_t doDE = kFALSE) {

  TFile file("plotAllPdfs.root", "RECREATE");
	
  for (int t = 0; t < BdkEvtTypes::NTYPES; ++t) {
    cout << "Making hists for " << BdkEvtTypes::name(t) << endl;

    BdkPdfProdAll * pdfN = pdfOnResDK.prodN(t);
    BdkPdfProdAll * pdfP = pdfOnResDK.prodP(t);
    
    RooDataSet * dataN = pdfN->generate(nEvents);
    RooDataSet * dataP = pdfP->generate(nEvents);
    
    TH2 * histN = dataN->createHistogram(*m12, *m13);
    histN->SetName(TString(BdkEvtTypes::name(t)) + ".N");
    histN->SetTitle(TString(BdkEvtTypes::name(t)) + ".N");

    TH2 * histP = dataP->createHistogram(*m12, *m13);
    histP->SetName(TString(BdkEvtTypes::name(t)) + ".P");
    histP->SetTitle(TString(BdkEvtTypes::name(t)) + ".P");

    histN->Write();
    histP->Write();
  }

  cout << "histograms written to " << file.GetName() << endl;
  file.Close();
}


void readPlotAllPdfsFile() {
  TFile * file = new TFile("plotAllPdfs.root");
  for (int t = 0; t < BdkEvtTypes::NTYPES; ++t) {
    TH2F * histN = (TH2F*)(file->Get(TString(BdkEvtTypes::name(t)) + ".N"));
    TH2F * histP = (TH2F*)(file->Get(TString(BdkEvtTypes::name(t)) + ".P"));

    TCanvas * can = new TCanvas(BdkEvtTypes::name(t), BdkEvtTypes::name(t), 
				1000, 500);
    can->Divide(2,1);
    can->cd(1);
    histN->Draw();
    can->cd(2);
    histP->Draw();
    can->Update();
  }
}  
  
