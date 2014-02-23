// $Id: badSigShape.cc,v 1.2 2006/07/11 21:44:01 fwinkl Exp $
// Get alternative Dalitz histograms for DKBadD

void badSigShape()
{
  TCut cut = cutSigReg;

  gROOT->cd();

  // Get bad signal from flat MC and split into B+/B-
  TTree* flatTree = sigFlatTree->CopyTree(cut+cutDKBadD);
  RooDataSet* flatBpData = read(flatTree,0,"Hdtrkchge>0",allVars);
  RooDataSet* flatBmData = read(flatTree,0,"Hdtrkchge<0",allVars);
  
  data = shapeBpBm(flatBpData, flatBmData);


  getDalitzParamsDKBadD(data);

}


RooDataSet* shapeBpBm(RooDataSet* flatBpData, RooDataSet* flatBmData)
{
  // Reshape these event with good(!) signal PDF
  RooDataSet* badBpData = generateFromFlat(*flatBpData,
                                           (BdkDKDalitz*)dalitzHolderP.sigGoodD0()->getPdf());
  RooDataSet* badBmData = generateFromFlat(*flatBmData,
                                           (BdkDKDalitz*)dalitzHolderN.sigGoodD0()->getPdf());

  Double_t effBp = (Double_t)badBpData->numEntries()/flatBpData->numEntries();
  Double_t effBm = (Double_t)badBmData->numEntries()/flatBmData->numEntries();
  
  cout << "Reweighting efficencies: B+ = "<<effBp<<endl
       << "                         B- = "<<effBm<<endl;    
 

  badBpData->append(*badBmData);
  return badBpData;
}


// Define the DKBadD histogram binning
void getDalitzParamsDKBadD(RooDataSet* data, TString postfix = "")
{
  // Define binning
  RooBinning bins12, bins13;
  bins12.setMin(0); bins12.setMax(3);
  bins13.setMin(0); bins13.setMax(3);

  // We use the same binning for m12 and m13
  Double_t d12[] = {0.25,0.5,0.75,1.6,
		    2.15,2.30,2.45,2.60,2.75};
  
  //  Double_t d13[] = {0.24,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
 
  setBinning(bins12,d12,sizeof(d12)/sizeof(Double_t));
  setBinning(bins13,d12,sizeof(d12)/sizeof(Double_t));

  TH2 *hist = getDalitzHist(&data->tree(),"",bins12,bins13);
  hist->SetName("hist_dkbadd");

  // Remove the palette since that gives warnings on reading the root file
  TList *l = hist->GetListOfFunctions();
  l->Remove(l->FindObject("palette"));

  const char* file = "hist_dkbadd"+postfix+".root";
  cout << "Writing DKBadD histogram PDF to "<<file<<endl;
  TFile f(file,"recreate");
  hist->Write();
  f.Close();
}
