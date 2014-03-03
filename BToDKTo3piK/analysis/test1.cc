#include "../BToDKTo3piK/globals/setup.cc"
#include "../BToDKTo3piK/utils/utils.cc"

void test1() {
  // Define files to read:
  //  TChain chain("h1");
  //  chain.Add("bch1.root");
  //
  
  // Read them:
  data = read(chainDK);

  // Fix everything and then unfix just what we want floating:
  pdfOnResDK.fixAll();

  RooRealVar * parToFloat = 
    (RooRealVar*)paramsOnResDK.find("pdfOnResDK.sigGoodD0NumEvts");

  parToFloat->setConstant(kFALSE);

  // Fit:
  cout<<"Fitting the variables"<<endl;
  pdfOnResDK.parametersFree().Print("V");
  RooFitResult * result = fit(pdfOnResDK, *data);
  result->Print();

  // Plot:  
  data->plotOn(mesFrame);
  pdfOnResDK.getPdf()->plotOn(mesFrame);

  data->plotOn(DeltaeFrame);
  pdfOnResDK.getPdf()->plotOn(DeltaeFrame);

  data->plotOn(d0massFrame);
  pdfOnResDK.getPdf()->plotOn(d0massFrame);

  TCanvas *c1 = new TCanvas("c1","test bch",900,600);
  c1->Divide(3,1);

  c1->cd(1);
  mesFrame->Draw();
  c1->cd(2);
  DeltaeFrame->Draw();
  c1->cd(3);
  d0massFrame->Draw();
}
