// makes a TGraph of nll values of a number of DP 
// fits with done different starting points.

void nllPlot(){
  Float_t num[64] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
		   11,12,13,14,15,16,17,18,19,20,
		   21,22,23,24,25,26,27,28,29,30,
		   31,32,33,34,35,36,37,38,39,40,
		   41,42,43,44,45,46,47,48,49,50,
		   51,52,53,54,55,56,57,58,59,60,
		   61,62,63,64};

  Float_t nll[64];
  TFile f21("~/public_html/forAbi/PhaseScan/DP_PhaseScan.root","READ");

  nll[0] = getName(f21, 0, 0, 0, 0, 0, 0);
  nll[1] = getName(f21, 1, 1, 1, 1, 1, 1);

  nll[2] = getName(f21, 1, 0, 0, 0, 0, 0);
  nll[3] = getName(f21, 0, 1, 0, 0, 0, 0);
  nll[4] = getName(f21, 0, 0, 1, 0, 0, 0);
  nll[5] = getName(f21, 0, 0, 0, 1, 0, 0);
  nll[6] = getName(f21, 0, 0, 0, 0, 1, 0);
  nll[7] = getName(f21, 0, 0, 0, 0, 0, 1);

  nll[8] = getName(f21, 1, 1, 0, 0, 0, 0);
  nll[9] = getName(f21, 0, 1, 1, 0, 0, 0);
  nll[10] = getName(f21, 0, 0, 1, 1, 0, 0);
  nll[11] = getName(f21, 0, 0, 0, 1, 1, 0);
  nll[12] = getName(f21, 0, 0, 0, 0, 1, 1);
  nll[13] = getName(f21, 1, 0, 0, 0, 0, 1);
  nll[14] = getName(f21, 1, 0, 1, 0, 0, 0);
  nll[15] = getName(f21, 0, 1, 0, 1, 0, 0);
  nll[16] = getName(f21, 0, 0, 1, 0, 1, 0);
  nll[17] = getName(f21, 0, 0, 0, 1, 0, 1);
  nll[18] = getName(f21, 1, 0, 0, 0, 1, 0);
  nll[19] = getName(f21, 0, 1, 0, 0, 0, 1);
  nll[20] = getName(f21, 1, 0, 0, 1, 0, 0);
  nll[21] = getName(f21, 0, 1, 0, 0, 1, 0);
  nll[22] = getName(f21, 0, 0, 1, 0, 0, 1);

  nll[23] = getName(f21, 1, 1, 1, 0, 0, 0);
  nll[24] = getName(f21, 0, 1, 1, 1, 0, 0);
  nll[25] = getName(f21, 0, 0, 1, 1, 1, 0);
  nll[26] = getName(f21, 0, 0, 0, 1, 1, 1);
  nll[27] = getName(f21, 1, 0, 0, 0, 1, 1);
  nll[28] = getName(f21, 1, 1, 0, 0, 0, 1);
  nll[29] = getName(f21, 1, 0, 1, 1, 0, 0);
  nll[30] = getName(f21, 0, 1, 0, 1, 1, 0);
  nll[31] = getName(f21, 0, 0, 1, 0, 1, 1);
  nll[32] = getName(f21, 1, 0, 0, 1, 0, 1);
  nll[33] = getName(f21, 1, 1, 0, 0, 1, 0);
  nll[34] = getName(f21, 0, 1, 1, 0, 0, 1);
  nll[35] = getName(f21, 1, 0, 0, 1, 1, 0);
  nll[36] = getName(f21, 0, 1, 0, 0, 1, 1);
  nll[37] = getName(f21, 1, 0, 1, 0, 0, 1);
  nll[38] = getName(f21, 1, 1, 0, 1, 0, 0);
  nll[39] = getName(f21, 0, 1, 1, 0, 1, 0);
  nll[40] = getName(f21, 0, 0, 1, 1, 0, 1);
  nll[41] = getName(f21, 0, 0, 0, 1, 1, 1);
  nll[42] = getName(f21, 1, 1, 0, 0, 0, 1);

  nll[43] = getName(f21, 1, 1, 1, 1, 0, 0);
  nll[44] = getName(f21, 0, 1, 1, 1, 1, 0);
  nll[45] = getName(f21, 0, 0, 1, 1, 1, 1);
  nll[46] = getName(f21, 1, 0, 0, 1, 1, 1);
  nll[47] = getName(f21, 1, 1, 0, 0, 1, 1);
  nll[48] = getName(f21, 1, 1, 1, 0, 0, 1);
  nll[49] = getName(f21, 1, 0, 1, 1, 1, 0);
  nll[50] = getName(f21, 0, 1, 0, 1, 1, 1);
  nll[51] = getName(f21, 1, 0, 1, 0, 1, 1);
  nll[52] = getName(f21, 1, 1, 0, 1, 0, 1);
  nll[53] = getName(f21, 1, 1, 1, 0, 1, 0);
  nll[54] = getName(f21, 0, 1, 1, 1, 0, 1);
  nll[55] = getName(f21, 1, 0, 1, 1, 0, 1);
  nll[56] = getName(f21, 1, 1, 0, 1, 1, 0);
  nll[57] = getName(f21, 0, 1, 1, 0, 1, 1);

  nll[58] = getName(f21, 0, 1, 1, 1, 1, 1);
  nll[59] = getName(f21, 1, 0, 1, 1, 1, 1);
  nll[60] = getName(f21, 1, 1, 0, 1, 1, 1);
  nll[61] = getName(f21, 1, 1, 1, 0, 1, 1);
  nll[62] = getName(f21, 1, 1, 1, 1, 0, 1);
  nll[63] = getName(f21, 1, 1, 1, 1, 1, 0);

TGraph* tgint = new TGraph(64, num, nll);
TCanvas* tcint = new TCanvas("tcint","");
tgint->SetMarkerColor(2);
tgint->SetMarkerStyle(21);
tgint->SetTitle("Plot of FCN vs Fit No.");
//tgint->GetXaxis()->SetLimits(0.5, 1.02);
tgint->GetXaxis()->SetTitle("Fit No.");
tgint->Draw("AP");
tcint->SaveAs("nll.gif");
tcint->SaveAs("nll.eps");

}


Double_t getName(TFile& resultFile, bool NR_180, bool NR_Sign, bool RhoM_180, 
		 bool RhoM_Sign, bool Rho0_180, bool Rho0_Sign) {

  TString name = "PhaseScan";
  if(NR_180 == true)  name += "_NR_180"; 
  if(NR_Sign == true)  name += "_NR_Sign"; 
  if(RhoM_180 == true)  name += "_Rho-_180"; 
  if(RhoM_Sign == true)  name += "_Rho-_Sign"; 
  if(Rho0_180 == true)  name += "_Rho0_180"; 
  if(Rho0_Sign == true)  name += "_Rho0_Sign"; 
  RooFitResult* Res = (RooFitResult*) resultFile.Get(name);
  Double_t nll;
  if(Res!=0) nll = Res->minNll();
  return(nll);
}
