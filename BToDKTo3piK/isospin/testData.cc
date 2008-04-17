{

TFile file("data-iso-inDP.root");
const TTree  * sig = (TTree*)file.Get("signal");
const TTree  * bgd = (TTree*)file.Get("background");


BdkTriBin::setCalcAreaNSteps(200);

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
BdkBinning b(5, sig, bgd);

b.calcWFs();


b.setBinNamesNEventsSR();
TCanvas * ccNSR = b.Draw(kTRUE);
ccNSR->SaveAs("testDataNSR.eps");

b.setBinNamesNEventsSB();
TCanvas * ccNSB = b.Draw(kTRUE, kTRUE, "SB");
ccNSB->SaveAs("testDataNSB.eps");

b.setBinNamesWaveFunc();
TCanvas * ccWF = b.Draw(kTRUE);
ccWF->SaveAs("testDataWF.eps");

b.setBinNamesWaveFuncErr();
TCanvas * ccWFE = b.Draw(kTRUE);
ccWFE->SaveAs("testDataWFE.eps");

TCanvas * cc = b.Draw();
cc->SaveAs("testDataSR.eps");

TCanvas * cc = b.Draw(kFALSE, kTRUE, "SB");
cc->SaveAs("testDataSB.eps");

int sumSig = 0;
for (int bin = 0; bin < b.nBins(); ++bin) {
  sumSig += b.bin(bin).nEventsSR();
}

cout << "Read " << sig->GetEntries() << " signal events" << endl
<< "BdkBinning has " << sumSig << endl;


}

