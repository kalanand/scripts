{

BdkTriBin::setCalcAreaNSteps(200);

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
BdkBinning b(5);

b.divideBin(BdkP2(0.5, 0.25));

b.setBinNamesIndex();
//b.setBinNamesInDalitzFrac();

TCanvas * cc = b.Draw(kTRUE, kTRUE);


cc->SaveAs("testDiv.eps");
}
