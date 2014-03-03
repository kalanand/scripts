{

BdkTriBin::setCalcAreaNSteps(200);

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
BdkBinning b(5);

b.setBinNamesIndex();

TCanvas * cc = b.Draw(kTRUE, kTRUE);

b.bin(0).Draw(kBlack, 3);
b.bin(0).R(kTRUE).Draw(kRed, 3);
b.bin(0).RR(kTRUE).Draw(kYellow, 3);
b.bin(0).E(kTRUE).Draw(kGreen, 3);
b.bin(0).ER(kTRUE).Draw(kBlue, 3);
b.bin(0).ERR(kTRUE).Draw(kMagenta, 3);

b.bin(0).writeName();
b.bin(0).R(kTRUE).writeName();
b.bin(0).RR(kTRUE).writeName();
b.bin(0).E(kTRUE).writeName();
b.bin(0).ER(kTRUE).writeName();
b.bin(0).ERR(kTRUE).writeName();


cc->SaveAs("testSym0.eps");
}
