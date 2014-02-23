{

BdkTriBin::setCalcAreaNSteps(200);

gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
BdkBinning b(5);

b.divideBin(BdkP2(0.5, 0.25));

b.buildWgtLists();
b.setBinNamesInDalitzFrac();

const BdkTriBin & bin = b.bin(BdkP2(0.5, 0.25));
const BdkBinWgtList & wgts = bin.wgtList(BdkIsoOp::P32);

TCanvas * cc = b.Draw(kTRUE);

wgts.i().bin()->Draw(kBlack, 3);
wgts.r().bin()->Draw(kRed, 3);
wgts.rr().bin()->Draw(kYellow, 3);
wgts.e().bin()->Draw(kGreen, 3);
wgts.er().bin()->Draw(kBlue, 3);
wgts.err().bin()->Draw(kMagenta, 3);

cout << "i=" << wgts.i().weight() << endl;
cout << "r=" << wgts.r().weight() << endl;
cout << "rr=" << wgts.rr().weight() << endl;
cout << "e=" << wgts.e().weight() << endl;
cout << "er=" << wgts.er().weight() << endl;
cout << "err=" << wgts.err().weight() << endl;




cc->SaveAs("testWgt.eps");
}
