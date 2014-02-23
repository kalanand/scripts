// Makes plots that look at the pion PID efficiencies

void piPidStudy(){
  
  TCut cutNoPid = cutSigReg + TCut("Pippidbit>-1");  
  RooDataSet * mcQQ = (RooDataSet *)dataSetFromTree(qqTree, cutNoPid);
  RooDataSet * mcBB = (RooDataSet *)dataSetFromTree(bbTree, cutNoPid);
  RooDataSet * mcSig = (RooDataSet *)dataSetFromTree(sigTree, cutNoPid);
  
  TCut theCut = cutNoPid;
  theCut += cutPiT;
  RooDataSet * mcQQpid = (RooDataSet *)dataSetFromTree(qqTree, theCut);
  RooDataSet * mcBBpid = (RooDataSet *)dataSetFromTree(bbTree, theCut);
  RooDataSet * mcSigpid = (RooDataSet *)dataSetFromTree(sigTree, theCut);
  
  RooPlot * qqPid = Pippidbit->frame();
  qqPid->SetTitle("pi Pid for QQ");
  mcQQ->plotOn(qqPid);
  
  RooPlot * bbPid = Pippidbit->frame();
  bbPid->SetTitle("pi Pid for BB");
  mcBB->plotOn(bbPid);
  
  RooPlot * sigPid = Pippidbit->frame();
  sigPid->SetTitle("pi Pid for SIG");
  mcSig->plotOn(sigPid);
  
  TCanvas * can = new TCanvas("can", "can", 1200, 400);
  can->Divide(3,1);
  can->cd(1);
  bbPid->Draw();
  can->cd(2);
  qqPid->Draw();
  can->cd(3);
  sigPid->Draw();
  
  cout << "eff sig = " << (double)mcSigpid->numEntries() / mcSig->numEntries()
       << "eff qq  = " << (double)mcQQpid->numEntries() / mcQQ->numEntries()
       << "eff bb  = " << (double)mcBBpid->numEntries() / mcBB->numEntries()
       << endl;
}

