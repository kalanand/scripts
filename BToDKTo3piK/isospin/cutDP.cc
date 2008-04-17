// Apply the Dalitz plot cut on Kalanand's trees

void cutDP() {

  TTree * sigTree2 = new TTree("signal", "signal");
  TTree * bgdTree2 = new TTree("background", "background");
  
  cutCopy(sigTree2);
  cutCopy(bgdTree2);

  TFile file2("data-iso-inDP.root", "RECREATE");
  sigTree2->Write();
  bgdTree2->Write();

  cout << "New ntuples written to " << file2.GetName() << endl;
  file2.Close();
}


void cutCopy(TTree * tree2) {
  const char * name = tree2->GetName();

  // We will put the tree info in these:
  double s12, s13;
 
  // Open the old tree and prepare to extract it into s12, s13:
  TFile file("data-kalanand.root");
  const TTree  * tree1 = (TTree*)file.Get(name);
  Long64_t nIn = tree1->GetEntries();
  
  TBranch * br12 = tree1->GetBranch("s12");
  TBranch * br13 = tree1->GetBranch("s13");

  br12->SetAddress(&s12);
  br13->SetAddress(&s13);

  // Point the branches of the new tree into the same doubles:
  tree2->Branch("s12", &s12, "s12/D");
  tree2->Branch("s13", &s13, "s13/D");

  // Apply the cut and copy:  
  for (int e = 0; e < tree1->GetEntries(); ++e) {
    br12->GetEntry(e);
    br13->GetEntry(e);
    
    BdkP2 point(s12, s13);
    if (point.inDalitz()) {
      tree2->Fill();
    }
  }

  Long64_t nOut = tree2->GetEntries();

  file.Close();

  cout << nOut << " " << name << " points passed out of " 
       << nIn << endl;
}

