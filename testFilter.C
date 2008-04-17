
TString dirName = "/uscms_data/d1/kalanand/trash/";
TString prefix  = "Summer09_7TeV_Summer09-7TeV-ZeeJets-Pt_";

void testFilter() {
  Process( 0,15 );
  Process( 15,20 );
  Process( 20,30 );
  Process( 30,50 );
  Process( 50,80 );
  Process(80, 120);
  Process( 120,170 );
  Process( 170,230, "/uscms_data/d1/kalanand/ZJet_new_1.root");
  Process( 230,300, "/uscms_data/d1/kalanand/ZJet_new_1.root");

}




void Process(int ptLow, int ptHigh, 
	     char* outFileName = "/uscms_data/d1/kalanand/ZJet_new.root" ){

  TString bin = Form("%d_%d", ptLow, ptHigh);
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$ pt bin : " << bin << endl;
  TString fName = dirName + prefix + bin + TString(".root");
  TFile f(fName);
  TTree* inputTree = (TTree*) f.Get("ZJet");
  inputTree->SetBranchStatus("*",1);

  Float_t mZee;
  Float_t ePlusPt;
  Float_t eMinusPt;
  Float_t ePlusEta;
  Float_t eMinusEta;
  Float_t ePlus_trackiso;	   
  Float_t eMinus_trackiso;
  Float_t ePlus_ecaliso;
  Float_t eMinus_ecaliso;
  Bool_t iseMinusLoose;
  Bool_t isePlusLoose;

  inputTree->SetBranchAddress("mZee",            &mZee);
  inputTree->SetBranchAddress("ePlusPt",         &ePlusPt);
  inputTree->SetBranchAddress("eMinusPt",        &eMinusPt);
  inputTree->SetBranchAddress("ePlusEta",        &ePlusEta);
  inputTree->SetBranchAddress("eMinusEta",       &eMinusEta);
  inputTree->SetBranchAddress("eMinus_ecaliso",  &eMinus_ecaliso);
  inputTree->SetBranchAddress("ePlus_ecaliso",   &ePlus_ecaliso);
  inputTree->SetBranchAddress("ePlus_trackiso",  &ePlus_trackiso);	   
  inputTree->SetBranchAddress("eMinus_trackiso", &eMinus_trackiso);
  inputTree->SetBranchAddress("iseMinusLoose",   &iseMinusLoose);
  inputTree->SetBranchAddress("isePlusLoose",    &isePlusLoose);


  inputTree->SetBranchStatus("mZeeGen",   0);
  inputTree->SetBranchStatus("Z_VxGen",   0);
  inputTree->SetBranchStatus("Z_VyGen",   0);
  inputTree->SetBranchStatus("Z_VzGen",   0);
  inputTree->SetBranchStatus("ePlusChargeGen",  0);
  inputTree->SetBranchStatus("ePlusVxGen",      0);
  inputTree->SetBranchStatus("ePlusVyGen",      0);
  inputTree->SetBranchStatus("ePlusVzGen",      0);
  inputTree->SetBranchStatus("ePlusClassificationGen", 0);      
  inputTree->SetBranchStatus("eMinusChargeGen", 0);
  inputTree->SetBranchStatus("eMinusVxGen",     0);
  inputTree->SetBranchStatus("eMinusVyGen",     0);
  inputTree->SetBranchStatus("eMinusVzGen",     0);
  inputTree->SetBranchStatus("eMinusClassificationGen", 0);
  inputTree->SetBranchStatus("NumGenJetAlgo", 0);
  inputTree->SetBranchStatus("NumGenJets",    0);   
  inputTree->SetBranchStatus("JetGenType",    0);
  inputTree->SetBranchStatus("JetGenInvisibleEnergy",  0);
  inputTree->SetBranchStatus("JetGenAuxiliaryEnergy",  0);


  inputTree->SetBranchStatus("NumJPTJetAlgo", 0);
  inputTree->SetBranchStatus("NumJPTJets",    0);   
  inputTree->SetBranchStatus("JetJPTType",    0);

  inputTree->SetBranchStatus("JetJPTMaxEInEmTowers",  0);
  inputTree->SetBranchStatus("JetJPTMaxEInHadTowers", 0);
  inputTree->SetBranchStatus("JetJPTHadEnergyInHB",  0);
  inputTree->SetBranchStatus("JetJPTHadEnergyInHO",  0);
  inputTree->SetBranchStatus("JetJPTHadEnergyInHE",  0);
  inputTree->SetBranchStatus("JetJPTHadEnergyInHF",  0);
  inputTree->SetBranchStatus("JetJPTEmEnergyInEB",   0);
  inputTree->SetBranchStatus("JetJPTEmEnergyInEE",   0);
  inputTree->SetBranchStatus("JetJPTEmEnergyInHF",   0);
  inputTree->SetBranchStatus("JetJPTTowersArea",     0);
  inputTree->SetBranchStatus("JetJPTN90",            0);
  inputTree->SetBranchStatus("JetJPTN60",            0);

  inputTree->SetBranchStatus("NumCorJetAlgo", 0);
  inputTree->SetBranchStatus("NumCorJets",    0);   
  inputTree->SetBranchStatus("JetCorType",    0);
  inputTree->SetBranchStatus("JetCorMaxEInEmTowers",  0);
  inputTree->SetBranchStatus("JetCorMaxEInHadTowers", 0);
  inputTree->SetBranchStatus("JetCorHadEnergyInHB",  0);
  inputTree->SetBranchStatus("JetCorHadEnergyInHO",  0);
  inputTree->SetBranchStatus("JetCorHadEnergyInHE",  0);
  inputTree->SetBranchStatus("JetCorHadEnergyInHF",  0);
  inputTree->SetBranchStatus("JetCorEmEnergyInEB",   0);
  inputTree->SetBranchStatus("JetCorEmEnergyInEE",   0);
  inputTree->SetBranchStatus("JetCorEmEnergyInHF",   0);
  inputTree->SetBranchStatus("JetCorTowersArea",     0);
  inputTree->SetBranchStatus("JetCorN90",            0);
  inputTree->SetBranchStatus("JetCorN60",            0);


  inputTree->SetBranchStatus("NumRecoJetAlgo", 0);
  inputTree->SetBranchStatus("NumRecoJets",    0);   
  inputTree->SetBranchStatus("JetRecoType",    0);
  inputTree->SetBranchStatus("JetRecoMaxEInEmTowers",  0);
  inputTree->SetBranchStatus("JetRecoMaxEInHadTowers", 0);
  inputTree->SetBranchStatus("JetRecoHadEnergyInHB",  0);
  inputTree->SetBranchStatus("JetRecoHadEnergyInHO",  0);
  inputTree->SetBranchStatus("JetRecoHadEnergyInHE",  0);
  inputTree->SetBranchStatus("JetRecoHadEnergyInHF",  0);
  inputTree->SetBranchStatus("JetRecoEmEnergyInEB",   0);
  inputTree->SetBranchStatus("JetRecoEmEnergyInEE",   0);
  inputTree->SetBranchStatus("JetRecoEmEnergyInHF",   0);
  inputTree->SetBranchStatus("JetRecoTowersArea",     0);
  inputTree->SetBranchStatus("JetRecoN90",            0);
  inputTree->SetBranchStatus("JetRecoN60",            0);

  inputTree->SetBranchStatus("run",  0);
  inputTree->SetBranchStatus("event",  0);
  inputTree->SetBranchStatus("eventsize", 0); 
  inputTree->SetBranchStatus("nZee",        0);   
  inputTree->SetBranchStatus("Z_Vx",        0);
  inputTree->SetBranchStatus("Z_Vy",        0);
  inputTree->SetBranchStatus("Z_Vz",        0);
  inputTree->SetBranchStatus("ePlusVx",     0);
  inputTree->SetBranchStatus("ePlusVy",     0);
  inputTree->SetBranchStatus("ePlusVz",      0);
  inputTree->SetBranchStatus("ePlusIsolation", 0);
  inputTree->SetBranchStatus("ePlus_sc_y", 0);
  inputTree->SetBranchStatus("ePlus_sc_z", 0);
  inputTree->SetBranchStatus("ePlus_sc_R", 0);
  inputTree->SetBranchStatus("ePlus_sc_Rt", 0);
  inputTree->SetBranchStatus("ePlus_EoverPout", 0);
  inputTree->SetBranchStatus("ePlus_EoverPin",  0);
  inputTree->SetBranchStatus("ePlus_InvEMinusInvP", 0);
  inputTree->SetBranchStatus("ePlus_BremFraction",  0);
  inputTree->SetBranchStatus("ePlus_DeltaEtaIn", 0);
  inputTree->SetBranchStatus("ePlus_DeltaPhiIn", 0);
  inputTree->SetBranchStatus("ePlus_DeltaPhiOut", 0);
  inputTree->SetBranchStatus("ePlus_DeltaEtaOut", 0);
  inputTree->SetBranchStatus("ePlus_Trackmom_calo", 0);
  inputTree->SetBranchStatus("ePlus_Trackmom_vtx", 0);	  
  inputTree->SetBranchStatus("ePlus_vx", 0);
  inputTree->SetBranchStatus("ePlus_vy", 0);
  inputTree->SetBranchStatus("ePlus_vz", 0);	  
  inputTree->SetBranchStatus("ePlus_td0", 0);
  inputTree->SetBranchStatus("ePlus_d0", 0);	 	  
  inputTree->SetBranchStatus("eMinusVx",     0);
  inputTree->SetBranchStatus("eMinusVy",     0);
  inputTree->SetBranchStatus("eMinusVz",      0);
  inputTree->SetBranchStatus("eMinusIsolation", 0);
  inputTree->SetBranchStatus("eMinus_sc_y", 0);
  inputTree->SetBranchStatus("eMinus_sc_z", 0);
  inputTree->SetBranchStatus("eMinus_sc_R", 0);
  inputTree->SetBranchStatus("eMinus_sc_Rt", 0);
  inputTree->SetBranchStatus("eMinus_EoverPout", 0);
  inputTree->SetBranchStatus("eMinus_EoverPin",  0);
  inputTree->SetBranchStatus("eMinus_InvEMinusInvP", 0);
  inputTree->SetBranchStatus("eMinus_BremFraction",  0);
  inputTree->SetBranchStatus("eMinus_DeltaEtaIn", 0);
  inputTree->SetBranchStatus("eMinus_DeltaPhiIn", 0);
  inputTree->SetBranchStatus("eMinus_DeltaPhiOut", 0);
  inputTree->SetBranchStatus("eMinus_DeltaEtaOut", 0);
  inputTree->SetBranchStatus("eMinus_Trackmom_calo", 0);
  inputTree->SetBranchStatus("eMinus_Trackmom_vtx", 0);	  
  inputTree->SetBranchStatus("eMinus_vx", 0);
  inputTree->SetBranchStatus("eMinus_vy", 0);
  inputTree->SetBranchStatus("eMinus_vz", 0);	  
  inputTree->SetBranchStatus("eMinus_td0", 0);
  inputTree->SetBranchStatus("eMinus_d0", 0);	 	  


  TFile f1(outFileName,"update");
  TTree* outputTree = inputTree->CloneTree(0); 
  outputTree->SetName(bin);


  float weight = 0.0, crosssection = 0.0;
  if( bin.CompareTo("0_15")==0 )     { weight = 2.06877;   crosssection  =  6430.0; }
  if( bin.CompareTo("15_20")==0 )    { weight = 0.0693239; crosssection  =  230.0 ; }
  if( bin.CompareTo("20_30")==0 )    { weight = 0.0750612; crosssection  =  211.0 ; }
  if( bin.CompareTo("30_50")==0 )    { weight = 0.0441168; crosssection  =  142.0 ; }
  if( bin.CompareTo("50_80")==0 )    { weight = 0.0273266; crosssection  =  56.8 ; }
  if( bin.CompareTo("80_120")==0 )   { weight = 0.00694582; crosssection  =  18.8 ; }
  if( bin.CompareTo("120_170")==0 )  { weight = 0.00215549; crosssection  =  5.4 ; }
  if( bin.CompareTo("170_230")==0 )  { weight = 0.000521085; crosssection  =  1.55 ; }
  if( bin.CompareTo("230_300")==0 ) { weight = 0.000175759; crosssection  =  0.45  ; }

  TVector2 v2( weight, crosssection);
  TString vname = TString("weight_") + bin;
  f1.cd();
  v2.Write(vname);


  // Make candidate selection and dump to new Ntuple

  for (Int_t entry(0); entry<inputTree->GetEntries(); entry++) {
    inputTree->GetEntry(entry); 
    if(entry%10000==0) std::cout<<"****************** Event # "<< entry <<std::endl;
    if( BoolCutResult(mZee, ePlusPt, ePlusEta, eMinusPt, 
		      eMinusEta, eMinus_trackiso, ePlus_trackiso, 
		      eMinus_ecaliso, ePlus_ecaliso) &&
	isePlusLoose && iseMinusLoose ) outputTree->Fill();
  }  



  f1.cd();
  outputTree->Write();
  f1.Close();
}





bool BoolCutResult( float mZ, float e1Pt, float e1Eta, 
		    float e2Pt, float e2Eta, float e1trackiso, 
		    float e2trackiso, 
		    float e1ecaliso, float e2ecaliso) {

bool result = true;

// Z mass cut
if( fabs(mZ-91.2) > 10 ) result = false;


// electron pT cut
if( e1Pt < 20.0 ) result = false;
if( e2Pt < 20.0 ) result = false;


// electron acceptance
if( !((fabs(e1Eta)<1.4442) || 
(fabs(e1Eta)>1.560 && fabs(e1Eta)<2.5)) ) result = false;
if( !((fabs(e2Eta)<1.4442) || 
(fabs(e2Eta)>1.560 && fabs(e2Eta)<2.5)) ) result = false;

// electron isolation: e1 barrel
 if( fabs(e1Eta)<1.4442) { 
   if(e1trackiso > 7.2 )  result = false; 
   if(e1ecaliso > 5.7 )  result = false; 
   //   if(e1hcaliso > 8.1 )  result = false; 
 }

 // e1 endcap
 if( fabs(e1Eta)>1.560 && fabs(e1Eta)<2.5 ) { 
   if(e1trackiso > 5.1 )  result = false; 
   if(e1ecaliso > 5.0 )  result = false; 
   // if(e1hcaliso > 3.4 )  result = false; 
 }

 // e2 barrel
 if( fabs(e2Eta)<1.4442) { 
   if(e2trackiso > 7.2 )  result = false; 
   if(e2ecaliso > 5.7 )  result = false;
   // if(e2hcaliso > 8.1 )  result = false; 
 }

 // e2 endcap
 if( fabs(e2Eta)>1.560 && fabs(e2Eta)<2.5 ) { 
   if(e2trackiso > 5.1 )  result = false; 
   if(e2ecaliso > 5.0 )  result = false; 
   //   if(e2hcaliso > 3.4 )  result = false; 
 }

return result;
}
