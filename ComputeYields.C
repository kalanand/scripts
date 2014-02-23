{


   // ****** Define your selection: "WP95" or "WP80" *********** //
   //char* selection = "WP95";
   char* selection = "WP80";

   const double MINMASS = 60.0;
   const double MAXMASS = 120.0;
   const double ETMIN = 20.0;

   // ********* Input and output file names ******************* //
   char* inputFileName = "../allTPtrees_35invpb.root";
   char* inputTreeName = "WP95ToHLT/fitter_tree";
   if(selection == "WP80") 
      inputTreeName = "WP95ToHLT/fitter_tree";


   // ********* List energy scale here ******************* //
//    const double ebScale = 1.0151;
//    const double eeScale = 1.0330;

   const double ebScale = 1.0045;
   const double eeScale = 1.045;


   // ********* List all cuts here ******************* //
   // WP95 cuts
   const float MAX_MissingHits_WP95      = 1.0;
   const float MIN_Dist_WP95             = 0.0;
   const float MIN_Dcot_WP95             = 0.0;

   const float cut_EB_trackRel03_WP95    = 0.15;
   const float cut_EB_ecalRel03_WP95     = 2.00;
   const float cut_EB_hcalRel03_WP95     = 0.12;
   const float cut_EB_sigmaIetaIeta_WP95 = 0.01;
   const float cut_EB_deltaPhi_WP95      = 0.8;
   const float cut_EB_deltaEta_WP95      = 0.007;
   const float cut_EB_HoverE_WP95        = 0.15;

   const float cut_EE_trackRel03_WP95    = 0.08;
   const float cut_EE_ecalRel03_WP95     = 0.06;
   const float cut_EE_hcalRel03_WP95     = 0.05;
   const float cut_EE_sigmaIetaIeta_WP95 = 0.03;
   const float cut_EE_deltaPhi_WP95      = 0.7;
   const float cut_EE_deltaEta_WP95      = 0.01; // = 1000000.0; no cut in 3.6.X
   const float cut_EE_HoverE_WP95        = 0.07;



   // ********* List all cuts here ******************* //
   // WP80 cuts
   const float MAX_MissingHits_WP80      = 0.0;
   const float MIN_Dist_WP80             = 0.02;
   const float MIN_Dcot_WP80             = 0.02;

   const float cut_EB_trackRel03_WP80    = 0.09;
   const float cut_EB_ecalRel03_WP80     = 0.07;
   const float cut_EB_hcalRel03_WP80     = 0.10;
   const float cut_EB_sigmaIetaIeta_WP80 = 0.01;
   const float cut_EB_deltaPhi_WP80      = 0.06;
   const float cut_EB_deltaEta_WP80      = 0.004;
   const float cut_EB_HoverE_WP80        = 0.04;

   const float cut_EE_trackRel03_WP80    = 0.04;
   const float cut_EE_ecalRel03_WP80     = 0.05;
   const float cut_EE_hcalRel03_WP80     = 0.025;
   const float cut_EE_sigmaIetaIeta_WP80 = 0.03;
   const float cut_EE_deltaPhi_WP80      = 0.03;
   const float cut_EE_deltaEta_WP80      = 0.007; // = 1000000.0;  no cut in 3.6.X
   const float cut_EE_HoverE_WP80        = 0.025;



   // ********* List all cuts here ******************* //
   float MAX_MissingHits      = MAX_MissingHits_WP95;
   float MIN_Dist             = MIN_Dist_WP95;
   float MIN_Dcot             = MIN_Dcot_WP95;

   float cut_EB_trackRel03    = cut_EB_trackRel03_WP95;
   float cut_EB_ecalRel03     = cut_EB_ecalRel03_WP95;
   float cut_EB_hcalRel03     = cut_EB_hcalRel03_WP95;
   float cut_EB_sigmaIetaIeta = cut_EB_sigmaIetaIeta_WP95;
   float cut_EB_deltaPhi      = cut_EB_deltaPhi_WP95;
   float cut_EB_deltaEta      = cut_EB_deltaEta_WP95;
   float cut_EB_HoverE        = cut_EB_HoverE_WP95;

   float cut_EE_trackRel03    = cut_EE_trackRel03_WP95;
   float cut_EE_ecalRel03     = cut_EE_ecalRel03_WP95;
   float cut_EE_hcalRel03     = cut_EE_hcalRel03_WP95;
   float cut_EE_sigmaIetaIeta = cut_EE_sigmaIetaIeta_WP95;
   float cut_EE_deltaPhi      = cut_EE_deltaPhi_WP95;
   float cut_EE_deltaEta      = cut_EE_deltaEta_WP95; 
   float cut_EE_HoverE        = cut_EE_HoverE_WP95;


   if(selection == "WP80") {
      MAX_MissingHits      = MAX_MissingHits_WP80;
      MIN_Dist             = MIN_Dist_WP80;
      MIN_Dcot             = MIN_Dcot_WP80;

      cut_EB_trackRel03    = cut_EB_trackRel03_WP80;
      cut_EB_ecalRel03     = cut_EB_ecalRel03_WP80;
      cut_EB_hcalRel03     = cut_EB_hcalRel03_WP80;
      cut_EB_sigmaIetaIeta = cut_EB_sigmaIetaIeta_WP80;
      cut_EB_deltaPhi      = cut_EB_deltaPhi_WP80;
      cut_EB_deltaEta      = cut_EB_deltaEta_WP80;
      cut_EB_HoverE        = cut_EB_HoverE_WP80;

      cut_EE_trackRel03    = cut_EE_trackRel03_WP80;
      cut_EE_ecalRel03     = cut_EE_ecalRel03_WP80;
      cut_EE_hcalRel03     = cut_EE_hcalRel03_WP80;
      cut_EE_sigmaIetaIeta = cut_EE_sigmaIetaIeta_WP80;
      cut_EE_deltaPhi      = cut_EE_deltaPhi_WP80;
      cut_EE_deltaEta      = cut_EE_deltaEta_WP80; 
      cut_EE_HoverE        = cut_EE_HoverE_WP80;
   }


   TFile fin(inputFileName, "read");
   TTree* in_tree = (TTree*) fin.Get( inputTreeName );

  Long64_t nentries = in_tree->GetEntries();
  Float_t mass, tag_sc_eta, probe_sc_eta, tag_sc_phi,probe_sc_et,probe_sc_phi, 
     tag_gsfEle_dcot, tag_gsfEle_dist,
     probe_gsfEle_dcot, probe_gsfEle_dist,
     tag_gsfEle_deltaPhi, tag_gsfEle_deltaEta, tag_gsfEle_HoverE,
     tag_gsfEle_trackiso_dr03, tag_gsfEle_et, tag_gsfEle_ecaliso_dr03, 
     tag_gsfEle_hcaliso_dr03, tag_sc_et, probe_sc_et,
     tag_gsfEle_ecalDrivenSeed, probe_gsfEle_ecalDrivenSeed,
     tag_gsfEle_missingHits, probe_gsfEle_missingHits, 
     tag_gsfEle_isEB, tag_gsfEle_isEE, probe_gsfEle_isEB,probe_gsfEle_isEE,
     probe_gsfEle_deltaPhi, probe_gsfEle_deltaEta, probe_gsfEle_HoverE,
     probe_gsfEle_trackiso_dr03, probe_gsfEle_et, probe_gsfEle_ecaliso_dr03, 
     probe_gsfEle_hcaliso_dr03, probe_gsfEle_sigmaIetaIeta, tag_gsfEle_sigmaIetaIeta; 



  Int_t probe_passingId80;
  UInt_t run,lumi,event;

  char outputFileName[50];
  sprintf(outputFileName, "ZeeEvents_corrected_35invpb_%s.root", selection);
  TFile newfile( outputFileName, "RECREATE");
  TTree *newtree = in_tree->CloneTree(0);

  in_tree->SetBranchAddress("mass", &mass);
  in_tree->SetBranchAddress("tag_sc_eta", &tag_sc_eta);
  in_tree->SetBranchAddress("probe_sc_eta", &probe_sc_eta);
  in_tree->SetBranchAddress("tag_gsfEle_ecalDrivenSeed", &tag_gsfEle_ecalDrivenSeed);
  in_tree->SetBranchAddress("probe_gsfEle_ecalDrivenSeed", &probe_gsfEle_ecalDrivenSeed);
 //  in_tree->SetBranchAddress("tag_dcot", &tag_gsfEle_dcot);
//   in_tree->SetBranchAddress("tag_dist", &tag_gsfEle_dist);
  in_tree->SetBranchAddress("probe_dcot", &probe_gsfEle_dcot);
  in_tree->SetBranchAddress("probe_dist", &probe_gsfEle_dist);
  in_tree->SetBranchAddress("tag_gsfEle_missingHits", &tag_gsfEle_missingHits);
  in_tree->SetBranchAddress("probe_gsfEle_missingHits", &probe_gsfEle_missingHits);
  in_tree->SetBranchAddress("tag_gsfEle_deltaPhi", &tag_gsfEle_deltaPhi);
  in_tree->SetBranchAddress("tag_gsfEle_deltaEta", &tag_gsfEle_deltaEta);
  in_tree->SetBranchAddress("tag_gsfEle_HoverE", &tag_gsfEle_HoverE);
  in_tree->SetBranchAddress("tag_gsfEle_trackiso", &tag_gsfEle_trackiso_dr03);
  in_tree->SetBranchAddress("tag_gsfEle_et", &tag_gsfEle_et);
  in_tree->SetBranchAddress("tag_gsfEle_ecaliso", &tag_gsfEle_ecaliso_dr03);
  in_tree->SetBranchAddress("tag_gsfEle_hcaliso", &tag_gsfEle_hcaliso_dr03);
  in_tree->SetBranchAddress("tag_sc_et", &tag_sc_et);
  in_tree->SetBranchAddress("probe_sc_et", &probe_sc_et);
  in_tree->SetBranchAddress("tag_gsfEle_isEB", &tag_gsfEle_isEB);
  in_tree->SetBranchAddress("tag_gsfEle_isEE", &tag_gsfEle_isEE);
  in_tree->SetBranchAddress("probe_gsfEle_isEB", &probe_gsfEle_isEB);
  in_tree->SetBranchAddress("probe_gsfEle_isEE", &probe_gsfEle_isEE);
  in_tree->SetBranchAddress("probe_gsfEle_deltaPhi", &probe_gsfEle_deltaPhi);
  in_tree->SetBranchAddress("probe_gsfEle_deltaEta", &probe_gsfEle_deltaEta);
  in_tree->SetBranchAddress("probe_gsfEle_HoverE", &probe_gsfEle_HoverE);
  in_tree->SetBranchAddress("probe_gsfEle_et", &probe_gsfEle_et);
  in_tree->SetBranchAddress("probe_gsfEle_trackiso", &probe_gsfEle_trackiso_dr03);
  in_tree->SetBranchAddress("probe_gsfEle_ecaliso", &probe_gsfEle_ecaliso_dr03);
  in_tree->SetBranchAddress("probe_gsfEle_hcaliso", &probe_gsfEle_hcaliso_dr03);
  in_tree->SetBranchAddress("probe_gsfEle_sigmaIetaIeta", &probe_gsfEle_sigmaIetaIeta);
  in_tree->SetBranchAddress("tag_gsfEle_sigmaIetaIeta", &tag_gsfEle_sigmaIetaIeta);

  in_tree->SetBranchAddress("run", &run);
  in_tree->SetBranchAddress("lumi", &lumi);
  in_tree->SetBranchAddress("event", &event);
  in_tree->SetBranchAddress("tag_sc_phi", &tag_sc_phi);
  in_tree->SetBranchAddress("probe_sc_et", &probe_sc_et);
  in_tree->SetBranchAddress("probe_sc_phi", &probe_sc_phi);


  int count = 0;
  int EBEBcount = 0;
  int EBEEcount = 0;
  int EEEEcount = 0;

  int lastEvent = 0;
  int lastRun = 0;

  for (Int_t i = 0; i < nentries; i++) {
    in_tree->GetEntry(i);

    if( abs(tag_gsfEle_ecalDrivenSeed)==0 ) continue;
    if( abs(probe_gsfEle_ecalDrivenSeed)==0 ) continue;
    if( tag_gsfEle_missingHits>MAX_MissingHits )  continue;
    if( probe_gsfEle_missingHits>MAX_MissingHits )  continue;    
//     if( fabs(tag_gsfEle_dcot)<MIN_Dcot && 
//     fabs(tag_gsfEle_dist)<MIN_Dist ) continue;
    if( fabs(probe_gsfEle_dcot)<MIN_Dcot && 
    fabs(probe_gsfEle_dist)<MIN_Dist ) continue;


    bool select = true;
    double trackiso=0.0;
    double ecaliso=0.0;
    double hcaliso=0.0;

// all tag cuts
    if(tag_gsfEle_isEB) {
       tag_sc_et = ebScale*tag_sc_et;
       if( tag_sc_et<ETMIN ) select = false; 
       tag_gsfEle_et = ebScale*tag_gsfEle_et;

       if( fabs(tag_gsfEle_deltaPhi)>cut_EB_deltaPhi ) select = false;
       if( fabs(tag_gsfEle_deltaEta)>cut_EB_deltaEta ) select = false;
       if( tag_gsfEle_HoverE>cut_EB_HoverE ) select = false;
       if( tag_gsfEle_trackiso_dr03/tag_gsfEle_et>cut_EB_trackRel03 ) select = false;
       if( tag_gsfEle_ecaliso_dr03/tag_gsfEle_et>cut_EB_ecalRel03 ) select = false;
       if( tag_gsfEle_hcaliso_dr03/tag_gsfEle_et>cut_EB_hcalRel03 ) select = false;
       if( tag_gsfEle_sigmaIetaIeta > cut_EB_sigmaIetaIeta ) select = false;
    }
    else if(tag_gsfEle_isEE) {
       tag_sc_et = eeScale*tag_sc_et;
       if( tag_sc_et<20.0 ) select = false; 
       tag_gsfEle_et = eeScale*tag_gsfEle_et;

       if( fabs(tag_gsfEle_deltaPhi)>cut_EE_deltaPhi ) select = false;
       if( fabs(tag_gsfEle_deltaEta)>cut_EE_deltaEta ) select = false;
       if( tag_gsfEle_HoverE>cut_EE_HoverE ) select = false;
       if( tag_gsfEle_trackiso_dr03/tag_gsfEle_et>cut_EE_trackRel03 ) select = false;
       if( tag_gsfEle_ecaliso_dr03/tag_gsfEle_et>cut_EE_ecalRel03 ) select = false;
       if( tag_gsfEle_hcaliso_dr03/tag_gsfEle_et>cut_EE_hcalRel03 ) select = false;
       if( tag_gsfEle_sigmaIetaIeta > cut_EE_sigmaIetaIeta) select = false;
    }
    else select = false;


// all probe cuts
    if(probe_gsfEle_isEB) {
       probe_sc_et = ebScale*probe_sc_et;
       if( probe_sc_et<ETMIN ) select = false; 
       probe_gsfEle_et = ebScale*probe_gsfEle_et;

       if( fabs(probe_gsfEle_deltaPhi)>cut_EB_deltaPhi ) select = false;
       if( fabs(probe_gsfEle_deltaEta)>cut_EB_deltaEta ) select = false;
       if( probe_gsfEle_HoverE>cut_EB_HoverE ) select = false;
       if( probe_gsfEle_trackiso_dr03/probe_gsfEle_et>cut_EB_trackRel03 ) select = false;
       if( probe_gsfEle_ecaliso_dr03/probe_gsfEle_et>cut_EB_ecalRel03 ) select = false;
       if( probe_gsfEle_hcaliso_dr03/probe_gsfEle_et>cut_EB_hcalRel03 ) select = false;
       if( probe_gsfEle_sigmaIetaIeta > cut_EB_sigmaIetaIeta ) select = false;
    }
    else if(probe_gsfEle_isEE) {
       probe_sc_et = eeScale*probe_sc_et;
       if( probe_sc_et<20.0 ) select = false; 
       probe_gsfEle_et = eeScale*probe_gsfEle_et;

       if( fabs(probe_gsfEle_deltaPhi)>cut_EE_deltaPhi ) select = false;
       if( fabs(probe_gsfEle_deltaEta)>cut_EE_deltaEta ) select = false;
       if( probe_gsfEle_HoverE>cut_EE_HoverE ) select = false;
       if( probe_gsfEle_trackiso_dr03/probe_gsfEle_et>cut_EE_trackRel03 ) select = false;
       if( probe_gsfEle_ecaliso_dr03/probe_gsfEle_et>cut_EE_ecalRel03 ) select = false;
       if( probe_gsfEle_hcaliso_dr03/probe_gsfEle_et>cut_EE_hcalRel03 ) select = false;
       if( probe_gsfEle_sigmaIetaIeta > cut_EE_sigmaIetaIeta) select = false;
    }
    else select = false;

// mass cut 
    if(tag_gsfEle_isEB && probe_gsfEle_isEB)  mass = ebScale*mass; 
    if(tag_gsfEle_isEE && probe_gsfEle_isEE)  mass = eeScale*mass; 
    if( (tag_gsfEle_isEB && probe_gsfEle_isEE) || 
    (tag_gsfEle_isEE && probe_gsfEle_isEB) ) mass = sqrt(ebScale*eeScale)*mass;    
    if( mass<MINMASS || mass>MAXMASS ) select = false; 


// now we have selected events, just do not double count the TT events
    if(select) {
       if( !(lastEvent==event && lastRun==run) )
       printf("%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", 
       run,lumi,event,mass,tag_sc_et,tag_sc_eta,tag_sc_phi,
       probe_sc_et,probe_sc_eta,probe_sc_phi);
        
       if( !(lastEvent==event && lastRun==run) ) count ++;
       if( !(lastEvent==event && lastRun==run) ) newtree->Fill();


       if(tag_gsfEle_isEB && probe_gsfEle_isEB && !(lastEvent==event && lastRun==run))
       EBEBcount ++;
       if(tag_gsfEle_isEE && probe_gsfEle_isEE && !(lastEvent==event && lastRun==run))
       EEEEcount ++;
       if( ((tag_gsfEle_isEB && probe_gsfEle_isEE) || 
       (tag_gsfEle_isEE && probe_gsfEle_isEB)) && !(lastEvent==event && lastRun==run))
       EBEEcount ++;

       lastEvent=event;
       lastRun = run;
    }
  }

  cout << "Z signal yield = " << count << endl;
  cout << "Z signal yield EBEB = " << EBEBcount << endl;
  cout << "Z signal yield EBEE = " << EBEEcount << endl;
  cout << "Z signal yield EEEE = " << EEEEcount << endl;
  cout << "Z signal yield total = " << EBEBcount+EBEEcount+EEEEcount << endl;

  newfile.cd();
  newtree->Write("tree");
  newfile .Close();

}
