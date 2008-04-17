//setup the neural networks output pdf holder

#include "../BToDKTo3piK/globals/pdfs.hh"
#include "../BToDKTo3piK/globals/vars.hh"

void setupPdfHolder(BdkPdfDKDalitz::COORD coord = BdkPdfDKDalitz::POLAR,
                    RooAbsReal* var12 = m12,
                    RooAbsReal* var13 = m13) {

  cout<<"----- START setupPdfHolder() -----"<<endl;

  // Inititalize efficiencies
  eff.init("eff", "GoodD efficiency", *var12, *var13);
  effOther.init("effOther", "BadD efficiency", *var12, *var13);
  effDKBadD.init("effDKBadD", "DKBadD efficiency", *var12, *var13);

  const char* effPar = "../BToDKTo3piK/params/eff.par";

  cout << "Reading efficiency parameters from "<<effPar<<endl;
  eff.parameters().readFromFile(effPar);
  effOther.parameters().readFromFile(effPar);
  effDKBadD.parameters().readFromFile("../BToDKTo3piK/params/syst/effDKBadD.par");
  
  nnHolder.init("nnHolder", 
                "N1 PDF holder", 
                nnout,
                "../BToDKTo3piK/params/nnout.root");
  
  bknnHolder.init("bknnHolder",
                  "N2 PDF holder", 
                  bknnout, 
                  "../BToDKTo3piK/params/bknnout.root");
  
 
  cout << "Using "<< (coord==BdkPdfDKDalitz::POLAR ? "POLAR" : "CARTESIAN")
       << " coordinates for signal PDF."<<endl;

  Bool_t flipQqGoodDFlavor = false;  // for systematic studies

  dalitzHolderN.init("dalitzHolderN", 
                     "dalitz PDF Holder forB-",
                     *var12,
                     *var13,
                     BdkDalitzBase::D0,
                     coord,
                     (BdkDalitzEff*)eff.getPdf(), 
                     (BdkDalitzEff*)effOther.getPdf(),
                     "../BToDKTo3piK/params/hist_dpix.root",
                     "../BToDKTo3piK/params/hist_dkbadd.root",
                     *blindMode,
                     0,
                     flipQqGoodDFlavor);

  dalitzHolderP.init("dalitzHolderP", 
                     "dalitz PDF Holder for B+",
                     *var12,
                     *var13,
                     BdkDalitzBase::D0BAR,
                     coord,
                     (BdkDalitzEff*)eff.getPdf(),
                     (BdkDalitzEff*)effOther.getPdf(),
                     "../BToDKTo3piK/params/hist_dpix.root",
                     "../BToDKTo3piK/params/hist_dkbadd.root",
                     *blindMode,
                     &dalitzHolderN,   // this makes them share the same dalitzAmps
                     flipQqGoodDFlavor);  		
	
  cout<< "------- End setupPdfHolder() ----" <<endl; 
}
