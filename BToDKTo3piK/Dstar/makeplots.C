#include <iomanip>
#include "TMatrix.h"

using namespace std;

void makeplots() {
  //nominal fit corresponds to i==3

  for(int i=3; i<4; i++) { Loop(i); }
}



void Loop(int GROUP) {
  gROOT->Reset();
  gStyle->Reset();
  gROOT->SetStyle("BABAR");

  bool doNorm = true;
  int numBINS = 30;


  cout<< "*************************** Now doing GROUP # " << GROUP << endl;


  RooRealVar fit_S23("fit_S23","m^{2}(#pi^{+}#pi^{0})",0,3); //This is pi+ pi0  
  RooRealVar fit_S31("fit_S31","m^{2}(#pi^{-}#pi^{0})",0,3); //This is pi- pi0  
  RooRealVar fit_S12("fit_S12","m^{2}(#pi^{-}#pi^{+})",0,3); //This is pi+ pi-
  fit_S23.setBins(numBINS);
  fit_S31.setBins(numBINS);
  fit_S12.setBins(numBINS);


  //read in the background dataset
  TFile hello("hhPi0.root","READ");
  TTree* suptree = (TTree*)hello.Get("pipipi0_data");
//   TTree* supBg = (TTree*)hello.Get("pipipi0_ccbar");
  gROOT->cd();


  TCut remKs("(sqrt(fit_S12)<0.489||sqrt(fit_S12)>0.508)&& abs(hh_SignedDistance)<2.0");
//   TTree* sigtree = suptree->CopyTree(remKs && "abs(Dmass-1.8637)<0.016");

  TTree* sigtree = suptree->CopyTree(remKs && "abs(Dmass-1.8637)<0.016");
  TTree* bkgtree = suptree->CopyTree(remKs && "Dmass>1.93 && Dmass<1.99");
  // if taking background shape from MC instead of data sideband: 
  // TTree* bkgtree = supBg->CopyTree(remKs && "abs(Dmass-1.864)<0.016 && Flag>2");


  RooArgList myList;
  myList.add(RooArgSet(fit_S23,fit_S31,fit_S12));
  RooDataSet *data = new RooDataSet("data","data",sigtree,myList);
  RooDataSet *d2 =  new RooDataSet("d2","d2",bkgtree,myList); 


  BdkDalitzCfg* dalitzCfg = new BdkDalitzCfg("dalitzCfg","dalitzCfg");
  dalitzCfg->getParameters(RooArgSet())->readFromFile("../BToDKTo3piK/Dstar/dalitzCfg.par");
  dalitzCfg->getParameters(RooArgSet())->Print("v");


  // define signal pdf
  BdkPdf2DpolyDalitz eff("eff","",fit_S23, fit_S31);
  eff.parameters().readFromFile("eff.par");
  //eff.parameters().readFromFile("eff_noPID.par");
  BdkPdfDDalitz dstar("dstar","dstar",fit_S23,fit_S31,BdkDalitzBase::D0);
  dstar.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());
  dstar.parameters().readFromFile("PPP_NominalFit.par");
  // BdkDDalitzAmp::setUseFixedWidth(kTRUE); // for fixed width

  // define bkg pdf
  TH2D bgHist("bgHist","2d bkgd hist from data sideband",6,0,3,6,0,3);
  bkgtree->Draw("fit_S23:fit_S31>>bgHist","","goff");
  ((BdkDalitzBase*) dstar.getPdf())->weightBins(&bgHist,true,0.001);
  BdkDalitzHist bkgpdf("bkgpdf", "bkgpdf", BdkDalitzBase::D0, 
		       BdkDalitzBase::PPP0, fit_S23, fit_S31, bgHist);



  // define signal+bkg pdf
  RooRealVar fsig("fsig","fsig",0.982);
  RooRealVar fbkg("fbkg","fbkg",1.0-fsig.getVal());
  RooAddPdf mypdf("mypdf","mypdf",RooArgList(*dstar.getPdf(),bkgpdf),
		  RooArgList(fsig,fbkg));

  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-6);
  cfg->setEpsRel(1E-6);
  mypdf.setIntegratorConfig(*cfg);


  RooArgSet allPars = dstar.parameters();
  RooRealVar* var1 = allPars.find("dstar.pdf.dalitzAmp.Nonres_amp");
  RooRealVar* var2 = allPars.find("dstar.pdf.dalitzAmp.Nonres_phase");
  RooRealVar* var3 = allPars.find("dstar.pdf.dalitzAmp.Rho+_amp");
  RooRealVar* var4 = allPars.find("dstar.pdf.dalitzAmp.Rho+_phase");
  RooRealVar* var5 = allPars.find("dstar.pdf.dalitzAmp.Rho-_amp");
  RooRealVar* var6 = allPars.find("dstar.pdf.dalitzAmp.Rho-_phase");
  RooRealVar* var7 = allPars.find("dstar.pdf.dalitzAmp.Rho0_amp");
  RooRealVar* var8 = allPars.find("dstar.pdf.dalitzAmp.Rho0_phase");

  RooRealVar* var9 = allPars.find("dstar.pdf.dalitzAmp.Rho2s+_amp");
  RooRealVar* var10 = allPars.find("dstar.pdf.dalitzAmp.Rho2s+_phase");
  RooRealVar* var11 = allPars.find("dstar.pdf.dalitzAmp.Rho2s-_amp");
  RooRealVar* var12 = allPars.find("dstar.pdf.dalitzAmp.Rho2s-_phase");
  RooRealVar* var13 = allPars.find("dstar.pdf.dalitzAmp.Rho2s0_amp");
  RooRealVar* var14 = allPars.find("dstar.pdf.dalitzAmp.Rho2s0_phase");
  RooRealVar* var15 = allPars.find("dstar.pdf.dalitzAmp.Sigma_amp");
  RooRealVar* var16 = allPars.find("dstar.pdf.dalitzAmp.Sigma_phase");
  RooRealVar* var17 = allPars.find("dstar.pdf.dalitzAmp.Omega_amp");
  RooRealVar* var18 = allPars.find("dstar.pdf.dalitzAmp.Omega_phase");

  RooRealVar* var19 = allPars.find("dstar.pdf.dalitzAmp.F0_amp");
  RooRealVar* var20 = allPars.find("dstar.pdf.dalitzAmp.F0_phase");
  RooRealVar* var21 = allPars.find("dstar.pdf.dalitzAmp.F2_amp");
  RooRealVar* var22 = allPars.find("dstar.pdf.dalitzAmp.F2_phase");
  RooRealVar* var23 = allPars.find("dstar.pdf.dalitzAmp.F0_1370_amp");
  RooRealVar* var24 = allPars.find("dstar.pdf.dalitzAmp.F0_1370_phase");
  RooRealVar* var25 = allPars.find("dstar.pdf.dalitzAmp.Rho1700+_amp");
  RooRealVar* var26 = allPars.find("dstar.pdf.dalitzAmp.Rho1700+_phase");
  RooRealVar* var27 = allPars.find("dstar.pdf.dalitzAmp.Rho1700-_amp");
  RooRealVar* var28 = allPars.find("dstar.pdf.dalitzAmp.Rho1700-_phase");
  RooRealVar* var29 = allPars.find("dstar.pdf.dalitzAmp.Rho17000_amp");
  RooRealVar* var30 = allPars.find("dstar.pdf.dalitzAmp.Rho17000_phase");

  RooRealVar* var31 = allPars.find("dstar.pdf.dalitzAmp.F0_1500_amp");
  RooRealVar* var32 = allPars.find("dstar.pdf.dalitzAmp.F0_1500_phase");
  RooRealVar* var33 = allPars.find("dstar.pdf.dalitzAmp.F0_1710_amp");
  RooRealVar* var34 = allPars.find("dstar.pdf.dalitzAmp.F0_1710_phase");
  RooRealVar* var35 = allPars.find("dstar.pdf.dalitzAmp.F2P1525_amp");
  RooRealVar* var36 = allPars.find("dstar.pdf.dalitzAmp.F2P1525_phase");

  RooRealVar* var37 = allPars.find("dstar.pdf.dalitzAmp.NRPW_PM_amp");
  RooRealVar* var38 = allPars.find("dstar.pdf.dalitzAmp.NRPW_PM_phase");
  RooRealVar* var39 = allPars.find("dstar.pdf.dalitzAmp.NRPW_M0_amp");
  RooRealVar* var40 = allPars.find("dstar.pdf.dalitzAmp.NRPW_M0_phase");
  RooRealVar* var41 = allPars.find("dstar.pdf.dalitzAmp.NRPW_0P_amp");
  RooRealVar* var42 = allPars.find("dstar.pdf.dalitzAmp.NRPW_0P_phase");

  RooRealVar* var43 = allPars.find("dstar.pdf.dalitzAmp.resRadius");


  RooRealVar* var51 = allPars.find("dstar.pdf.dalitzAmp.Rho+_mass");
  RooRealVar* var52 = allPars.find("dstar.pdf.dalitzAmp.Rho+_width");
  RooRealVar* var53 = allPars.find("dstar.pdf.dalitzAmp.Rho2s+_mass");
  RooRealVar* var54 = allPars.find("dstar.pdf.dalitzAmp.Rho2s+_width");
  RooRealVar* var55 = allPars.find("dstar.pdf.dalitzAmp.Rho1700+_mass");
  RooRealVar* var56 = allPars.find("dstar.pdf.dalitzAmp.Rho1700+_width");
  RooRealVar* var57 = allPars.find("dstar.pdf.dalitzAmp.Sigma_mass");
  RooRealVar* var58 = allPars.find("dstar.pdf.dalitzAmp.Sigma_width");
  RooRealVar* var59 = allPars.find("dstar.pdf.dalitzAmp.Omega_mass");
  RooRealVar* var60 = allPars.find("dstar.pdf.dalitzAmp.Omega_width");
  RooRealVar* var61 = allPars.find("dstar.pdf.dalitzAmp.F0_mass");
  RooRealVar* var62 = allPars.find("dstar.pdf.dalitzAmp.F0_width");
  RooRealVar* var63 = allPars.find("dstar.pdf.dalitzAmp.F0_1370_mass");
  RooRealVar* var64 = allPars.find("dstar.pdf.dalitzAmp.F0_1370_width");
  RooRealVar* var65 = allPars.find("dstar.pdf.dalitzAmp.F0_1500_mass");
  RooRealVar* var66 = allPars.find("dstar.pdf.dalitzAmp.F0_1500_width");
  RooRealVar* var67 = allPars.find("dstar.pdf.dalitzAmp.F0_1710_mass");
  RooRealVar* var68 = allPars.find("dstar.pdf.dalitzAmp.F0_1710_width");
  RooRealVar* var69 = allPars.find("dstar.pdf.dalitzAmp.F2_mass");
  RooRealVar* var70 = allPars.find("dstar.pdf.dalitzAmp.F2_width");
  RooRealVar* var71 = allPars.find("dstar.pdf.dalitzAmp.F2P1525_mass");
  RooRealVar* var72 = allPars.find("dstar.pdf.dalitzAmp.F2P1525_width");


  var1->setConstant(false);
  var2->setConstant(false);
  var5->setConstant(false);
  var6->setConstant(false);
  var7->setConstant(false);
  var8->setConstant(false);




  // Removing  nothing
  if(GROUP==0) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setConstant(false);
    var16->setConstant(false);
    var17->setConstant(false);
    var18->setConstant(false);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setConstant(false);
    var22->setConstant(false);
    var23->setConstant(false);
    var24->setConstant(false);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setConstant(false);
    var32->setConstant(false);
    var33->setConstant(false);
    var34->setConstant(false);
    var35->setConstant(false);
    var36->setConstant(false);
    var37->setConstant(false);
    var38->setConstant(false);
    var39->setConstant(false);
    var40->setConstant(false);
    var41->setConstant(false);
    var42->setConstant(false);
  }

  // Removing Omega and f2'(1525)
  if(GROUP==1) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setConstant(false);
    var16->setConstant(false);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setConstant(false);
    var22->setConstant(false);
    var23->setConstant(false);
    var24->setConstant(false);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setConstant(false);
    var32->setConstant(false);
    var33->setConstant(false);
    var34->setConstant(false);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setConstant(false);
    var38->setConstant(false);
    var39->setConstant(false);
    var40->setConstant(false);
    var41->setConstant(false);
    var42->setConstant(false);
  }


  // Removing above and sigma
  if(GROUP==2) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setVal(0.0);
    var16->setVal(0.0);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setConstant(false);
    var22->setConstant(false);
    var23->setConstant(false);
    var24->setConstant(false);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setConstant(false);
    var32->setConstant(false);
    var33->setConstant(false);
    var34->setConstant(false);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setConstant(false);
    var38->setConstant(false);
    var39->setConstant(false);
    var40->setConstant(false);
    var41->setConstant(false);
    var42->setConstant(false);
  }

  // ************ This is our nominal fit
  // Removing above and p-wave NR but not sigma 
  if(GROUP==3) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setConstant(false);
    var16->setConstant(false);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setConstant(false);
    var22->setConstant(false);
    var23->setConstant(false);
    var24->setConstant(false);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setConstant(false);
    var32->setConstant(false);
    var33->setConstant(false);
    var34->setConstant(false);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
  }




  // Removing above and f0(1370), f0(1500),
  // f0(1710) and f2
  if(GROUP==4) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setVal(0.0);
    var16->setVal(0.0);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setVal(0.0);
    var22->setVal(0.0);
    var23->setVal(0.0);
    var24->setVal(0.0);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setVal(0.0);
    var32->setVal(0.0);
    var33->setVal(0.0);
    var34->setVal(0.0);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
  }


  // Removing above and Rho(1700)
  if(GROUP==5) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setVal(0.0);
    var16->setVal(0.0);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setVal(0.0);
    var22->setVal(0.0);
    var23->setVal(0.0);
    var24->setVal(0.0);
    var25->setVal(0.0);
    var26->setVal(0.0);
    var27->setVal(0.0);
    var28->setVal(0.0);
    var29->setVal(0.0);
    var30->setVal(0.0);
    var31->setVal(0.0);
    var32->setVal(0.0);
    var33->setVal(0.0);
    var34->setVal(0.0);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
  }


  // Removing above and Rho(1450)
  if(GROUP==6) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setVal(0.0);
    var10->setVal(0.0);
    var11->setVal(0.0);
    var12->setVal(0.0);
    var13->setVal(0.0);
    var14->setVal(0.0);
    var15->setVal(0.0);
    var16->setVal(0.0);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setVal(0.0);
    var22->setVal(0.0);
    var23->setVal(0.0);
    var24->setVal(0.0);
    var25->setVal(0.0);
    var26->setVal(0.0);
    var27->setVal(0.0);
    var28->setVal(0.0);
    var29->setVal(0.0);
    var30->setVal(0.0);
    var31->setVal(0.0);
    var32->setVal(0.0);
    var33->setVal(0.0);
    var34->setVal(0.0);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
  }


  // Removing above and f0
  if(GROUP==7) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setVal(0.0);
    var10->setVal(0.0);
    var11->setVal(0.0);
    var12->setVal(0.0);
    var13->setVal(0.0);
    var14->setVal(0.0);
    var15->setVal(0.0);
    var16->setVal(0.0);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setVal(0.0);
    var20->setVal(0.0);
    var21->setVal(0.0);
    var22->setVal(0.0);
    var23->setVal(0.0);
    var24->setVal(0.0);
    var25->setVal(0.0);
    var26->setVal(0.0);
    var27->setVal(0.0);
    var28->setVal(0.0);
    var29->setVal(0.0);
    var30->setVal(0.0);
    var31->setVal(0.0);
    var32->setVal(0.0);
    var33->setVal(0.0);
    var34->setVal(0.0);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
  }



  // nominal fit (same as GROUP 3) but Changing the resRadius to 0.0
  if(GROUP==8) {
    var1->setConstant(false);
    var2->setConstant(false);
    var5->setConstant(false);
    var6->setConstant(false);
    var7->setConstant(false);
    var8->setConstant(false);
    var9->setConstant(false);
    var10->setConstant(false);
    var11->setConstant(false);
    var12->setConstant(false);
    var13->setConstant(false);
    var14->setConstant(false);
    var15->setConstant(false);
    var16->setConstant(false);
    var17->setVal(0.0);
    var18->setVal(0.0);
    var19->setConstant(false);
    var20->setConstant(false);
    var21->setConstant(false);
    var22->setConstant(false);
    var23->setConstant(false);
    var24->setConstant(false);
    var25->setConstant(false);
    var26->setConstant(false);
    var27->setConstant(false);
    var28->setConstant(false);
    var29->setConstant(false);
    var30->setConstant(false);
    var31->setConstant(false);
    var32->setConstant(false);
    var33->setConstant(false);
    var34->setConstant(false);
    var35->setVal(0.0);
    var36->setVal(0.0);
    var37->setVal(0.0);
    var38->setVal(0.0);
    var39->setVal(0.0);
    var40->setVal(0.0);
    var41->setVal(0.0);
    var42->setVal(0.0);
    var43->setVal(0.0);
  }

   if(doNorm==true) {               
     dstar.parameters(); 
     BdkDDalitzAmp::normalizeAll(); 
   }


  RooFitResult* res = mypdf.fitTo(*data,"mr");
  RooArgSet fitfractions = (RooArgSet&) dstar.fitFractions();
  fitfractions.Print("v");

//   printFitResult(res);


  ofstream of;
  of.open("cov.txt", ios_base::out);
  printCovMatrix(res,of);
  of.close();
                       

  ofstream of2;
  of2.open("corr.txt", ios_base::out);
  printCorMatrix(res,of2,kFALSE,"%.4f");
  of2.close();
             


//   dstar.parameters().writeToFile(TString("dstar_PPP_Group")+GROUP+TString(".par")); 



//   // Machinary to plot the third variable
//   Double_t md0 = 1.8645;
//   Double_t mpi0=0.1349766;
//   Double_t mpi = 0.13957;
//   Double_t total=md0*md0+mpi0*mpi0+mpi*mpi+mpi*mpi;
//   RooRealVar totalm("totalm","totalm",total);
//   RooFormulaVar fit_S31a("fit_S31a","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S12));
//   BdkPdfDDalitz dstar12("dstar12","dstar12",fit_S23,fit_S31a,BdkDalitzBase::D0); 
//   dstar12.setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());
//   dstar12.getPdf()->setIntegratorConfig(*cfg);
//   RooArgSet allPars12 = dstar12.parameters();
//   TIterator * iter12 = allPars12.createIterator();
//   RooRealVar * var;
//   while(0!= (var= (RooRealVar*)iter12->Next())) {
//     TString vname = TString(var->GetName()).ReplaceAll("dstar12","dstar");
//     if(vname.Contains("dstar.")) {
//       RooRealVar* arvar = allPars.find(vname);
//       var->setVal(arvar->getVal());
//     }
//   }
  
//   delete iter12;



//   RooAddPdf mypdf12("mypdf12","mypdf12",RooArgList(*dstar12.getPdf(),bkgpdf),
// 		    RooArgList(fsig,fbkg));




//   //Make the plots
//   //******************************************************
//   RooPlot* xframe = fit_S23.frame(0,3);
//   data->plotOn(xframe,MarkerSize(0.1),DrawOption("z"));
//   mypdf.plotOn(xframe);
//   xframe->getAttLine()->SetLineWidth(1);
//   xframe->getAttLine()->SetLineStyle(1);
//   xframe->SetTitleSize(0.08, "Y" );
//   xframe->SetTitleOffset(1.4, "Y");
//   xframe->GetYaxis()->SetTitle("");

//   RooPlot* yframe = fit_S31.frame(0,3);
//   data->plotOn(yframe,MarkerSize(0.1),DrawOption("z"));
//   mypdf.plotOn(yframe); 
//   yframe->getAttLine()->SetLineWidth(1);
//   yframe->getAttLine()->SetLineStyle(1);
//   yframe->SetTitleSize(0.08, "Y" );
//   yframe->SetTitleOffset(1.4, "Y");
//   yframe->GetYaxis()->SetTitle("");

//   RooPlot* zframe = fit_S12.frame(0,3);
//   data->plotOn(zframe,MarkerSize(0.1),DrawOption("z"));
//   mypdf12.plotOn(zframe);
//   zframe->getAttLine()->SetLineWidth(1);
//   zframe->getAttLine()->SetLineStyle(1); 
//   zframe->SetTitleSize(0.08, "Y" );
//   zframe->SetTitleOffset(1.4, "Y");
//   zframe->GetYaxis()->SetTitle("");




//   TCanvas *canvas = new TCanvas("canvas","allevents",930,400);
//   canvas->Divide(3,1);
//   canvas->cd(1);xframe->Draw();
//   canvas->cd(2);yframe->Draw();
//   canvas->cd(3);
//   zframe->Draw();
//   canvas->SaveAs(TString("ppp_Group")+GROUP+TString(".gif"));
//   canvas->SaveAs(TString("ppp_Group")+GROUP+TString(".eps"));
//   canvas->SaveAs(TString("ppp_Group")+GROUP+TString(".root"));
  //  delete canvas;


  //
  // Plot chi^2
  //

//   TH2D h("h","m^{2}(#pi^{-}#pi^{0}) : m^{2}(#pi^{+}#pi^{0})",
// 	 numBINS,0,3,numBINS,0,3);
//   sigtree->Draw("fit_S23:fit_S31>>h","","goff");
//   RooDataSet *fitdata = mypdf.generate(RooArgSet(fit_S23,fit_S31),data->numEntries());
//   RooFormulaVar fit_S12a("fit_S12","@0-@1-@2",RooArgSet(totalm,fit_S23,fit_S31));
//   fitdata->addColumn(fit_S12a);

//   TH2F* hFit = fitdata->createHistogram(fit_S31,fit_S23);
//   h.Sumw2();
//   hFit->Sumw2();
//   Double_t hintegral = h.Integral();
//   Double_t hFitintegral = hFit->Integral();
//   Double_t Normalization = hintegral/hFitintegral;
//   hFit->Scale(Normalization);

//   TH2D *hchi2 = new TH2D("hchi2","",numBINS,0,3,numBINS,0,3);  
//   TH1D* pullhist = new TH1D("pullhist","",25,-4,4); 
//   double chi2Total =0.0;
//   int ndof = numBINS*numBINS;
//   ndof--;   // histos are normalized to each other

//   Double_t chi2Total = 0;    // total chi2
//   for (int x=1; x <= numBINS; x++) {
//     for (int y=1; y <= numBINS; y++) {
      
//       Double_t bin1 =  h.GetBinContent(x,y);
//       Double_t bin2 = hFit->GetBinContent(x,y); 

//       if (bin2==0 && bin1==0) ndof--;    // no data -> one less dof
//       else {
// 	Double_t sqrtChi2 = (bin1-bin2);
// 	sqrtChi2 /= sqrt(h.GetBinError(x,y)**2 + hFit->GetBinError(x,y)**2);
// 	hchi2->SetBinContent(x,y,sqrtChi2);
// 	chi2Total += sqrtChi2**2;
// 	pullhist->Fill(sqrtChi2);
//       }
//     }
//   }


//   Double_t chi2Prob = TMath::Prob(chi2Total, ndof);
//   cout << "Chi2          = " << chi2Total << endl;
//   cout << "ndof          = " << ndof << endl;
//   cout << "Chi2/ndof     = " << chi2Total/ndof << endl;
//   cout << "Chi2 prob     = " << chi2Prob << endl;


//   TCanvas* aCanvas = new TCanvas("c3", "c3", 900, 600);
//   gStyle->SetOptStat(0);
//   aCanvas->Divide(2,2);
//   aCanvas->cd(1);
//   h.Draw("colz");
//   aCanvas->cd(2);
//   hFit->Draw("colz");
//   aCanvas->cd(3);
//   hchi2->SetMaximum(4);
//   hchi2->SetMinimum(-4);
//   hchi2->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
//   hchi2->Draw("colz");
//   aCanvas->cd(4);  
//   pullhist->SetTitle("#sqrt{#chi^{2}}: m^{2}(#pi^{-}#pi^{0}), m^{2}(#pi^{+}#pi^{0})");
//   pullhist->Fit("gaus");
//   gStyle->SetStatX(0.99);
//   gStyle->SetStatY(0.9);
//   gStyle->SetOptFit(11);
//   pullhist->Draw();
//   aCanvas->Draw();
//   aCanvas->SaveAs(TString("ppp_chi_Group")+GROUP+TString(".gif"));
//   aCanvas->SaveAs(TString("ppp_chi_Group")+GROUP+TString(".eps"));
//   aCanvas->SaveAs(TString("ppp_chi_Group")+GROUP+TString(".root"));
}  //end the macro



/* ******************Printing Results *************************/

void printCovMatrix(RooFitResult * res, ostream & stream);
void getCovCorMatrix(RooFitResult* res, TMatrix& cov, TMatrix& cor);

void printCorMatrix(const TMatrix& cor, ostream& stream = cout,
                    Bool_t printLatex = kFALSE, const char* fmt = "%.3f");

void printCorMatrix(RooFitResult* res, ostream& stream = cout,
                    Bool_t printLatex = kFALSE, const char* fmt = "%.3f");

//--------------------------------------------------------------------
void printChanges(RooFitResult * res, ostream & os = cout) 
{
  for (int i = 0 ; i < res->floatParsFinal().getSize(); ++i) {
  
    RooRealVar *final = (RooRealVar*)res->floatParsFinal().at(i);
    
    os << "Change in " << setw(20) << final->GetName();
    const double init  = ((RooRealVar*)res->floatParsInit().at(i))->getVal();
    double diff = final->getVal() - init;
    
    double error = -1;
    if (final->hasError())
      error = ((RooRealVar*)res->floatParsFinal().at(i))->getError();
    else if (final->hasAsymError()) {
      if (diff<0) error = fabs(final->getAsymErrorLo());
      else error = fabs(final->getAsymErrorHi());
    }
    if (error > 0) diff /= error;
    else diff = -999;

    os << " = " << diff << endl;
  }
}

//--------------------------------------------------------------------
void printFitResult(RooFitResult * res) {
  if (0 == res) {
    return;
  }

  TMatrix cov;
  TMatrix cor;

  res->Print();
  cout << endl << "Fit status = "<<res->status()<<endl<<endl;

  getCovCorMatrix(res, cov, cor);
  cout << "covariance matrix" << endl;
  cov.Print();
  cout << "correlation matrix" << endl;
  cor.Print();

  printChanges(res);
}


//--------------------------------------------------------------------
void printFitResult(RooFitResult * res, ostream & stream) {
  if (0 == res) {
    return;
  }

  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);

  printCovMatrix(res, stream);
  printCorMatrix(cor, stream);

  res->printToStream(stream, RooPrintable::Verbose);
  stream << endl << "Fit status = "<<res->status()<<endl<<endl;
  printChanges(res, stream);

}

//--------------------------------------------------------------------
// Print the correlation matrix
void printCorMatrix(RooFitResult* res, ostream& stream, 
                    Bool_t printLatex, const char* fmt)
{
  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);
  printCorMatrix(cor, stream, printLatex, fmt);
}


void printCorMatrix(const TMatrix& cor, ostream& stream, 
                    Bool_t printLatex, const char* fmt)
{
  for (int i = 0; i < cor.GetNrows(); ++i) {
    for (int j = 0; j < cor.GetNcols(); ++j) {
      TString s;
      if (cor(i,j)==1) s = "1";
      else s.Form(fmt,(double)cor(i,j));

      if (printLatex) {
        if (j>0) stream << " & ";
        stream << s;      
      }
      else stream << setw(13) << s;
    }
    if (printLatex) stream << "\\\\";
    stream << endl;
  }
}

//--------------------------------------------------------------------
// Print the covariance matrix
void printCovMatrix(RooFitResult * res, ostream & stream) {

  TMatrix cov, cor;  
  getCovCorMatrix(res, cov, cor);

  stream << "//$" << endl;
  for (int i = 0; i < cov.GetNrows(); ++i) {
    stream << "// ";
    for (int j = 0; j < cov.GetNcols(); ++j) {
      stream << setprecision(8) << setw(15) << cov(i,j);
    }
    stream << endl;
  }
  stream << "//$" << endl << endl;
}

//--------------------------------------------------------------------
// Fill the correlation and covariance matrix from the fit result
void getCovCorMatrix(RooFitResult* res, TMatrix& cov, TMatrix& cor)
{
  if (!res) return;
  const int size = res->floatParsFinal().getSize();
  cov.ResizeTo(size, size);
  cor.ResizeTo(size, size);
  for (int i = 0; i < size; ++i) {
    RooRealVar * vari = (RooRealVar *)(res->floatParsFinal().at(i));
    for (int j = 0; j < size; ++j) {
      RooRealVar * varj = (RooRealVar *)(res->floatParsFinal().at(j));
      cor(i,j) = res->correlation(*vari, *varj);
      cov(i,j) = cor(i,j) * vari->getError() * varj->getError();
    }
  }
}




