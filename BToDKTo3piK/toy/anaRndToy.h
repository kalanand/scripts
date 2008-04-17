//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 11 06:47:43 2007 by ROOT version 4.04/02b
// from TTree fitParData/Fit Parameters DataSet
// found on file: rndToy-1.root
//////////////////////////////////////////////////////////

#ifndef anaRndToy_h
#define anaRndToy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class anaRndToy : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
   Double_t        dalitzHolderN_sigGoodD0_rho;
   Double_t        dalitzHolderN_sigGoodD0_rho_err;
   Double_t        dalitzHolderN_sigGoodD0_rho_aerr_lo;
   Double_t        dalitzHolderN_sigGoodD0_rho_aerr_hi;
   Double_t        dalitzHolderN_sigGoodD0_theta;
   Double_t        dalitzHolderN_sigGoodD0_theta_err;
   Double_t        dalitzHolderN_sigGoodD0_theta_aerr_lo;
   Double_t        dalitzHolderN_sigGoodD0_theta_aerr_hi;
   Double_t        dalitzHolderP_sigGoodD0_rho;
   Double_t        dalitzHolderP_sigGoodD0_rho_err;
   Double_t        dalitzHolderP_sigGoodD0_rho_aerr_lo;
   Double_t        dalitzHolderP_sigGoodD0_rho_aerr_hi;
   Double_t        dalitzHolderP_sigGoodD0_theta;
   Double_t        dalitzHolderP_sigGoodD0_theta_err;
   Double_t        dalitzHolderP_sigGoodD0_theta_aerr_lo;
   Double_t        dalitzHolderP_sigGoodD0_theta_aerr_hi;
   Double_t        pdfOnResDK_DPiXFrac;
   Double_t        pdfOnResDK_DPiXFrac_err;
   Double_t        pdfOnResDK_DPiXFrac_aerr_lo;
   Double_t        pdfOnResDK_DPiXFrac_aerr_hi;
   Double_t        pdfOnResDK_DpiGoodD0NumEvts;
   Double_t        pdfOnResDK_DpiGoodD0NumEvts_err;
   Double_t        pdfOnResDK_DpiGoodD0NumEvts_aerr_lo;
   Double_t        pdfOnResDK_DpiGoodD0NumEvts_aerr_hi;
   Double_t        pdfOnResDK_asymSigGoodD;
   Double_t        pdfOnResDK_asymSigGoodD_err;
   Double_t        pdfOnResDK_asymSigGoodD_aerr_lo;
   Double_t        pdfOnResDK_asymSigGoodD_aerr_hi;
   Double_t        pdfOnResDK_qqBadD0NumEvts;
   Double_t        pdfOnResDK_qqBadD0NumEvts_err;
   Double_t        pdfOnResDK_qqBadD0NumEvts_aerr_lo;
   Double_t        pdfOnResDK_qqBadD0NumEvts_aerr_hi;
   Double_t        pdfOnResDK_sigGoodD0NumEvts;
   Double_t        pdfOnResDK_sigGoodD0NumEvts_err;
   Double_t        pdfOnResDK_sigGoodD0NumEvts_aerr_lo;
   Double_t        pdfOnResDK_sigGoodD0NumEvts_aerr_hi;
   Double_t        pdfOnResDK_totBBNumEvts;
   Double_t        pdfOnResDK_totBBNumEvts_err;
   Double_t        pdfOnResDK_totBBNumEvts_aerr_lo;
   Double_t        pdfOnResDK_totBBNumEvts_aerr_hi;
   Double_t        NLL;
   Double_t        NLL_err;
   Double_t        NLL_aerr_lo;
   Double_t        NLL_aerr_hi;
   Double_t        xyNLL;
   Double_t        xyNLL_err;
   Double_t        xyNLL_aerr_lo;
   Double_t        xyNLL_aerr_hi;
   Double_t        rhoPgen;
   Double_t        rhoNgen;
   Double_t        thetaPgen;
   Double_t        thetaNgen;
   Double_t        cov_00;
   Double_t        cov_10;
   Double_t        cov_11;
   Double_t        cov_20;
   Double_t        cov_21;
   Double_t        cov_22;
   Double_t        cov_30;
   Double_t        cov_31;
   Double_t        cov_32;
   Double_t        cov_33;
   Double_t        dalitzHolderN_sigGoodD0_rhoerr;
   Double_t        dalitzHolderN_sigGoodD0_rhopull;
   Double_t        dalitzHolderN_sigGoodD0_thetaerr;
   Double_t        dalitzHolderN_sigGoodD0_thetapull;
   Double_t        dalitzHolderP_sigGoodD0_rhoerr;
   Double_t        dalitzHolderP_sigGoodD0_rhopull;
   Double_t        dalitzHolderP_sigGoodD0_thetaerr;
   Double_t        dalitzHolderP_sigGoodD0_thetapull;
   Double_t        pdfOnResDK_DPiXFracerr;
   Double_t        pdfOnResDK_DPiXFracpull;
   Double_t        pdfOnResDK_DpiGoodD0NumEvtserr;
   Double_t        pdfOnResDK_DpiGoodD0NumEvtspull;
   Double_t        pdfOnResDK_asymSigGoodDerr;
   Double_t        pdfOnResDK_asymSigGoodDpull;
   Double_t        pdfOnResDK_qqBadD0NumEvtserr;
   Double_t        pdfOnResDK_qqBadD0NumEvtspull;
   Double_t        pdfOnResDK_sigGoodD0NumEvtserr;
   Double_t        pdfOnResDK_sigGoodD0NumEvtspull;
   Double_t        pdfOnResDK_totBBNumEvtserr;
   Double_t        pdfOnResDK_totBBNumEvtspull;

   // List of branches
   TBranch        *b_dalitzHolderN_sigGoodD0_rho;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_rho_err;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_rho_aerr_lo;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_rho_aerr_hi;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_theta;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_theta_err;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_theta_aerr_lo;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_theta_aerr_hi;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rho;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rho_err;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rho_aerr_lo;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rho_aerr_hi;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_theta;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_theta_err;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_theta_aerr_lo;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_theta_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_DPiXFrac;   //!
   TBranch        *b_pdfOnResDK_DPiXFrac_err;   //!
   TBranch        *b_pdfOnResDK_DPiXFrac_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_DPiXFrac_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvts;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvts_err;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvts_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvts_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodD;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodD_err;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodD_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodD_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvts;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvts_err;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvts_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvts_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvts;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvts_err;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvts_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvts_aerr_hi;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvts;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvts_err;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvts_aerr_lo;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvts_aerr_hi;   //!
   TBranch        *b_NLL;   //!
   TBranch        *b_NLL_err;   //!
   TBranch        *b_NLL_aerr_lo;   //!
   TBranch        *b_NLL_aerr_hi;   //!
   TBranch        *b_xyNLL;   //!
   TBranch        *b_xyNLL_err;   //!
   TBranch        *b_xyNLL_aerr_lo;   //!
   TBranch        *b_xyNLL_aerr_hi;   //!
   TBranch        *b_rhoPgen;   //!
   TBranch        *b_rhoNgen;   //!
   TBranch        *b_thetaPgen;   //!
   TBranch        *b_thetaNgen;   //!
   TBranch        *b_cov_00;   //!
   TBranch        *b_cov_10;   //!
   TBranch        *b_cov_11;   //!
   TBranch        *b_cov_20;   //!
   TBranch        *b_cov_21;   //!
   TBranch        *b_cov_22;   //!
   TBranch        *b_cov_30;   //!
   TBranch        *b_cov_31;   //!
   TBranch        *b_cov_32;   //!
   TBranch        *b_cov_33;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_rhoerr;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_rhopull;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_thetaerr;   //!
   TBranch        *b_dalitzHolderN_sigGoodD0_thetapull;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rhoerr;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_rhopull;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_thetaerr;   //!
   TBranch        *b_dalitzHolderP_sigGoodD0_thetapull;   //!
   TBranch        *b_pdfOnResDK_DPiXFracerr;   //!
   TBranch        *b_pdfOnResDK_DPiXFracpull;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvtserr;   //!
   TBranch        *b_pdfOnResDK_DpiGoodD0NumEvtspull;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodDerr;   //!
   TBranch        *b_pdfOnResDK_asymSigGoodDpull;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvtserr;   //!
   TBranch        *b_pdfOnResDK_qqBadD0NumEvtspull;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvtserr;   //!
   TBranch        *b_pdfOnResDK_sigGoodD0NumEvtspull;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvtserr;   //!
   TBranch        *b_pdfOnResDK_totBBNumEvtspull;   //!

   anaRndToy(TTree *tree=0) { }
   virtual ~anaRndToy() { }
   virtual Int_t   Version() const {return 1;}
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) {fInput = input;}
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(anaRndToy,0);
};

#endif

#ifdef anaRndToy_cxx
void anaRndToy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rho",&dalitzHolderN_sigGoodD0_rho);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rho_err",&dalitzHolderN_sigGoodD0_rho_err);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rho_aerr_lo",&dalitzHolderN_sigGoodD0_rho_aerr_lo);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rho_aerr_hi",&dalitzHolderN_sigGoodD0_rho_aerr_hi);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.theta",&dalitzHolderN_sigGoodD0_theta);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.theta_err",&dalitzHolderN_sigGoodD0_theta_err);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.theta_aerr_lo",&dalitzHolderN_sigGoodD0_theta_aerr_lo);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.theta_aerr_hi",&dalitzHolderN_sigGoodD0_theta_aerr_hi);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rho",&dalitzHolderP_sigGoodD0_rho);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rho_err",&dalitzHolderP_sigGoodD0_rho_err);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rho_aerr_lo",&dalitzHolderP_sigGoodD0_rho_aerr_lo);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rho_aerr_hi",&dalitzHolderP_sigGoodD0_rho_aerr_hi);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.theta",&dalitzHolderP_sigGoodD0_theta);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.theta_err",&dalitzHolderP_sigGoodD0_theta_err);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.theta_aerr_lo",&dalitzHolderP_sigGoodD0_theta_aerr_lo);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.theta_aerr_hi",&dalitzHolderP_sigGoodD0_theta_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFrac",&pdfOnResDK_DPiXFrac);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFrac_err",&pdfOnResDK_DPiXFrac_err);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFrac_aerr_lo",&pdfOnResDK_DPiXFrac_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFrac_aerr_hi",&pdfOnResDK_DPiXFrac_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvts",&pdfOnResDK_DpiGoodD0NumEvts);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvts_err",&pdfOnResDK_DpiGoodD0NumEvts_err);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvts_aerr_lo",&pdfOnResDK_DpiGoodD0NumEvts_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvts_aerr_hi",&pdfOnResDK_DpiGoodD0NumEvts_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodD",&pdfOnResDK_asymSigGoodD);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodD_err",&pdfOnResDK_asymSigGoodD_err);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodD_aerr_lo",&pdfOnResDK_asymSigGoodD_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodD_aerr_hi",&pdfOnResDK_asymSigGoodD_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvts",&pdfOnResDK_qqBadD0NumEvts);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvts_err",&pdfOnResDK_qqBadD0NumEvts_err);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvts_aerr_lo",&pdfOnResDK_qqBadD0NumEvts_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvts_aerr_hi",&pdfOnResDK_qqBadD0NumEvts_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvts",&pdfOnResDK_sigGoodD0NumEvts);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvts_err",&pdfOnResDK_sigGoodD0NumEvts_err);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvts_aerr_lo",&pdfOnResDK_sigGoodD0NumEvts_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvts_aerr_hi",&pdfOnResDK_sigGoodD0NumEvts_aerr_hi);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvts",&pdfOnResDK_totBBNumEvts);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvts_err",&pdfOnResDK_totBBNumEvts_err);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvts_aerr_lo",&pdfOnResDK_totBBNumEvts_aerr_lo);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvts_aerr_hi",&pdfOnResDK_totBBNumEvts_aerr_hi);
   fChain->SetBranchAddress("NLL",&NLL);
   fChain->SetBranchAddress("NLL_err",&NLL_err);
   fChain->SetBranchAddress("NLL_aerr_lo",&NLL_aerr_lo);
   fChain->SetBranchAddress("NLL_aerr_hi",&NLL_aerr_hi);
   fChain->SetBranchAddress("xyNLL",&xyNLL);
   fChain->SetBranchAddress("xyNLL_err",&xyNLL_err);
   fChain->SetBranchAddress("xyNLL_aerr_lo",&xyNLL_aerr_lo);
   fChain->SetBranchAddress("xyNLL_aerr_hi",&xyNLL_aerr_hi);
   fChain->SetBranchAddress("rhoPgen",&rhoPgen);
   fChain->SetBranchAddress("rhoNgen",&rhoNgen);
   fChain->SetBranchAddress("thetaPgen",&thetaPgen);
   fChain->SetBranchAddress("thetaNgen",&thetaNgen);
   fChain->SetBranchAddress("cov_00",&cov_00);
   fChain->SetBranchAddress("cov_10",&cov_10);
   fChain->SetBranchAddress("cov_11",&cov_11);
   fChain->SetBranchAddress("cov_20",&cov_20);
   fChain->SetBranchAddress("cov_21",&cov_21);
   fChain->SetBranchAddress("cov_22",&cov_22);
   fChain->SetBranchAddress("cov_30",&cov_30);
   fChain->SetBranchAddress("cov_31",&cov_31);
   fChain->SetBranchAddress("cov_32",&cov_32);
   fChain->SetBranchAddress("cov_33",&cov_33);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rhoerr",&dalitzHolderN_sigGoodD0_rhoerr);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.rhopull",&dalitzHolderN_sigGoodD0_rhopull);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.thetaerr",&dalitzHolderN_sigGoodD0_thetaerr);
   fChain->SetBranchAddress("dalitzHolderN.sigGoodD0.thetapull",&dalitzHolderN_sigGoodD0_thetapull);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rhoerr",&dalitzHolderP_sigGoodD0_rhoerr);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.rhopull",&dalitzHolderP_sigGoodD0_rhopull);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.thetaerr",&dalitzHolderP_sigGoodD0_thetaerr);
   fChain->SetBranchAddress("dalitzHolderP.sigGoodD0.thetapull",&dalitzHolderP_sigGoodD0_thetapull);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFracerr",&pdfOnResDK_DPiXFracerr);
   fChain->SetBranchAddress("pdfOnResDK.DPiXFracpull",&pdfOnResDK_DPiXFracpull);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvtserr",&pdfOnResDK_DpiGoodD0NumEvtserr);
   fChain->SetBranchAddress("pdfOnResDK.DpiGoodD0NumEvtspull",&pdfOnResDK_DpiGoodD0NumEvtspull);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodDerr",&pdfOnResDK_asymSigGoodDerr);
   fChain->SetBranchAddress("pdfOnResDK.asymSigGoodDpull",&pdfOnResDK_asymSigGoodDpull);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvtserr",&pdfOnResDK_qqBadD0NumEvtserr);
   fChain->SetBranchAddress("pdfOnResDK.qqBadD0NumEvtspull",&pdfOnResDK_qqBadD0NumEvtspull);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvtserr",&pdfOnResDK_sigGoodD0NumEvtserr);
   fChain->SetBranchAddress("pdfOnResDK.sigGoodD0NumEvtspull",&pdfOnResDK_sigGoodD0NumEvtspull);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvtserr",&pdfOnResDK_totBBNumEvtserr);
   fChain->SetBranchAddress("pdfOnResDK.totBBNumEvtspull",&pdfOnResDK_totBBNumEvtspull);
}

Bool_t anaRndToy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_dalitzHolderN_sigGoodD0_rho = fChain->GetBranch("dalitzHolderN.sigGoodD0.rho");
   b_dalitzHolderN_sigGoodD0_rho_err = fChain->GetBranch("dalitzHolderN.sigGoodD0.rho_err");
   b_dalitzHolderN_sigGoodD0_rho_aerr_lo = fChain->GetBranch("dalitzHolderN.sigGoodD0.rho_aerr_lo");
   b_dalitzHolderN_sigGoodD0_rho_aerr_hi = fChain->GetBranch("dalitzHolderN.sigGoodD0.rho_aerr_hi");
   b_dalitzHolderN_sigGoodD0_theta = fChain->GetBranch("dalitzHolderN.sigGoodD0.theta");
   b_dalitzHolderN_sigGoodD0_theta_err = fChain->GetBranch("dalitzHolderN.sigGoodD0.theta_err");
   b_dalitzHolderN_sigGoodD0_theta_aerr_lo = fChain->GetBranch("dalitzHolderN.sigGoodD0.theta_aerr_lo");
   b_dalitzHolderN_sigGoodD0_theta_aerr_hi = fChain->GetBranch("dalitzHolderN.sigGoodD0.theta_aerr_hi");
   b_dalitzHolderP_sigGoodD0_rho = fChain->GetBranch("dalitzHolderP.sigGoodD0.rho");
   b_dalitzHolderP_sigGoodD0_rho_err = fChain->GetBranch("dalitzHolderP.sigGoodD0.rho_err");
   b_dalitzHolderP_sigGoodD0_rho_aerr_lo = fChain->GetBranch("dalitzHolderP.sigGoodD0.rho_aerr_lo");
   b_dalitzHolderP_sigGoodD0_rho_aerr_hi = fChain->GetBranch("dalitzHolderP.sigGoodD0.rho_aerr_hi");
   b_dalitzHolderP_sigGoodD0_theta = fChain->GetBranch("dalitzHolderP.sigGoodD0.theta");
   b_dalitzHolderP_sigGoodD0_theta_err = fChain->GetBranch("dalitzHolderP.sigGoodD0.theta_err");
   b_dalitzHolderP_sigGoodD0_theta_aerr_lo = fChain->GetBranch("dalitzHolderP.sigGoodD0.theta_aerr_lo");
   b_dalitzHolderP_sigGoodD0_theta_aerr_hi = fChain->GetBranch("dalitzHolderP.sigGoodD0.theta_aerr_hi");
   b_pdfOnResDK_DPiXFrac = fChain->GetBranch("pdfOnResDK.DPiXFrac");
   b_pdfOnResDK_DPiXFrac_err = fChain->GetBranch("pdfOnResDK.DPiXFrac_err");
   b_pdfOnResDK_DPiXFrac_aerr_lo = fChain->GetBranch("pdfOnResDK.DPiXFrac_aerr_lo");
   b_pdfOnResDK_DPiXFrac_aerr_hi = fChain->GetBranch("pdfOnResDK.DPiXFrac_aerr_hi");
   b_pdfOnResDK_DpiGoodD0NumEvts = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvts");
   b_pdfOnResDK_DpiGoodD0NumEvts_err = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvts_err");
   b_pdfOnResDK_DpiGoodD0NumEvts_aerr_lo = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvts_aerr_lo");
   b_pdfOnResDK_DpiGoodD0NumEvts_aerr_hi = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvts_aerr_hi");
   b_pdfOnResDK_asymSigGoodD = fChain->GetBranch("pdfOnResDK.asymSigGoodD");
   b_pdfOnResDK_asymSigGoodD_err = fChain->GetBranch("pdfOnResDK.asymSigGoodD_err");
   b_pdfOnResDK_asymSigGoodD_aerr_lo = fChain->GetBranch("pdfOnResDK.asymSigGoodD_aerr_lo");
   b_pdfOnResDK_asymSigGoodD_aerr_hi = fChain->GetBranch("pdfOnResDK.asymSigGoodD_aerr_hi");
   b_pdfOnResDK_qqBadD0NumEvts = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvts");
   b_pdfOnResDK_qqBadD0NumEvts_err = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvts_err");
   b_pdfOnResDK_qqBadD0NumEvts_aerr_lo = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvts_aerr_lo");
   b_pdfOnResDK_qqBadD0NumEvts_aerr_hi = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvts_aerr_hi");
   b_pdfOnResDK_sigGoodD0NumEvts = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvts");
   b_pdfOnResDK_sigGoodD0NumEvts_err = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvts_err");
   b_pdfOnResDK_sigGoodD0NumEvts_aerr_lo = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvts_aerr_lo");
   b_pdfOnResDK_sigGoodD0NumEvts_aerr_hi = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvts_aerr_hi");
   b_pdfOnResDK_totBBNumEvts = fChain->GetBranch("pdfOnResDK.totBBNumEvts");
   b_pdfOnResDK_totBBNumEvts_err = fChain->GetBranch("pdfOnResDK.totBBNumEvts_err");
   b_pdfOnResDK_totBBNumEvts_aerr_lo = fChain->GetBranch("pdfOnResDK.totBBNumEvts_aerr_lo");
   b_pdfOnResDK_totBBNumEvts_aerr_hi = fChain->GetBranch("pdfOnResDK.totBBNumEvts_aerr_hi");
   b_NLL = fChain->GetBranch("NLL");
   b_NLL_err = fChain->GetBranch("NLL_err");
   b_NLL_aerr_lo = fChain->GetBranch("NLL_aerr_lo");
   b_NLL_aerr_hi = fChain->GetBranch("NLL_aerr_hi");
   b_xyNLL = fChain->GetBranch("xyNLL");
   b_xyNLL_err = fChain->GetBranch("xyNLL_err");
   b_xyNLL_aerr_lo = fChain->GetBranch("xyNLL_aerr_lo");
   b_xyNLL_aerr_hi = fChain->GetBranch("xyNLL_aerr_hi");
   b_rhoPgen = fChain->GetBranch("rhoPgen");
   b_rhoNgen = fChain->GetBranch("rhoNgen");
   b_thetaPgen = fChain->GetBranch("thetaPgen");
   b_thetaNgen = fChain->GetBranch("thetaNgen");
   b_cov_00 = fChain->GetBranch("cov_00");
   b_cov_10 = fChain->GetBranch("cov_10");
   b_cov_11 = fChain->GetBranch("cov_11");
   b_cov_20 = fChain->GetBranch("cov_20");
   b_cov_21 = fChain->GetBranch("cov_21");
   b_cov_22 = fChain->GetBranch("cov_22");
   b_cov_30 = fChain->GetBranch("cov_30");
   b_cov_31 = fChain->GetBranch("cov_31");
   b_cov_32 = fChain->GetBranch("cov_32");
   b_cov_33 = fChain->GetBranch("cov_33");
   b_dalitzHolderN_sigGoodD0_rhoerr = fChain->GetBranch("dalitzHolderN.sigGoodD0.rhoerr");
   b_dalitzHolderN_sigGoodD0_rhopull = fChain->GetBranch("dalitzHolderN.sigGoodD0.rhopull");
   b_dalitzHolderN_sigGoodD0_thetaerr = fChain->GetBranch("dalitzHolderN.sigGoodD0.thetaerr");
   b_dalitzHolderN_sigGoodD0_thetapull = fChain->GetBranch("dalitzHolderN.sigGoodD0.thetapull");
   b_dalitzHolderP_sigGoodD0_rhoerr = fChain->GetBranch("dalitzHolderP.sigGoodD0.rhoerr");
   b_dalitzHolderP_sigGoodD0_rhopull = fChain->GetBranch("dalitzHolderP.sigGoodD0.rhopull");
   b_dalitzHolderP_sigGoodD0_thetaerr = fChain->GetBranch("dalitzHolderP.sigGoodD0.thetaerr");
   b_dalitzHolderP_sigGoodD0_thetapull = fChain->GetBranch("dalitzHolderP.sigGoodD0.thetapull");
   b_pdfOnResDK_DPiXFracerr = fChain->GetBranch("pdfOnResDK.DPiXFracerr");
   b_pdfOnResDK_DPiXFracpull = fChain->GetBranch("pdfOnResDK.DPiXFracpull");
   b_pdfOnResDK_DpiGoodD0NumEvtserr = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvtserr");
   b_pdfOnResDK_DpiGoodD0NumEvtspull = fChain->GetBranch("pdfOnResDK.DpiGoodD0NumEvtspull");
   b_pdfOnResDK_asymSigGoodDerr = fChain->GetBranch("pdfOnResDK.asymSigGoodDerr");
   b_pdfOnResDK_asymSigGoodDpull = fChain->GetBranch("pdfOnResDK.asymSigGoodDpull");
   b_pdfOnResDK_qqBadD0NumEvtserr = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvtserr");
   b_pdfOnResDK_qqBadD0NumEvtspull = fChain->GetBranch("pdfOnResDK.qqBadD0NumEvtspull");
   b_pdfOnResDK_sigGoodD0NumEvtserr = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvtserr");
   b_pdfOnResDK_sigGoodD0NumEvtspull = fChain->GetBranch("pdfOnResDK.sigGoodD0NumEvtspull");
   b_pdfOnResDK_totBBNumEvtserr = fChain->GetBranch("pdfOnResDK.totBBNumEvtserr");
   b_pdfOnResDK_totBBNumEvtspull = fChain->GetBranch("pdfOnResDK.totBBNumEvtspull");

   return kTRUE;
}

#endif // #ifdef anaRndToy_cxx
