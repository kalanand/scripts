// fits the true Dalitz vars and the measured ones to study resolution 
// effect for different values of rho and theta in the range allowed by
// the Kspipi result.

#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooAbsData.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRandom.hh"

#include "BToDKTo3piK/BdkPdfDKDalitz.hh"
#include "BToDKTo3piK/BdkDKDalitz.hh"
#include "BToDKTo3piK/BdkPdf2DpolyDalitz.hh"

#include "../BToDKTo3piK/globals/globals.hh"
#include "../BToDKTo3piK/utils/read.cc"
#include "../BToDKTo3piK/utils/fit.cc"

// forward declarations
RooDataSet * generateFromFlat(const RooAbsData & flat, BdkDKDalitz * pdf);

void setEff(BdkPdfDKDalitz & wrapper, Bool_t applyEff);

void fillFromTruth(const RooArgSet * event, 
		   Bool_t & errMsgPrinted);


// Generate shaped signal from flat signal and fit it, looping over
// the signal PDF parameters:
void fitShapedFlat(const char * outFileName = "fitShapedFlat.dat",
		   double xMStart = -0.022, 
		   double xMStep = 0.063, 
		   int nstepsXM = 3,
		   double yMStart = -0.019,
		   double yMStep = 0.075,
		   int nstepsYM = 3,
		   double xPStart = -0.151, 
		   double xPStep = 0.079, 
		   int nstepsXP = 3,
		   double yPStart = -0.102,
		   double yPStep = 0.069,
		   int nstepsYP = 3,
		   Bool_t toy = kFALSE) {

  char * canName;
  RooAbsData * dataFlatP = 0;
  RooAbsData * dataFlatM = 0;

  if (toy) {
    // test with a truly flat dataset:
    BdkPdf2DpolyDalitz flat("flat", "flat", *m12, *m13);
    flat.parameters().readFromFile("test.par");
    dataFlatP = flat.generate(17500);    
    dataFlatN = flat.generate(17500);    
    canName = "Toy";
  }
  else {
    // test with the real signal MC:
    readCut = cutSigReg+cutDKGoodD;
    dataFlatP = read(sigFlatTree, 0, cutPlus + readCut, allVars);    
    dataFlatN = read(sigFlatTree, 0, cutMinus + readCut, allVars);    
    canName = "Flat Sig MC";
  }  

  // Set all bgd yields to 0:
  pdfOnResDK.parameters().
    readFromFile("../BToDKTo3piK/params/numEvts-sigOnly.par");

  // easy pointer to signal PDF:
  BdkPdfDKDalitz & pdfP = dalitzHolderP.sigGoodD0Type();
  BdkPdfDKDalitz & pdfN = dalitzHolderN.sigGoodD0Type();
  ofstream outFile(outFileName);
  cout << "Will write results to " << outFileName << endl;

  for (int ixM = 0; ixM < nstepsXM; ++ixM) {
    double xM = xMStart + ixM * xMStep;
    for (int iyM = 0; iyM < nstepsYM; ++iyM) {
      double yM = yMStart + iyM * yMStep;
      for (int ixP = 0; ixP < nstepsXP; ++ixP) {
	double xP = xPStart + ixP * xPStep;
	for (int iyP = 0; iyP < nstepsYP; ++iyP) {
	  double yP = yPStart + iyP * yPStep;

	  // remove eff func before shaping data:
	  setEff(pdfP, kFALSE);
	  setEff(pdfN, kFALSE);
	  
	  // Get the polar coords:
	  RooRealVar* thetaMinus = 
	    (RooRealVar*)(pdfOnResDK.parameters().
			  find("dalitzHolderN.sigGoodD0.theta"));
	  
	  RooRealVar* rhoMinus = 
	    (RooRealVar*)(pdfOnResDK.parameters().
			  find("dalitzHolderN.sigGoodD0.rho"));
	  
	  RooRealVar* thetaPlus = 
	    (RooRealVar*)(pdfOnResDK.parameters().
			  find("dalitzHolderP.sigGoodD0.theta"));
	  
	  RooRealVar* rhoPlus = 
	    (RooRealVar*)(pdfOnResDK.parameters().
			  find("dalitzHolderP.sigGoodD0.rho"));

	  // calculate their values given the values of x and y: 
	  const double thetaM = thetaFromCart(xM, yM);
	  const double rhoM = rhoFromCart(xM, yM);
	  
	  const double thetaP = thetaFromCart(xP, yP);
	  const double rhoP = rhoFromCart(xP, yP);

	  // set pdf parameters, recovering from the effects of setEff:
	  thetaMinus->setVal(thetaM);
	  rhoMinus->setVal(rhoM);
	  
	  thetaPlus->setVal(thetaP);
	  rhoPlus->setVal(rhoP);
	  
	  // Set the expected yields:
	  pdfOnResDK.setNsigAsymFromXY();
	  
	  // shape the data:
	  RooDataSet * dSet = generateFromFlat(*dataFlatP, 
					       (BdkDKDalitz*)(pdfP.getPdf()));
	  
	  cout << "Shaped " << dSet->numEntries() << " B+ events" << endl;

	  RooDataSet * dSetN = generateFromFlat(*dataFlatN, 
						(BdkDKDalitz*)(pdfN.getPdf()));
	  
	  cout << "Shaped " << dSetN->numEntries() << " B- events" << endl;

	  dSet->append(*dSetN);
	  
	  cout << "Total=" << dSet->numEntries() << " B+/- events" << endl;

	  // reset eff func before fitting:
	  setEff(pdfP, kTRUE);
	  setEff(pdfN, kTRUE);
	  
	  if (doFit){
	    pdfOnResDK.fit(*dSet, 0.0);
	    RooFitResult * fitRes = pdfOnResDK.xyFitResult();
	    
	    outFile << fitRes->status() << " " 
		    << rhoM << " " 
		    << thetaM << " " 
		    << rhoP << " " 
		    << thetaP << " " 
		    << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
		    << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
		    << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
		    << thetaPlus->getVal() << " " << thetaPlus->getError() << " ";
	    
	    cout << "Wrote to " << outFileName << " "
		 << fitRes->status() << " " 
		 << rhoM << " " 
		 << thetaM << " " 
		 << rhoP << " " 
		 << thetaP << " " 
		 << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
		 << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
		 << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
		 << thetaPlus->getVal() << " " << thetaPlus->getError() << " ";

	    // Fit the MC truth DP vars. 
	    // Now create a dataset for the true Dalitz vars:
	    RooDataSet * trueData = 
	      (RooDataSet*)dSet->emptyClone("true", "true");
	    
	    // and fill it with events from dSet:
	    cout << "True evt" << endl;
	    for (int e = 0; e < dSet->numEntries(); ++e){
	      // Copy the dSet event:
	      const RooArgSet * event = dSet->get(e);
 	      RooArgSet * trueEvt = event->clone(TString(event->GetName()) + "-copy");

	      // Copy the mc values into the "measured" RRVs:
	      const RooRealVar * mass12old = 
		(RooRealVar*)(event->find(mass12mc->GetName()));

	      const RooRealVar * mass13old = 
		(RooRealVar*)(event->find(mass13mc->GetName()));

	      RooRealVar * mass12new = 
		(RooRealVar*)(trueEvt->find(mass12->GetName()));
	      
	      RooRealVar * mass13new = 
		(RooRealVar*)(trueEvt->find(mass13->GetName()));
	      
	      mass12new->setVal(mass12old->getVal());
	      mass13new->setVal(mass13old->getVal());

	      // and put the event into the trueData:
	      trueData->add(*trueEvt);
	    }

	    // Make sure this replacement works:
	    dSet->get(0)->find(mass12->GetName())->Print();
	    dSet->get(0)->find(mass13->GetName())->Print();
	    dSet->get(0)->find(mass12mc->GetName())->Print();
	    dSet->get(0)->find(mass13mc->GetName())->Print();

	    trueData->get(0)->find(mass12->GetName())->Print();
	    trueData->get(0)->find(mass13->GetName())->Print();
	    trueData->get(0)->find(mass12mc->GetName())->Print();
	    trueData->get(0)->find(mass13mc->GetName())->Print();


	    pdfOnResDK.fit(*trueData, 0.0);
	    fitRes = pdfOnResDK.xyFitResult();	    

	    outFile << fitRes->status() << " " 
		    << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
		    << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
		    << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
		    << thetaPlus->getVal() << " " << thetaPlus->getError() 
		    << endl;

	    cout << fitRes->status() << " " 
		 << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
		 << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
		 << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
		 << thetaPlus->getVal() << " " << thetaPlus->getError() 
		 << endl;
	    
	  } // end of if (doFit)
	} // end loop on yP
      } // end loop on xP
    } // end loop on yM
  } //end loop on xM
}

    
    
// Given a flat dataset, generates a dataset weighted according to a given PDF.
// Caller owns the returned dataset:

RooDataSet * generateFromFlat(const RooAbsData & flat, 
			      BdkDKDalitz * pdf) {
  double max = 0;
  double avg = 0;
  int e;
  RooArgSet * nullSet = 0;
  Bool_t errMsg = kFALSE;


  // find max:
  for (e = 0; e < flat.numEntries(); ++e){
    const RooArgSet * event = flat.get(e);

    fillFromTruth(event, errMsg);

    double val = pdf->getVal(nullSet);
    avg += val;
    if (val > max) {
      max = val;
    }
  }

  avg /= flat.numEntries();
  //  cout << "average PDF value = " << avg <<". max PDF value = " << max << endl;

  // init a RDS with the same vars as flat:
  RooDataSet * result = (RooDataSet*)flat.emptyClone(TString(flat.GetName()) + ".shaped",
                                                     TString(flat.GetTitle()) + ".shaped");

  // accept/reject:
  for (e = 0; e < flat.numEntries(); ++e){
    const RooArgSet * event = flat.get(e);

    fillFromTruth(event, errMsg);

    double val = pdf->getVal(nullSet);
    double rand = max * RooRandom::uniform();

    if (val > rand) {
      // then keep the event:
      result->add(*event);
    }
  }

  return result;
}



// fill m12 and m13 with the values of m12mc/m13mc or, if not
// found in the event, with those of m12/m13:
void fillFromTruth(const RooArgSet * event, 
		   Bool_t & errMsgPrinted) {
  // get the true m12, m13:
  RooRealVar * tempM12mc = (RooRealVar*)(event->find(m12mc->GetName()));
  RooRealVar * tempM13mc = (RooRealVar*)(event->find(m13mc->GetName()));
  
  if (0 == tempM12mc || 0 == tempM13mc) {
    // if can't find them, use the measured ones:
    if (kFALSE == errMsgPrinted) {
      cout << "******************" << endl
	   << "** WARNING ** generateFromFlat:"
	   << " Dataset doesn't include m12mc or m13mc. "
	   << " Using m12 and m13 instead." << endl
	   << "******************" << endl;
      errMsgPrinted = kTRUE;
    }
    
    tempM12mc = (RooAbsReal*)(event->find(m12->GetName()));
    tempM13mc = (RooAbsReal*)(event->find(m13->GetName()));
  }

  // fill them into the official m12, m13:
  m12->setVal(tempM12mc->getVal());
  m13->setVal(tempM13mc->getVal());
  
}


// apply or remove the efficiency function, changing the normalization
// coefficients appropriately:
void setEff(BdkPdfDKDalitz & wrapper, Bool_t applyEff) {
  BdkDKDalitz * pdf = (BdkDKDalitz*)(wrapper.getPdf());
  
  if (applyEff) {
    pdf->setEfficiencyFunc((BdkDalitzEff*)eff.getPdf());    
    wrapper.parameters().readFromFile("../BToDKTo3piK/params/dalitzNorm.par");
    cout << "---- setting efficiency function" << endl;
  }
  else {
    // If the PDF has an efficiency function, remove it:
    pdf->setEfficiencyFunc(0);
    wrapper.parameters().readFromFile("../BToDKTo3piK/params/dalitzNormFlatEff.par");
    cout << "---- removing efficiency function" << endl;
  }
}
