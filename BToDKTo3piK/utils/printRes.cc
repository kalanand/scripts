// $Id: printRes.cc,v 1.2 2006/07/03 17:36:24 fwinkl Exp $
// Print a LaTex table of all resonances in amp
// By default only prints amplitudes and phases

void printRes(const BdkDDalitzAmp& amp,
              Bool_t printMass = kFALSE,
              Bool_t printWidth = kFALSE,
              Int_t precision = 1)
{
  RooArgList* ampRes = renameRes(amp.ampRes(),amp);
  RooArgList* phaseRes = renameRes(amp.phaseRes(),amp);
  RooArgList* massRes = renameRes(amp.massRes(),amp);
  RooArgList* gammaRes = renameRes(amp.gammaRes(),amp);

  RooArgList* fitFrac = renameRes(new RooArgList(amp.fitFractions()),amp);

  RooCmdArg cmd1;
  RooCmdArg cmd2;
  if (printMass) cmd1 = Sibling(*massRes);
  if (printWidth) cmd2 = Sibling(*gammaRes);

  ampRes->printLatex(Sibling(*phaseRes),Sibling(*fitFrac),cmd1,cmd2,Format("NE",VerbatimName(true),AutoPrecision(precision)));

}


// Make a copy of list and assign all RooRealVars the resonance name.
// Deletes list!
RooArgList* renameRes(RooArgList* list, const BdkDDalitzAmp& amp)
{  
  RooArgList* newlist = (RooArgList*)list->snapshot(false);
  delete list;

  for (int i=0; i<newlist->getSize(); i++) {
    newlist->at(i)->SetName(amp.nameRes(i));
  }
  return newlist;
}
