// Declares useful global flags. These can be changed from anywhere and
// affect the behavior of all functions that use them:

#ifndef FLAGS_HH
#define FLAGS_HH

Bool_t doFit = kTRUE;      // whether to do a fit (used by fit.cc)
TString fitOption = "mer"; // the fit option (used by fit.cc)
TString optOption = "c";   // "c" is the default option for the const optimizer
Bool_t doPlot = kTRUE;     // if to plot or not

#endif

