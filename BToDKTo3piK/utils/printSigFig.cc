
#include "TString.h"

// Print number with given number of significant digits
TString printSigFig(Double_t num, Int_t sigfig)
{
    if (sigfig <= 0) return TString("");
    if (num==0) return TString("0");

    //get significant digits of integer part of num
    Int_t lognum = (Int_t)log10(fabs(num));

    if (lognum>=0) ++lognum; //adjust for log10 results
    sigfig -= lognum; //subtract for decimal precision calculation

    if(sigfig < 0) sigfig = 0;

    //construct format string
    TString format = "%.";
    format += sigfig;
    format += "f";

    //use format string to get desired results
    TString s;
    s.Form(format.Data(),num);

    return s;
}
