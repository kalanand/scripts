// $Id: printOneLine.cc,v 1.2 2006/06/28 02:00:49 fwinkl Exp $
// Print all values for coll in one line 
// format is a list of standard printf format specifiers separated by ":" 
//   e.g.:  "%.4f:%.2f:%1.2f"
// If there are less format specifiers than variables in coll, the last
// specifier will be used for the remaining variables.

TString printOneLine(const RooAbsCollection& coll, TString format,
		     Bool_t printError = kFALSE, Bool_t printLatex = kTRUE)
{
  TString sep = " ";
  TString pm = " +/- ";
  TString enclose = "";
  if (printLatex) sep = " & ";
  if (printLatex) pm = "\\pm";
  if (printLatex) enclose = "$";

  // Get array of format specifiers
  TObjArray* afmt = format.Tokenize(":");

  TString result;
  TIterator* iter = coll.createIterator();  
  RooRealVar* r;
  TString fmt = "";
  int i = 0;
  while (r = (RooRealVar*)iter->Next()) {

    // If there is a format string left, use it
    if (i < afmt->GetEntries()) fmt = ((TObjString*)afmt->At(i))->GetString();

    TString s;    
    if (printError) {
      TString fmtStr = enclose+TString(fmt)+pm+TString(fmt)+enclose;
      s.Form(fmtStr.Data(),r->getVal(),r->getError());
    }
    else {
      TString fmtStr = enclose+TString(fmt)+enclose;
      s.Form(fmtStr.Data(),r->getVal());
    }

    if (i>0) result += sep;
    result += s;
    i++;
  }
  delete iter;
  return result;
}
