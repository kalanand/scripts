// $Id: runAnaRndToy.cc,v 1.1 2007/05/11 14:06:16 fwinkl Exp $
// Script to run the analysis selector for rndToy
// Simply loads the file and runs the anaRndToy.cc TSelector on it
// All analysis code goes in Process() of anaRndToy

void runAnaRndToy(const char* filename)
{
  TFile f(filename);
  if (f.IsZombie()) {
    cout << "Cannot open file " << filename << endl;
    return;
  }

  RooDataSet* data = (RooDataSet*)f.Get("fitParData");
  TTree& tree = data->tree();

  tree.Process("anaRndToy.cc");
  f.Close();
}
