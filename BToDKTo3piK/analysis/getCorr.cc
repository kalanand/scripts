// gets correlations between the vars:

#include <math.h>

void getCorrAll(int todo = 0) {
  if (0 == todo || 1 == todo) {
    // DK good D:
    readCut = cutGoodD;
    getCorr(chainDK);
  }

  if (0 == todo || 2 == todo) {
    // DK bad D:
    readCut = cutBadD;
    getCorr(chainDK);
  }

  if (0 == todo || 3 == todo) {
    // Dpi good D:
    readCut = cutGoodD;
    getCorr(chainDpi);
  }

  if (0 == todo || 4 == todo) {
    // Dpi bad D:
    readCut = cutBadD;
    getCorr(chainDpi);
  }

  if (0 == todo || 5 == todo) {
    // charmless:
    readCut = "";
    getCorr(chainCharmless);
  }

  if (0 == todo || 6 == todo) {
    // BB bad D:
    readCut = cutBadD;
    getCorr(chainBBComb);
  }

  if (0 == todo || 7 == todo) {
    // BB good D:
    readCut = cutGoodD;
    getCorr(chainBBComb);
  }

  if (0 == todo || 8 == todo) {
    // uds:
    readCut = cutBadD;
    getCorr(chainQq);
  }

  if (0 == todo || 9 == todo) {
    // cc bad D:
    readCut = cutGoodD;
    getCorr(chainQq);
  }


  if (0 == todo || 10 == todo) {
    // off resonance:
    readCut = "";
    getCorr(chainOffData);
  }

  if (0 == todo || 11 == todo) {
    //dpipi
    readCut = "";
    getCorr(chainDpipi);
  }
}


//----------------------------------------------------------
// Correlations for all vars of this chain:
void getCorr(TChain & chain) {

  cout << "============================================================"
       << "=== Correlation for chain " << chain->GetTitle()  
       << "============================================================"
       << endl;

  if (chain.GetEntries() == 0) {
    cout << "empty chain" << endl;
    return;
  }

  // Read data:
  data = read(chain);
  if (0 == data || data->numEntries() == 0) {
    cout << "emptry data set" << endl;
    return;
  }

  const int NVARS = 5;
  char * names[NVARS] = {"mes", "Deltae", "d0mass", "nnout", "bknnout"};

  double corrs[NVARS][NVARS] = {{0,0,0,0,0},
				{0,0,0,0,0},
				{0,0,0,0,0},
				{0,0,0,0,0},
				{0,0,0,0,0}};

  int i, j;
  for (i = 0; i < NVARS; ++i){
    for (j = 0; j < NVARS; ++j){
      corrs[i][j] = getCorr(names[i], names[j]);
    }
  }


  cout << "=== variable order: mes Deltae d0mass annout bknnout" << endl;

  for (i = 0; i < NVARS; ++i){
    for (j = 0; j < NVARS; ++j){
      cout << "\t" << corrs[i][j];
    }
    cout << endl;
  }
  cout << endl;
}

//---------------------------------------------------------
// Calculate the correlations:
double getCorr(const char * nam1, const char * nam2) {

  double sum1 = 0;
  double sum2 = 0;
  double sum12 = 0;
  double sum11 = 0;
  double sum22 = 0;

  int size = data->numEntries();
  for (int e = 0; e < size; ++e) {
    RooArgSet * vars = data->get(e);
    double var1 = ((RooRealVar*)vars->find(nam1))->getVal();
    double var2 = ((RooRealVar*)vars->find(nam2))->getVal();

    sum1 += var1;
    sum2 += var2;
    sum12 += var1 * var2;
    sum11 += var1 * var1;
    sum22 += var2 * var2;
    
    /*       cout
    <<    " var1=" << var1
    <<    " var2=" << var2
    <<    " sum1=" << sum1
    <<    " sum2=" << sum2
    <<    " sum12=" << sum12
    <<    " sum11=" << sum11
    <<    " sum22=" << sum22 
    << endl;*/
  }

  sum1 /=  size;
  sum2 /=  size;
  sum12 /= size;
  sum11 /= size;
  sum22 /= size;


  double v12 = sum12 - sum1 * sum2;
  double v11 = sum11 - sum1 * sum1;
  double v22 = sum22 - sum2 * sum2;
  
  //  cout << v12 << " " << v11 << " " << v22 << endl << v12 / sqrt(v11 * v22) << endl;

  return v12 / sqrt(v11 * v22);
}
