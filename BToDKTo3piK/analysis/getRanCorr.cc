// gets correlations between the vars:

#include <math.h>

void getRanCorrAll(int todo = 0) {
/*  nnout->setFitMin(0);
  bknnout->setFitMin(0);
  Deltae->setFitMax(0.14);
  d0mass->setFitMin(1.805);
  d0mass->setFitMax(1.924);
  mes->setFitMin(5.2);
  mes->setFitMax(5.3);

  cutGoodD = "Trued0flg>0&&d0recdec==3";
  cutBadD  = "Trued0flg<=0&&d0recdec!=3";
*/

  if (0 == todo || 1 == todo) {
    // DK good D:
    cout<<" ===========D K good D =========="<<endl;
    readCut = cutGoodD;
    data = read(chainRanDK);
    getRanCorr(data);
  }

  if (0 == todo || 2 == todo) {
    // DK bad D:
    cout<<" ===========D K bad D =========="<<endl;
    readCut = cutBadD;
    data = read(chainRanDK);
    getRanCorr(data);
  }

  if (0 == todo || 3 == todo) {
    // Dpi good D:
    cout<<" ============D0 pi Good D ============="<<endl;
    readCut = cutGoodD;
    data = read(chainRanBchDpi);
    getRanCorr(data);
  }

  if (0 == todo || 4 == todo) {
    // Dpi bad D:
    cout<<" ============D0 pi BAD D ============="<<endl;
    readCut = cutBadD;
    data = read(chainRanBchDpi);
    getRanCorr(data);
  }

  if (0 == todo || 5 == todo) {
    cout<<"=================  D0 pi pi( D(*) pi/rho )  ============="<<endl;
    RooAbsData * data0 = 0; 
    TCut DpiOthCut = "((B1decmode>99&&B1decmode<109)||B1decmode>149)||((B2decmode<109&&B2decmode>99)||B2decmode>149)";
    int numDpiOther = addChunk(data0, chainRanBchComb,DpiOthCut+cutBadD, 1);
    numDpiOther += addChunk(data0, chainRanB0Comb,DpiOthCut+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " D pi pi (D(*) pi/rho ) has "<< numDpiOther << " entries."<<endl; 
    getRanCorr(data0);
  }

  if (0 == todo || 6 == todo) {
    // charmless:
    cout<<"=============== charmless( D* K*)  ================="<<endl;
    RooAbsData * data0 = 0;
    // DKother:
    TCut dkCut1 = "((B1decmode>9&&B1decmode<19)||(109<B1decmode&&B1decmode<149))||((B2decmode<19&&B2decmode>9)||(109<B2decmode&&B2decmode<149))";
    int numDKother = addChunk(data0, chainRanBchComb,dkCut1+cutBadD, 1);
    numDKother += addChunk(data0, chainRanB0Comb,dkCut1+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " charm less(D* K*) has "<< numDKother <<" entries" <<endl;
    getRanCorr(data0);
    delete data0;
  }

  
  if (0 == todo || 7 == todo) {
    // BB bad D:
    cout<<"============= BB Bad D =============="<<endl;
    RooAbsData * data1 = 0;
    TCut bbo = "B1decmode<=0&&B2decmode<=0";
    int numBBcomBad = addChunk(data1, chainRanBchComb, bbo+cutBadD, 1);
    numBBcomBad += addChunk(data1, chainRanB0Comb,  bbo+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " BB comb bad has "<<numBBcomBad<<" entries " <<endl;
    getRanCorr(data1);
    delete data1;
  }

  if (0 == todo || 8 == todo) {
    // BB good D:
    cout<<" ==============BB good D ========"<<endl;
    RooAbsData * data2 = 0;
    int numBBgood = addChunk(data2, chainRanBchComb, cutGoodD, 1);
    numBBgood += addChunk(data2, chainRanB0Comb, cutGoodD, (WEIGHT_B0/WEIGHT_BCH));
    cout<< " BB comb good has "<<numBBgood<<" entries " <<endl;
    getRanCorr(data2);
    delete data2;
  }

  if (0 == todo || 9 == todo) {
    // continuum bad D
    cout<<" ============== continuum bad =============="<<endl; 
    RooAbsData * data3 = 0;
    int numQQbad = addChunk(data3, chainRanCc, cutBadD, 1);
    numQQbad += addChunk(data3, chainRanUds, cutBadD, (WEIGHT_UDS/WEIGHT_CC));
    cout << " continuum bad D has " << numQQbad << " entries " <<endl;
    getRanCorr(data3);
    delete data3;
  }

  if (0 == todo || 10 == todo) {
     cout<<"============ continuum good D "<<endl;
    // continuum good D:
    RooAbsData * data4 = 0;
    int numQQgood = addChunk(data4, chainRanCc, cutGoodD, 1);
    numQQgood += addChunk(data4, chainRanUds, cutGoodD, (WEIGHT_UDS/WEIGHT_CC));
    cout << " continuum good D has " << numQQgood << " entries " <<endl;
    getRanCorr(data4);
    delete data4;
  }
}


//----------------------------------------------------------
// Correlations for all vars of this chain:
void getRanCorr(RooAbsData * dat) {

  if (0 == dat || dat->numEntries() == 0) {
    cout << "emptry data set" << endl;
    return;
  }

  data = dat;

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
      corrs[i][j] = getRanCorr(names[i], names[j]);
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
double getRanCorr(const char * nam1, const char * nam2) {

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
