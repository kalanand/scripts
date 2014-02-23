// Returns a data set with the requested charge cut:

RooDataSet * chargeCut(RooAbsData * data, int charge = 1) {
  TList * sets = data->split(*Hdtrkchge);
  switch(charge) {
  case 1:
    return (RooDataSet *)sets->At(1);
    break;
  case -1:
    return (RooDataSet *)sets->At(0);
    break;
  }

  cout << "chargeCut() given invalid charge " << charge << endl;  
  return 0;
}
  
