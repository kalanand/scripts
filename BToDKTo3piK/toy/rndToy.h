#ifndef RNDTOY_H
#define RNDTOY_H

typedef struct {
  Double_t rhoPi, rhoPf;
  Double_t rhoNi, rhoNf;
  Double_t thetaPi, thetaPf;
  Double_t thetaNi, thetaNf;
} RndToyData;

const char* branchList = "rhoPi/D:rhoPf:rhoNi:rhoNf:thetaPi:thetaPf:thetaNi:thetaNf";
const Int_t matrixSize = 4;    // size of error matrix

#endif
