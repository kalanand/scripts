void printX0Y0()
{
  printX0Y0(dalitzHolderN.sigGoodD0Type().dalitzAmp(),false);
  // Since we flip phases in setupAuxPdfs.cc we need to renormalize for KsPiPi
  printX0Y0(ksppAmp,true);  
  printX0Y0(kkpAmp,false);
  printX0Y0(kskkAmp,false);
  printX0Y0(kskpAmp,false);
}


void printX0Y0(BdkAbsDDalitzAmp* amp, Bool_t recalc)
{
  if (!amp) return;

  if (recalc) amp->calDDbarNorm();
  printf("%20s: (x0,y0) = (%.3f,%.3f)\n",amp->GetTitle(),amp->x0(),amp->y0());
} 
