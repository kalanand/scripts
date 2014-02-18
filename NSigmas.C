void NSigmas(const double bkg, const double snl)
{
  int i,N;
  double sum,tmp,tot,high,sigmas,dsum,tmp1,tmp2;
  tot = bkg+snl;
  high = 2*tot;
  sum = TMath::Exp(-tot);
  i = 0;
  dsum = 10000.;
  while ((i<high) || (dsum>1.e-8*sum))
    {
      i++;
      tmp1 = TMath::Poisson(1.*i,tot);
      tmp2 = TMath::Gamma(1.*i+1,bkg);
      dsum = tmp1*tmp2;
      sum+=dsum;
      if ((i==high) || (dsum<1.e-3*sum))
        break;
    }
  sigmas = -sqrt(2)*TMath::ErfInverse(2*sum-1);
  cout<<sum<<" "<<sigmas<<endl;
}
