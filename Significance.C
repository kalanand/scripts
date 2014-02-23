void Significance()
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(111); 
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetPalette(1);
  const int NBound = 51;
  const int NBins = 50;
  const double MassBoundaries[NBound] = {80.5137, 95.8941, 113.044, 132.071, 153.087, 176.208, 201.558, 229.267, 259.469, 
                                         292.306, 327.929, 366.493, 408.164, 453.112, 501.520, 553.576, 609.481, 669.442, 733.679, 
                                         802.422, 875.911, 954.400, 1038.15, 1127.45, 1222.58, 1323.84, 1431.57, 1546.09, 1667.76, 
                                         1796.95, 1934.03, 2079.43, 2233.55, 2396.86, 2569.81, 2752.89, 2946.61, 3151.52, 3368.16, 
                                         3597.14, 3839.07, 4094.59, 4364.38, 4649.15, 4949.64, 5266.64, 5600.94, 5953.41, 6324.94, 
                                         6716.47, 7128.96};

  const double bkg[NBins] = {736929,   479034,  59141.6, 30727.8, 134589,  8994.03, 5027.02, 3107.43, 1664.71,  1024.7,
                             31399.3,  20531.7, 12379.1, 8199.17, 4967.13, 3200.73, 2136.93, 53658.8, 34743.2,  22980.2,
                             14958.5,  10346.9, 6819.22, 48064.3, 31524.7, 22221.3, 14628.3, 10424.7, 6855.34,  4808.43,
                             3142.4,   2168.28, 1504.8,  971.011, 688.547, 439.375, 289.499, 193.462, 131.77,   80.5653,
                             51.3084,  32.4631, 19.2031, 11.8424, 7.58899, 4.06498, 2.37414, 1.14408, 0.586682, 0.258657};

  const double bkg10[NBins] = {736929,   479034,   59141.6,  30727.8,  134589,    8994.03,   5027.02,   3107.43,    1664.71,    1024.7,
                               31399.3,  20531.7,  12379.1,  8199.17,  4967.13,   3200.73,   2136.93,   5365.88,    3474.32,    2298.02,
                               1495.85,  1034.69,  681.922,  480.643,  315.247,   222.213,   146.283,   104.247,    68.5534,    48.0843,
                               31.424,   21.6828,  15.048,   9.71011,  6.88547,   4.39375,   2.89499,   1.93462,    1.3177,     0.805653,
                               0.513084, 0.324631, 0.192031, 0.118424, 0.0758899, 0.0406498, 0.0237414, 0.01114408, 0.00586682, 0.00258657};

  const double bkg100[NBins] = {736929,   479034,   59141.6,  30727.8,  134589,    8994.03,   5027.02,   3107.43,    1664.71,    1024.7,
                               31399.3,  20531.7,  12379.1,  8199.17,  4967.13,   3200.73,   2136.93,   53658.8,    34743.2,    22980.2,
                               14958.5,  10346.9,  6819.22,  4806.43,  3152.47,   2222.13,   1462.83,   1042.47,    685.534,    480.843,
                               314.24,   216.828,  150.48,   97.1011,  68.8547,   43.9375,   28.9499,   19.3462,    13.177,     8.05653,
                               5.13084, 3.24631, 1.92031, 1.18424, 0.758899, 0.406498, 0.237414, 0.1114408, 0.0586682, 0.0258657};  
   
  const double snl0[NBins] = {0, 0, 0, 0, 0, 0.374747, 0.0936867, 0.374747, 0.28106, 0.468433, 28.106,  32.7903, 81.9759, 100.713, 
                              140.53, 238.901, 332.588, 20892.1, 12460.3, 2623.23, 1498.99, 655.807, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const double snl0_10[NBins] = {0, 0, 0, 0, 0, 0.374747, 0.0936867, 0.374747, 0.28106, 0.468433, 28.106, 32.7903, 81.9759, 100.713, 
                                 140.53, 238.901, 332.588, 2089.21, 1246.03, 262.323, 149.899, 65.5807, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const double snl0_100[NBins] = {0, 0, 0, 0, 0, 0.374747, 0.0936867, 0.374747, 0.28106, 0.468433, 28.106, 32.7903, 81.9759, 100.713, 
                                  140.53, 238.901, 332.588, 20892.1, 12460.3, 2623.23, 1498.99, 655.807, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  
  const double snl1[NBins] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0588954, 0, 0, 10.6012, 2.35582, 3.53373, 0, 22.3803, 10.6012, 
                              353.373, 235.582, 388.71, 471.163, 1013, 1071.9, 2002.44, 1813.98, 906.99, 176.686, 58.8954, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const double snl1_10[NBins] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0588954, 0, 0, 1.06012, 0.235582, 0.353373, 0, 2.23803, 1.06012, 
                                 3.53373, 2.35582, 3.8871, 4.71163, 10.13, 10.719, 20.0244, 18.1398, 9.0699, 1.76686, 0.588954, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const double snl1_100[NBins] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0588954, 0, 0, 10.6012, 2.35582, 3.53373, 0, 22.3803, 10.6012, 
                                  35.3373, 23.5582, 38.871, 47.1163, 101.3, 107.19, 200.244, 181.398, 90.699, 17.6686, 5.88954, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const double snl2[NBins] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.290504, 0.268157, 0.446929, 0.31285, 
                              0.268157, 0.31285, 0.290504, 0.379889, 0.424582, 0.424582, 0.424582,  0.558661, 0.826818, 0.782125, 1.00559, 1.34079,
                              2.07822, 4.69275, 2.5028, 0.469275, 0.0893857, 0, 0, 0};

  const double snl2_10[NBins] ={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00290504, 0.00268157, 0.00446929,
                                0.0031285, 0.00268157, 0.0031285, 0.00290504, 0.00379889, 0.00424582, 0.00424582, 0.00424582, 0.00558661, 0.00826818,
                                0.00782125, 0.0100559, 0.0134079, 0.0207822, 0.0469275, 0.025028, 0.00469275, 0.000893857, 0, 0, 0};

  const double snl2_100[NBins] ={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0290504, 0.0268157, 0.0446929, 
                                0.031285, 0.0268157, 0.031285, 0.0290504, 0.0379889, 0.0424582, 0.0424582, 0.0424582, 0.0558661, 0.0826818, 0.0782125,
                                0.100559, 0.134079, 0.207822, 0.469275, 0.25028, 0.0469275, 0.00893857, 0, 0, 0};

  int bin,j;
  double snl[3][NBins];
  double background[NBins];
  for(j=0;j<3;j++)
    {
      for(bin=0;bin<NBins;bin++)
        {
	  background[bin] = bkg[bin]; 
          if (j==0)
            snl[j][bin] = snl0[bin];
          if (j==1)
            snl[j][bin] = snl1[bin]; 
          if (j==2)
            snl[j][bin] = snl2[bin];
        }
    }
  double significance[NBins],observed_significance,estimated_significance;
  double b;
  double s;
  double max[3];
  int max_bin[3];
  char Resonance[3][100] = {"700","2000","5000"};
  char name[100];
  TH1F *hSign[3];
  TH1F *hSignMerged[3];
  TH1F *hSignObserved[3];
  TH1F *hSignEstimated[3];
  TCanvas *c[3];
  for(j=0;j<3;j++)
    {
      cout<<"Resonance: "<<Resonance[j]<<" GeV"<<endl;
      max[j] = 0;
      b = 0;
      s = 0;
      sprintf(name,"Resonance_%s_Significance",Resonance[j]);
      hSign[j] = new TH1F(name,name,NBins,MassBoundaries);
      sprintf(name,"Resonance_%s_Significance_Merged",Resonance[j]);
      hSignMerged[j] = new TH1F(name,name,NBins,MassBoundaries); 
      sprintf(name,"Resonance_%s_ObservedSignificance",Resonance[j]);
      hSignObserved[j] = new TH1F(name,name,NBins,MassBoundaries);
      sprintf(name,"Resonance_%s_EstimatedSignificance",Resonance[j]);
      hSignEstimated[j] = new TH1F(name,name,NBins,MassBoundaries);
      
      for(bin=0;bin<NBins;bin++)
        {
          if (snl[j][bin]==0)
            significance[bin] = 0;
          else 
            significance[bin] = NSigmas(background[bin],snl[j][bin]);
          if (significance[bin]>max[j])
            {
              max[j] = significance[bin];
              max_bin[j] = bin; 
            }
          observed_significance = NSigmasObserved(background[bin],snl[j][bin]);
          estimated_significance = NSigmasEstimated(background[bin],snl[j][bin]);
          hSign[j]->SetBinContent(bin+1,significance[bin]); 
          hSignObserved[j]->SetBinContent(bin+1,observed_significance);
          hSignEstimated[j]->SetBinContent(bin+1,estimated_significance); 
        }
      cout<<"bin = "<<max_bin[j]<<", significance = "<<max[j]<<", Mass bin: "<<MassBoundaries[max_bin[j]]<<" - "<<MassBoundaries[max_bin[j]+1]<<endl;
      for(bin=0;bin<3;bin++)
        {
          hSignMerged[j]->SetBinContent(max_bin[j]+bin,hSign[j]->GetBinContent(max_bin[j]+bin));
          b+=background[max_bin[j]+bin-1];
          s+=snl[j][max_bin[j]+bin-1];
          cout<<"Resonance "<<Resonance[j]<<" GeV: "<<max_bin[j]+bin-1<<" "<<b<<" "<<s<<" "<<NSigmas(b,s)<<endl;
        } 
	
      TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
      sprintf(name,"Resonance_%s_SignificanceComparison",Resonance[j]); 
      c[j] = new TCanvas(name,name,900,600);
      hSignObserved[j]->Draw();
      sprintf(name,"Resonance: %s GeV",Resonance[j]);
      hSignObserved[j]->SetTitle(name); 
      hSignObserved[j]->GetXaxis()->SetTitle("Dijet Mass (GeV)");
      hSignObserved[j]->GetYaxis()->SetTitle("Number of #sigma");
      hSignObserved[j]->SetLineWidth(2);
      hSign[j]->SetLineColor(2);
      hSign[j]->SetLineWidth(2);
      hSign[j]->Draw("same");
      hSignMerged[j]->SetLineColor(2);
      hSignMerged[j]->SetFillColor(2);
      hSignMerged[j]->SetFillStyle(3003);
      hSignMerged[j]->Draw("same"); 
      hSignEstimated[j]->SetLineColor(4);
      hSignEstimated[j]->Draw("same");
      hSignEstimated[j]->SetLineWidth(2);
      leg->AddEntry(hSignObserved[j],"#frac{s}{#sqrt{b}}","l");
      leg->AddEntry(hSignEstimated[j],"#frac{s}{#sqrt{b+s}}","l");
      leg->AddEntry(hSign[j],"Expected significance","l");
      leg->AddEntry(hSignMerged[j],"Expected significance: merged bins","F"); 
      leg->SetFillColor(0); 
      leg->Draw();
      TPaveLabel *pave = new TPaveLabel(0.2,0.5,0.4,0.6,"1fb^{-1}","NDC");
      pave->SetFillColor(0);
      pave->Draw();
  
      TPaveText *ptext = new TPaveText(0.7,0.5,0.9,0.6,"NDC");
      sprintf(name,"Single bin significance: %1.1f#sigma",NSigmasObserved(b,s));
      TText *t1 = ptext->AddText(name);
      t1->SetTextColor(1);
      sprintf(name,"Single bin significance: %1.1f#sigma",NSigmasEstimated(b,s));
      TText *t2 = ptext->AddText(name);
      t2->SetTextColor(4);
      sprintf(name,"Single bin significance: %1.1f#sigma",NSigmas(b,s));
      TText *t3 = ptext->AddText(name);
      t3->SetTextColor(2);
      ptext->SetFillColor(0); 
      ptext->Draw();
   }// resonance loop
}

double NSigmasObserved(double bkg, double snl)
{
  double sigmas;
  sigmas = snl/sqrt(bkg);
  return sigmas;
}

double NSigmasEstimated(double bkg, double snl)
{
  double tot,sigmas;
  tot = bkg+snl;
  sigmas = snl/sqrt(tot);
  return sigmas;
}

double NSigmas(double nb, double ns)
{
  int i;
  double sum,n,high,sigmas,dsum,w,p;
  n = nb+ns;
  high = 2*n;
  sum = TMath::Exp(-n);
  i = 0;
  dsum = 1.e+9;
  while ((i<high) || (dsum>1.e-8*sum))
    {
      i++;
      w = TMath::Exp(1.*i*TMath::Log(n)-n-TMath::LnGamma(1.*i+1.));//probability to observe i number of events when you expect n 
      p = TMath::Gamma(1.*i+1.,nb);//tail probability under the background hypothesis
      dsum = p*w;
      sum+=dsum;
      if ((i==high) || (dsum<1.e-8*sum))
        break;
    }
  if (sum<1.e-316)
    {
      sigmas = 0;
      cout<<"WARNING: tail probability out of limits!!!!"<<endl;
    }
  else if (sum>=0.5)
    {
      sigmas = 0;
      cout<<"WARNING: tail probability > 0.5 !!!!"<<endl;
    }
  else  
    sigmas = InvertTailGaus(sum);
  return sigmas;
}

double TailGaus(double x)
{
  // Function to integrate the tail of the Normal Gaussian
  // It uses the standard Trapezoid method.
  double xmax = 38;
  double sum = 0.; 
  double dx = 0.01;
  double y = x;
  double tmp;
  int i,N;
  N = (xmax-x)/dx;
  if (x>xmax)
    return 0.;
  else
    {
      if (x<=0)
        return 0.5;
      else
        { 
          sum = 0.5*(TMath::Gaus(x,0,1,kTRUE)+TMath::Gaus(xmax,0,1,kTRUE));
          for(i=1;i<N;i++)
            {
              tmp = TMath::Gaus(x+i*dx,0,1,kTRUE);
              sum+=tmp;
            }
          return dx*sum;
        }
    } 
}

double InvertTailGaus(double x)
{
  // Function to invert the tail probability of the Normal Gaussian
  // It uses the Bisection method to solve numerically the equation y-f(x)=0 for a given y.
  // The Newton method is faster but it does not work for x<1e-16.
  int n = 0;
  double a = 0;
  double b = 38;
  double y,tmp,result,fa,fb,f,e;
  if (x==0.5)
    result=0;
  else
    {
      e = 100;
      y = 0.5*(a+b);
      while ((e>1.e-3) && (n<100))
        {  
          f = x-TailGaus(y);
          fa = x-TailGaus(a);
          fb = x-TailGaus(b);
          if (f==0)
            {
              result=y;
              break;
            }
          if (fa==0)
            {
              result=a;
              break;
            } 
          if (fb==0)
            {
              result=b;
              break;
            }
          // Do not change the expression to fa*fb<0 because it will work only down to f,fa = 1e-158          
          if ((f>0 && fa<0) || (f<0 && fa>0))
            b = y;
          if ((f>0 && fb<0) || (f<0 && fb>0))
            a = y;
          tmp = y;
          y = 0.5*(a+b); 
          e = fabs(tmp-y);
          n++;
          if ((e<1.e-3) || (n>=100))
            {
              result=y;
              break;
            } 
         }      
    }
  return result; 
}
