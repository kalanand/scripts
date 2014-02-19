void ExpectedSignificance(double nb, double ns)
{
  cout<<"Average background    = "<<nb<<endl;
  cout<<"Average signal        = "<<ns<<endl;
  cout<<"Expected Significance = "<<NSigmas(nb,ns)<<" sigma"<<endl; 
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
