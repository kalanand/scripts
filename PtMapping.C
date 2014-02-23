void PtMapping(const double corPt)
{
  TF1 *Response = new TF1("Response","[0]-[1]/(pow(log10(x),[2])+[3])+[4]/x",4,5000);
  Response->SetParameter(0,0.976811);
  Response->SetParameter(1,14.2444);
  Response->SetParameter(2,4.47607);
  Response->SetParameter(3,18.482);
  Response->SetParameter(4,0.717231);  
  double r = Response->Eval(corPt);
  double Pt = corPt*r;
  cout<<"Corrected jet Pt = "<<corPt<<" GeV, Response = "<<r<<", Uncorrected jet Pt = "<<Pt<<" GeV"<<endl;                  
}


