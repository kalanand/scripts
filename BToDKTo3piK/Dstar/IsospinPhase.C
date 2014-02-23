const double DEGTORAD = TMath::Pi()/180.0;
const double B_NR     =  0.0885001;
const double B_rhoPlus   =  1.00011;
const double B_rho0   =  0.74043;
const double B_rhoMinus   =  0.84136;
const double B_F      =  0.058538;

void IsospinPhase(){

  double NonRes_amp = 0.089;   
  double NonRes_amp_stat = 0.011;
  double NonRes_amp_syst = 0.012; 
  double NonRes_phase = -11.0*DEGTORAD;
  double NonRes_phase_stat = 4.0*DEGTORAD;
  double NonRes_phase_syst = 2.0*DEGTORAD;

  double Rho770Plus_amp = 0.835;

  double Rho7700_amp = 0.491;
  double Rho7700_amp_stat = 0.005;
  double Rho7700_amp_syst =  0.002;
  double Rho7700_phase =  16.2*DEGTORAD; 
  double Rho7700_phase_stat =  0.6*DEGTORAD;
  double Rho7700_phase_syst =  0.4*DEGTORAD;

  double Rho770Minus_amp = 0.596; 
  double Rho770Minus_amp_stat = 0.007;
  double Rho770Minus_amp_syst =  0.002;
  double Rho770Minus_phase = -2.0*DEGTORAD;
  double Rho770Minus_phase_stat = 0.6*DEGTORAD;
  double Rho770Minus_phase_syst = 0.6*DEGTORAD;

  double Rho1450Plus_amp = 0.033;
  double Rho1450Plus_amp_stat = 0.009;
  double Rho1450Plus_amp_syst = 0.020;
  double Rho1450Plus_phase = -146.0*DEGTORAD; 
  double Rho1450Plus_phase_stat = 18.0*DEGTORAD;
  double Rho1450Plus_phase_syst = 24.0*DEGTORAD;

  double Rho14500_amp = 0.052; 
  double Rho14500_amp_stat = 0.010;
  double Rho14500_amp_syst = 0.006;
  double Rho14500_phase = 10.0*DEGTORAD;
  double Rho14500_phase_stat = 8.0*DEGTORAD;
  double Rho14500_phase_syst = 13.0*DEGTORAD;

  double Rho1450Minus_amp = 0.129; 
  double Rho1450Minus_amp_stat = 0.008; 
  double Rho1450Minus_amp_syst = 0.006;
  double Rho1450Minus_phase = 16.0*DEGTORAD;
  double Rho1450Minus_phase_stat = 3.0*DEGTORAD;
  double Rho1450Minus_phase_syst = 3.0*DEGTORAD;

  double Rho1700Plus_amp = 0.197;
  double Rho1700Plus_amp_stat = 0.016;
  double Rho1700Plus_amp_syst = 0.012;
  double Rho1700Plus_phase = -17.0*DEGTORAD;
  double Rho1700Plus_phase_stat = 2.0*DEGTORAD;
  double Rho1700Plus_phase_syst = 3.0*DEGTORAD;

  double Rho17000_amp = 0.223;
  double Rho17000_amp_stat = 0.013;
  double Rho17000_amp_syst =  0.012;
  double Rho17000_phase =  -17.0*DEGTORAD;
  double Rho17000_phase_stat = 2.0*DEGTORAD;
  double Rho17000_phase_syst = 2.0*DEGTORAD;

  double Rho1700Minus_amp = 0.175;
  double Rho1700Minus_amp_stat =  0.010;
  double Rho1700Minus_amp_syst = 0.006; 
  double Rho1700Minus_phase = -50.0*DEGTORAD;
  double Rho1700Minus_phase_stat = 3.0*DEGTORAD;
  double Rho1700Minus_phase_syst = 3.0*DEGTORAD;

  double F980_amp = 0.047;
  double F980_amp_stat =  0.004;
  double F980_amp_syst =  0.005;
  double F980_phase = -59.0*DEGTORAD; 
  double F980_phase_stat = 5.0*DEGTORAD;
  double F980_phase_syst = 4.0*DEGTORAD;

  double F1370_amp = 0.060;
  double F1370_amp_stat = 0.008;
  double F1370_amp_syst = 0.008;
  double F1370_phase = 156.0*DEGTORAD;
  double F1370_phase_stat = 9.0*DEGTORAD;
  double F1370_phase_syst = 6.0*DEGTORAD;

  double F1500_amp = 0.062;
  double F1500_amp_stat = 0.006;
  double F1500_amp_syst = 0.006;
  double F1500_phase = 12.0*DEGTORAD;
  double F1500_phase_stat = 9.0*DEGTORAD;
  double F1500_phase_syst = 4.0*DEGTORAD;

  double F1710_amp = 0.060;
  double F1710_amp_stat = 0.008;
  double F1710_amp_syst = 0.009;
  double F1710_phase = 51.0*DEGTORAD;
  double F1710_phase_stat = 8.0*DEGTORAD;
  double F1710_phase_syst = 7.0*DEGTORAD;

  double F2_amp = 0.109;
  double F2_amp_stat =  0.004;
  double F2_amp_syst =  0.022;
  double F2_phase = -171.0*DEGTORAD;
  double F2_phase_stat = 3.0*DEGTORAD;
  double F2_phase_syst = 4.0*DEGTORAD; 

  double sigma_amp = 0.087; 
  double sigma_amp_stat = 0.007; 
  double sigma_amp_syst = 0.014;
  double sigma_phase = 8.0*DEGTORAD;
  double sigma_phase_stat = 4.0*DEGTORAD;
  double sigma_phase_syst = 8.0*DEGTORAD;

  // Now calculate the phases

  double M_3_2_phase = NonRes_phase / DEGTORAD;

  double M_2_1_phase 
    = M21Phase(Rho770Plus_amp, 0.0, Rho1450Plus_amp, 
	       Rho1450Plus_phase, Rho1700Plus_amp, Rho1700Plus_phase,
	       Rho7700_amp, Rho7700_phase, Rho14500_amp, 
	       Rho14500_phase, Rho17000_amp, Rho17000_phase,
	       Rho770Minus_amp, Rho770Minus_phase, Rho1450Minus_amp, 
	       Rho1450Minus_phase, Rho1700Minus_amp, Rho1700Minus_phase);

  double M_1_1_phase 
    = M11Phase(Rho770Plus_amp, 0.0, Rho1450Plus_amp, 
	       Rho1450Plus_phase, Rho1700Plus_amp, Rho1700Plus_phase,
	       Rho770Minus_amp, Rho770Minus_phase, Rho1450Minus_amp, 
	       Rho1450Minus_phase, Rho1700Minus_amp, Rho1700Minus_phase);

  double M_0_1_phase 
    = M01Phase(Rho770Plus_amp, 0.0, Rho1450Plus_amp, 
	       Rho1450Plus_phase, Rho1700Plus_amp, Rho1700Plus_phase,
	       Rho7700_amp, Rho7700_phase, Rho14500_amp, 
	       Rho14500_phase, Rho17000_amp, Rho17000_phase,
	       Rho770Minus_amp, Rho770Minus_phase, Rho1450Minus_amp, 
	       Rho1450Minus_phase, Rho1700Minus_amp, Rho1700Minus_phase);

  double M_1_0_phase 
    = M10Phase(F980_amp, F980_phase, F1370_amp, F1370_phase,  
	       F1500_amp, F1500_phase, F1710_amp, F1710_phase, 
	       F2_amp, F2_phase, sigma_amp, sigma_phase);


  cout<< "M_3_2_phase = " << M_3_2_phase << endl;
  cout<< "M_2_1_phase = " << M_2_1_phase << endl;
  cout<< "M_1_1_phase = " << M_1_1_phase << endl;
  cout<< "M_0_1_phase = " << M_0_1_phase << endl;
  cout<< "M_1_0_phase = " << M_1_0_phase << endl;

}




double M21Phase(double p1amp, double p1phase, double p2amp, 
		double p2phase, double p3amp, double p3phase,
		double o1amp, double o1phase, double o2amp, 
		double o2phase, double o3amp, double o3phase,
		double m1amp, double m1phase, double m2amp, 
		double m2phase, double m3amp, double m3phase)
{
  double x = (p1amp*cos(p1phase) + p2amp*cos(p2phase) 
	      + p3amp*cos(p3phase))/B_rhoPlus - 
    (o1amp*cos(o1phase) + o2amp*cos(o2phase) 
     + o3amp*cos(o3phase))*2.0/B_rho0 + 
    (m1amp  * cos(m1phase)+ m2amp * cos(m2phase) 
     + m3amp * cos(m3phase))/B_rhoMinus; 
  
  double y = (p1amp*sin(p1phase) + p2amp*sin(p2phase) 
	      + p3amp*sin(p3phase))/B_rhoPlus - 
    (o1amp*sin(o1phase) + o2amp*sin(o2phase) 
     + o3amp*sin(o3phase))*2.0/B_rho0 + 
    (m1amp  * sin(m1phase)+ m2amp * sin(m2phase) 
     + m3amp * sin(m3phase))/B_rhoMinus; 

  return TMath::ATan(y/x)/ DEGTORAD;
}






double M11Phase(double p1amp, double p1phase, double p2amp, 
		double p2phase, double p3amp, double p3phase,
		double m1amp, double m1phase, double m2amp, 
		double m2phase, double m3amp, double m3phase)
{
  double x = (p1amp*cos(p1phase) + p2amp*cos(p2phase) 
	      + p3amp*cos(p3phase))/B_rhoPlus - 
    (m1amp  * cos(m1phase)+ m2amp * cos(m2phase) 
     + m3amp * cos(m3phase))/B_rhoMinus; 
  
  double y = (p1amp*sin(p1phase) + p2amp*sin(p2phase) 
	      + p3amp*sin(p3phase))/B_rhoPlus - 
    (m1amp  * sin(m1phase)+ m2amp * sin(m2phase) 
     + m3amp * sin(m3phase))/B_rhoMinus; 

  return TMath::ATan(y/x)/ DEGTORAD;
}





double M01Phase(double p1amp, double p1phase, double p2amp, 
		double p2phase, double p3amp, double p3phase,
		double o1amp, double o1phase, double o2amp, 
		double o2phase, double o3amp, double o3phase,
		double m1amp, double m1phase, double m2amp, 
		double m2phase, double m3amp, double m3phase)
{
  double x = (p1amp*cos(p1phase) + p2amp*cos(p2phase) 
	      + p3amp*cos(p3phase))/B_rhoPlus + 
    (o1amp*cos(o1phase) + o2amp*cos(o2phase) 
     + o3amp*cos(o3phase))/B_rho0 + 
    (m1amp  * cos(m1phase)+ m2amp * cos(m2phase) 
     + m3amp * cos(m3phase))/B_rhoMinus; 
  
  double y = (p1amp*sin(p1phase) + p2amp*sin(p2phase) 
	      + p3amp*sin(p3phase))/B_rhoPlus + 
    (o1amp*sin(o1phase) + o2amp*sin(o2phase) 
     + o3amp*sin(o3phase))/B_rho0 + 
    (m1amp  * sin(m1phase)+ m2amp * sin(m2phase) 
     + m3amp * sin(m3phase))/B_rhoMinus; 

  return TMath::ATan(y/x)/ DEGTORAD;
}





double M10Phase(double F980_amp, double F980_phase, double F1370_amp, 
		double F1370_phase, double F1500_amp, 
		double F1500_phase, double F1710_amp, double F1710_phase, 
                double F2_amp, double F2_phase, double sigma_amp, 
		double sigma_phase)
{

  double x = F980_amp*cos(F980_phase) + F1370_amp*cos(F1370_phase) + 
    F1500_amp*cos(F1500_phase) + F1710_amp*cos(F1710_phase) 
    + F2_amp*cos(F2_phase) + sigma_amp*cos(sigma_phase);
  double y = F980_amp*sin(F980_phase) + F1370_amp*sin(F1370_phase) + 
    F1500_amp*sin(F1500_phase) + F1710_amp*sin(F1710_phase) 
    + F2_amp*sin(F2_phase) + sigma_amp*sin(sigma_phase);

  return TMath::ATan(y/x)/ DEGTORAD;
}
