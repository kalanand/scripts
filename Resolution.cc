#include "Resolution.h"
#include <cmath>

using std::sqrt;
using std::abs;

Resolution::Resolution() {}

Resolution::~Resolution() {}

bool
Resolution::electronResolution(const double et, const double eta, 
                               double& etRes, double& etaRes, double& phiRes)
{
// Check that eta is in range

  if(abs(eta)>2.5)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.174) {
    aEt  = 0.01188;   bEt  = 0.045;     cEt  = 0.29;
    aEta = 0.0004763; bEta = 0.00059;   cEta = 0.0;
    aPhi = 0.0;       bPhi = 0.0014437; cPhi = 0.0;
  } else if(0.174<=abs(eta) && abs(eta)<0.261) {
    aEt  = 0.01256;   bEt  = 0.0564;   cEt  = 0.0;
    aEta = 0.0003963; bEta = 0.000848; cEta = 0.0;
    aPhi = 8.8e-05;   bPhi = 0.001193; cPhi = 0.0041;
  } else if(0.261<=abs(eta) && abs(eta)<0.348) {
    aEt  = 0.01129;  bEt  = 0.0703;   cEt  = 0.0;
    aEta = 0.000348; bEta = 0.00091;  cEta = 0.0;
    aPhi = 9.5e-05;  bPhi = 0.001192; cPhi = 0.00437;
  } else if(0.348<=abs(eta) && abs(eta)<0.435) {
    aEt  = 0.01275;   bEt  = 0.0621;  cEt  = 0.0;
    aEta = 0.0003152; bEta = 0.00096; cEta = 0.0;
    aPhi = 5.5e-05;   bPhi = 0.00143; cPhi = 0.00293;
  } else if(0.435<=abs(eta) && abs(eta)<0.522) {
    aEt  = 0.01256;   bEt  = 0.0678;   cEt  = 0.0;
    aEta = 0.0003111; bEta = 0.00093;  cEta = 0.0;
    aPhi = 7.4e-05;   bPhi = 0.001391; cPhi = 0.00326;
  } else if(0.522<=abs(eta) && abs(eta)<0.609) {
    aEt  = 0.01139;   bEt  = 0.0729;   cEt  = 0.0;
    aEta = 0.0003167; bEta = 0.00088;  cEta = 0.0;
    aPhi = 0.000114;  bPhi = 0.001294; cPhi = 0.00392;
  } else if(0.609<=abs(eta) && abs(eta)<0.696) {
    aEt  = 0.01285;   bEt  = 0.0599;   cEt  = 0.0;
    aEta = 0.0003251; bEta = 0.00102;  cEta = 0.0;
    aPhi = 7.8e-05;   bPhi = 0.001452; cPhi = 0.00304;
  } else if(0.696<=abs(eta) && abs(eta)<0.783) {
    aEt  = 0.01147;   bEt  = 0.0784;   cEt  = 0.0;
    aEta = 0.0003363; bEta = 0.001;    cEta = 0.0;
    aPhi = 0.000108;  bPhi = 0.001513; cPhi = 0.00293;
  } else if(0.783<=abs(eta) && abs(eta)<0.870) {
    aEt  = 0.01374;  bEt  = 0.0761;   cEt  = 0.0;
    aEta = 0.000324; bEta = 0.00106;  cEta = 0.0;
    aPhi = 0.000127; bPhi = 0.001556; cPhi = 0.00294;
  } else if(0.870<=abs(eta) && abs(eta)<0.957) {
    aEt  = 0.01431;   bEt  = 0.0754;  cEt  = 0.0;
    aEta = 0.0003081; bEta = 0.001;   cEta = 0.0;
    aPhi = 0.000164;  bPhi = 0.00149; cPhi = 0.00411;
  } else if(0.957<=abs(eta) && abs(eta)<1.044) {
    aEt  = 0.01196;   bEt  = 0.1066;   cEt  = 0.0;
    aEta = 0.0003212; bEta = 0.001;    cEta = 0.0;
    aPhi = 0.0001111; bPhi = 0.001933; cPhi = 0.0;
  } else if(1.044<=abs(eta) && abs(eta)<1.131) {
    aEt  = 0.01613;   bEt  = 0.1164; cEt   = 0.0;
    aEta = 0.0003348; bEta = 0.0011; cEta  = 0.0;
    aPhi = 0.000164;  bPhi = 0.00195; cPhi = 0.0022;
  } else if(1.131<=abs(eta) && abs(eta)<1.218) {
    aEt  = 0.0227;    bEt  = 0.1091;  cEt  = 0.0;
    aEta = 0.0003474; bEta = 0.00109; cEta = 0.0;
    aPhi = 0.000191;  bPhi = 0.00216; cPhi = 0.0026;
  } else if(1.218<=abs(eta) && abs(eta)<1.305) {
    aEt  = 0.0158;    bEt  = 0.1718;  cEt  = 0.0;
    aEta = 0.0003354; bEta = 0.00102; cEta = 0.0;
    aPhi = 0.000274;  bPhi = 0.00208; cPhi = 0.0028;
  } else if(1.305<=abs(eta) && abs(eta)<1.392) {
    aEt  = 0.0176;   bEt  = 0.1718;   cEt  = 0.0;
    aEta = 0.000332; bEta = 0.00109;  cEta = 0.0;
    aPhi = 0.000253; bPhi = 0.002472; cPhi = 0.0;
  } else if(1.392<=abs(eta) && abs(eta)<1.479) {
    aEt  = 0.0077;   bEt  = 0.2288;   cEt  = 0.0;
    aEta = 0.000317; bEta = 0.001049; cEta = 0.0;
    aPhi = 0.000285; bPhi = 0.00255;  cPhi = 0.003;
  } else if(1.479<=abs(eta) && abs(eta)<1.653) {
    aEt  = 0.047;     bEt  = 0.158;   cEt  = 0.0;
    aEta = 0.0003479; bEta = 0.0;     cEta = 0.0036;
    aPhi = 0.000333;  bPhi = 0.00277; cPhi = 0.0;
  } else if(1.653<=abs(eta) && abs(eta)<1.740) {
    aEt  = 0.0;     bEt  = 0.2;     cEt  = 0.0;
    aEta = 0.0;     bEta = 0.0;     cEta = 0.0;
    aPhi = 0.00038; bPhi = 0.00282; cPhi = 0.0;
  } else if(1.740<=abs(eta) && abs(eta)<1.830) {
    aEt  = 0.04019;  bEt  = 0.0;     cEt  = 0.0;
    aEta = 0.00033;  bEta = 0.0009;  cEta = 0.0019;
    aPhi = 0.000269; bPhi = 0.00324; cPhi = 0.0;
  } else if(1.830<=abs(eta) && abs(eta)<1.930) {
    aEt  = 0.0;      bEt  = 0.0;     cEt  = 0.0;
    aEta = 0.000348; bEta = 0.00096; cEta = 0.0016;
    aPhi = 0.000271; bPhi = 0.00369; cPhi = 0.0;
  } else if(1.930<=abs(eta) && abs(eta)<2.043) {
    aEt  = 0.038;     bEt  = 0.096;  cEt  = 0.0;
    aEta = 0.0003786; bEta = 0.0;    cEta = 0.00424;
    aPhi = 0.00028;   bPhi = 0.0031; cPhi = 0.0;
  } else if(2.043<=abs(eta) && abs(eta)<2.172) {
    aEt  = 0.0382;   bEt  = 0.076;   cEt  = 0.28;
    aEta = 0.000389; bEta = 0.00106; cEta = 0.0;
    aPhi = 0.000401; bPhi = 0.0025;  cPhi = 0.0114;
  } else if(2.172<=abs(eta) && abs(eta)<2.322) {
    aEt  = 0.035;    bEt  = 0.11;    cEt  = 0.0;
    aEta = 0.000486; bEta = 0.0002;  cEta = 0.0052;
    aPhi = 0.0;      bPhi = 0.00432; cPhi = 0.0088;
  } else if(2.322<=abs(eta) && abs(eta)<2.500) {
    aEt  = 0.0354;   bEt  = 0.123; cEt  = 0.1;
    aEta = 0.000568; bEta = 0.0;   cEta = 0.00734;
    aPhi = 0.000671; bPhi = 0.0;   cPhi = 0.0158;
  } else {
    return false;
  }
  etRes  = et * (sqrt(square(aEt)  + square(bEt/sqrt(et))  + square(cEt/et))),
  etaRes =       sqrt(square(aEta) + square(bEta/sqrt(et)) + square(cEta/et));
  phiRes =       sqrt(square(aPhi) + square(bPhi/sqrt(et)) + square(cPhi/et));
  return true;
}

bool
Resolution::muonResolution(const double et, const double eta, 
                           double& etRes, double& etaRes, double& phiRes)
{ 
// Check that eta is in range

  if(abs(eta)>2.4)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0; 
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0; 

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.100) {
    aEt  = 0.00475;   bEt  = 0.0002365; cEt  = 0.0;
    aEta = 0.0004348; bEta = 0.001063;  cEta = 0.0;
    aPhi = 6.28e-05;  bPhi = 0.0;       cPhi = 0.004545;
  } else if(0.100<=abs(eta) && abs(eta)<0.200) {
    aEt  = 0.00509;   bEt  = 0.0002298; cEt  = 0.0;
    aEta = 0.0004348; bEta = 0.001063;  cEta = 0.0;
    aPhi = 5.53e-05;  bPhi = 0.0;       cPhi = 0.004763;
  } else if(0.200<=abs(eta) && abs(eta)<0.300) {
    aEt  = 0.005942;  bEt  = 0.0002138; cEt  = 0.0;
    aEta = 0.0003412; bEta = 0.000857;  cEta = 0.00147;
    aPhi = 5.39e-5;   bPhi = 0.0;       cPhi = 0.004842;
  } else if(0.300<=abs(eta) && abs(eta)<0.400) {
    aEt  = 0.006989;  bEt  = 0.0002003; cEt  = 0.0;
    aEta = 0.0003208; bEta = 0.000604;  cEta = 0.00187;
    aPhi = 5.63e-5;   bPhi = 0.0;       cPhi = 0.00494;
  } else if(0.400<=abs(eta) && abs(eta)<0.500) {
    aEt  = 0.007227;  bEt  = 0.0001996; cEt  = 0.0;
    aEta = 0.0002908; bEta = 0.000733;  cEta = 0.00151;
    aPhi = 5.58e-5;   bPhi = 0.0;       cPhi = 0.00501;
  } else if(0.500<=abs(eta) && abs(eta)<0.600) {
    aEt  = 0.007528; bEt  = 0.0001935; cEt  = 0.0;
    aEta = 0.000289; bEta = 0.00076;   cEta = 0.00154;
    aPhi = 5.65e-5;  bPhi = 0.0;       cPhi = 0.005082;
  } else if(0.600<=abs(eta) && abs(eta)<0.700) {
    aEt  = 0.007909; bEt  = 0.0001863; cEt  = 0.0;
    aEta = 0.000309; bEta = 0.000667;  cEta = 0.00194;
    aPhi = 5.58e-5;  bPhi = 0.0;       cPhi = 0.005241;
  } else if(0.700<=abs(eta) && abs(eta)<0.800) {
    aEt  = 0.008298;  bEt  = 0.000185; cEt  = 0.0;
    aEta = 0.0002887; bEta = 0.000876; cEta = 0.00179;
    aPhi = 5.97e-5;   bPhi = 0.0;      cPhi = 0.005085;
  } else if(0.800<=abs(eta) && abs(eta)<0.900) {
    aEt  = 0.00918;   bEt  = 0.0001911; cEt  = 0.0;
    aEta = 0.0002956; bEta = 0.000752;  cEta = 0.00208;
    aPhi = 5.9e-5;    bPhi = 0.0;       cPhi = 0.005506;
  } else if(0.900<=abs(eta) && abs(eta)<1.000) {
    aEt  = 0.01096;   bEt  = 0.0001899; cEt  = 0.0;
    aEta = 0.0002734; bEta = 0.000967;  cEta = 0.00134;
    aPhi = 7.48e-5;   bPhi = 0.0;       cPhi = 0.005443;
  } else if(1.000<=abs(eta) && abs(eta)<1.100) {
    aEt  = 0.01262;   bEt  = 0.0001614; cEt  = 0.0;
    aEta = 0.0002831; bEta = 0.000968;  cEta = 0.00166;
    aPhi = 7.81e-5;   bPhi = 0.0;       cPhi = 0.005585;
  } else if(1.100<=abs(eta) && abs(eta)<1.200) {
    aEt  = 0.01379;  bEt  = 0.0001618; cEt  = 0.0;
    aEta = 0.000293; bEta = 0.000942;  cEta = 0.002;
    aPhi = 8.19e-5;  bPhi = 0.0;       cPhi = 0.005921;
  } else if(1.200<=abs(eta) && abs(eta)<1.300) {
    aEt  = 0.01485;   bEt  = 0.0001574; cEt  = 0.0;
    aEta = 0.0002907; bEta = 0.000832;  cEta = 0.002;
    aPhi = 7.89e-5;   bPhi = 0.00039;   cPhi = 0.00593;
  } else if(1.300<=abs(eta) && abs(eta)<1.400) {
    aEt  = 0.0152;    bEt  = 0.0001719; cEt  = 0.0;
    aEta = 0.0002937; bEta = 0.000839;  cEta = 0.00232;
    aPhi = 5.9e-5;    bPhi = 0.000724;  cPhi = 0.005664;
  } else if(1.400<=abs(eta) && abs(eta)<1.500) {
    aEt  = 0.01471;   bEt  = 0.0001828; cEt  = 0.0;
    aEta = 0.0002999; bEta = 0.000864;  cEta = 0.00229;
    aPhi = 4.7e-5;    bPhi = 0.000834;  cPhi = 0.00527;
  } else if(1.500<=abs(eta) && abs(eta)<1.600) {
    aEt  = 0.01337;   bEt  = 0.0002375; cEt  = 0.0;
    aEta = 0.0003035; bEta = 0.000746;  cEta = 0.00258;
    aPhi = 8.16e-5;   bPhi = 0.000757;  cPhi = 0.005558;
  } else if(1.600<=abs(eta) && abs(eta)<1.700) {
    aEt  = 0.01308;   bEt  = 0.000285; cEt  = 0.0;
    aEta = 0.0002967; bEta = 0.000798; cEta = 0.00263;
    aPhi = 6.2e-5;    bPhi = 0.001025; cPhi = 0.00523;
  } else if(1.700<=abs(eta) && abs(eta)<1.800) {
    aEt  = 0.01302;   bEt  = 0.0003797; cEt  = 0.0;
    aEta = 0.0003063; bEta = 0.000776;  cEta = 0.00278;
    aPhi = 0.000107;  bPhi = 0.001011;  cPhi = 0.00554;
  } else if(1.800<=abs(eta) && abs(eta)<1.900) {
    aEt  = 0.0139;    bEt  = 0.000492; cEt  = 0.0;
    aEta = 0.0003285; bEta = 0.00077;  cEta = 0.00292;
    aPhi = 0.000119 ; bPhi = 0.001163; cPhi = 0.00519;
  } else if(1.900<=abs(eta) && abs(eta)<2.000) {
    aEt  = 0.01507;   bEt  = 0.000581; cEt  = 0.0;
    aEta = 0.0003365; bEta = 0.00084;  cEta = 0.00323;
    aPhi = 0.000193;  bPhi = 0.00067;  cPhi = 0.00613;
  } else if(2.000<=abs(eta) && abs(eta)<2.100) {
    aEt  = 0.01711;   bEt  = 0.000731; cEt  = 0.0;
    aEta = 0.0003504; bEta = 0.00078;  cEta = 0.00365;
    aPhi = 0.000217;  bPhi = 0.00121;  cPhi = 0.00558;
  } else if(2.100<=abs(eta) && abs(eta)<2.200) {
    aEt  = 0.01973;  bEt  = 0.000823; cEt  = 0.0;
    aEta = 0.000381; bEta = 0.00088; cEta = 0.00369;
    aPhi = 0.000283; bPhi = 0.00082;  cPhi = 0.00608;
  } else if(2.200<=abs(eta) && abs(eta)<2.300) {
    aEt  = 0.02159;  bEt  = 0.0001052; cEt  = 0.0;
    aEta = 0.00042;  bEta = 0.00097;   cEta = 0.00393;
    aPhi = 0.000304; bPhi = 0.00149;   cPhi = 0.00549;
  } else if(2.300<=abs(eta) && abs(eta)<2.400) {
    aEt  = 0.02155;  bEt  = 0.001346; cEt  = 0.0;
    aEta = 0.000403; bEta = 0.00153;  cEta = 0.00403;
    aPhi = 0.000331; bPhi = 0.00183;  cPhi = 0.00585;
  } else {
    return false;
  }
  etRes  = et * (aEt + bEt * et);
  etaRes = sqrt(square(aEta) + square(bEta/sqrt(et)) + square(cEta/et));
  phiRes = sqrt(square(aPhi) + square(bPhi/sqrt(et)) + square(cPhi/et));
  return true;
}

bool
Resolution::caloMETResolution(const double et, double& etRes, double& etaRes, double& phiRes)
{
  etRes  = et * (sqrt(square(1.462/sqrt(et)) + square(18.19/et)));
  etaRes = 0.0; 
  phiRes =       sqrt(square(1.237/sqrt(et)) + square(18.702/et));
  return true;
}

bool
Resolution::PFMETResolution(const double et, double& etRes, double& etaRes, double& phiRes)
{
  etRes  = et * (sqrt(square(0.05469) +                          square(10.549/et)));
  etaRes =       0.0;
  phiRes =       sqrt(                  square(0.164/sqrt(et)) + square(11.068/et));
  return true;
}


bool
Resolution::udscCaloJetResolution(const double et, const double eta,
                                  double& etRes, double& etaRes, double& phiRes)
{ 
// Check that eta is in range

  if(abs(eta)>3.0)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.087) {
    aEt  = 0.031; bEt  = 1.236; cEt  = 4.44;
    aEta = 0.00836; bEta = 0.0; cEta = .4036;
    aPhi = 0.00858; bPhi = 0.0; cPhi = 2.475;
  } else if(0.087<=abs(eta) && abs(eta)<0.174) {
    aEt  = 0.0446; bEt  = 1.185; cEt  = 5.03;
    aEta = 0.00792; bEta = 0.0; cEta = 1.4432;
    aPhi = 0.00734; bPhi = 0.0; cPhi = 2.547;
  } else if(0.174<=abs(eta) && abs(eta)<0.261) {
    aEt  = 0.0478; bEt  = 1.172; cEt  = 5.23;
    aEta = 0.00807; bEta = 0.0; cEta = 1.4603;
    aPhi = 0.00912; bPhi = 0.0; cPhi = 2.502;
  } else if(0.261<=abs(eta) && abs(eta)<0.348) {
    aEt  = 0.0438; bEt  = 1.169; cEt  = 5.21;
    aEta = 0.00755; bEta = 0.0; cEta = 1.4781;
    aPhi = 0.00742; bPhi = 0.0; cPhi = 2.513;
  } else if(0.348<=abs(eta) && abs(eta)<0.435) {
    aEt  = 0.0443; bEt  = 1.163; cEt  = 5.14;
    aEta = 0.00772; bEta = 0.0; cEta = 1.5064;
    aPhi = 0.00828; bPhi = 0.0; cPhi = 2.529;
  } else if(0.435<=abs(eta) && abs(eta)<0.522) {
    aEt  = 0.0499; bEt  = 1.142; cEt  = 5.06;
    aEta = 0.00793; bEta = 0.0; cEta = 1.4902;
    aPhi = 0.00676; bPhi = 0.0; cPhi = 2.534;
  } else if(0.522<=abs(eta) && abs(eta)<0.609) {
    aEt  = 0.0536; bEt  = 1.121; cEt  = 5.24;
    aEta = 0.00803; bEta = 0.0; cEta = 1.4472;
    aPhi = 0.00659; bPhi = 0.0; cPhi = 2.498;
  } else if(0.609<=abs(eta) && abs(eta)<0.696) {
    aEt  = 0.0487; bEt  = 1.129; cEt  = 5.26;
    aEta = 0.00831; bEta = 0.0; cEta = 1.4409;
    aPhi = 0.00812; bPhi = 0.0; cPhi = 2.465;
  } else if(0.696<=abs(eta) && abs(eta)<0.783) {
    aEt  = 0.0434; aEt  = 1.194; aEt  = 4.64;
    aEta = 0.00844; bEta = 0.0; cEta = 1.4536;
    aPhi = 0.00706; bPhi = 0.0; cPhi = 2.504;
  } else if(0.783<=abs(eta) && abs(eta)<0.870) {
    aEt  = 0.0447; bEt  = 1.23; cEt  = 4.37;
    aEta = 0.00777; bEta = 0.0; cEta = 1.5148;
    aPhi = 0.00688; bPhi = 0.0; cPhi = 2.535;
  } else if(0.870<=abs(eta) && abs(eta)<0.957) {
    aEt  = 0.0383; bEt  = 1.263; cEt  = 4.45;
    aEta = 0.00753; bEta = 0.0; cEta = 1.5043;
    aPhi = 0.00698; bPhi = 0.0; cPhi = 2.512;
  } else if(0.957<=abs(eta) && abs(eta)<1.044) {
    aEt  = 0.0471; bEt  = 1.198; cEt  = 5.1;
    aEta = 0.00756; bEta = 0.0; cEta = 1.5162;
    aPhi = 0.00731; bPhi = 0.0; cPhi = 2.519;
  } else if(1.044<=abs(eta) && abs(eta)<1.131) {
    aEt  = 0.0485; bEt  = 1.245; cEt  = 4.88;
    aEta = 0.00737; bEta = 0.0; cEta = 1.5445;
    aPhi = 0.00755; bPhi = 0.0; cPhi = 2.526;
  } else if(1.131<=abs(eta) && abs(eta)<1.218) {
    aEt  = 0.043; bEt  = 1.271; cEt  = 5.0;
    aEta = 0.00779; bEta = 0.0; cEta = 1.56;
    aPhi = 0.00668; bPhi = 0.0; cPhi = 2.574;
  } else if(1.218<=abs(eta) && abs(eta)<1.305) {
    aEt  = 0.0361; bEt  = 1.323; cEt  = 4.63;
    aEta = 0.0084; bEta = 0.0; cEta = 1.622;
    aPhi = 0.0073; bPhi = 0.0; cPhi = 2.61;
  } else if(1.305<=abs(eta) && abs(eta)<1.392) {
    aEt  = 0.0449; bEt  = 1.319; cEt  = 5.24;
    aEta = 0.01231; bEta = 0.0; cEta = 1.653;
    aPhi = 0.00773; bPhi = 0.0; cPhi = 2.646;
  } else if(1.392<=abs(eta) && abs(eta)<1.479) {
    aEt  = 0.0; bEt  = 1.423; cEt  = 4.42;
    aEta = 0.01187; bEta = 0.0; cEta = 1.668;
    aPhi = 0.00789; bPhi = 0.0; cPhi = 2.823;
  } else if(1.479<=abs(eta) && abs(eta)<1.566) {
    aEt  = 0.0; bEt  = 1.341; cEt  = 5.48;
    aEta = 0.01267; bEta = 0.0; cEta = 1.647;
    aPhi = 0.0084; bPhi = 0.0; cPhi = 2.813;
  } else if(1.566<=abs(eta) && abs(eta)<1.653) {
    aEt  = 0.0; bEt  = 1.242; cEt  = 5.75;
    aEta = 0.00941; bEta = 0.0; cEta = 1.584;
    aPhi = 0.00523; bPhi = 0.0; cPhi = 2.672;
  } else if(1.653<=abs(eta) && abs(eta)<1.740) {
    aEt  = 0.0; bEt  = 1.1864; cEt  = 5.461;
    aEta = 0.00891; bEta = 0.0; cEta = 1.647;
    aPhi = 0.00773; bPhi = 0.0; cPhi = 2.487;
  } else if(1.740<=abs(eta) && abs(eta)<1.830) {
    aEt  = 0.028; bEt  = 1.115; cEt  = 5.5;
    aEta = 0.01023; bEta = 0.0; cEta = 1.649;
    aPhi = 0.00953; bPhi = 0.0; cPhi = 2.394;
  } else if(1.830<=abs(eta) && abs(eta)<1.930) {
    aEt  = 0.016; bEt  = 1.101; cEt  = 4.92;
    aEta = 0.01151; bEta = 0.0; cEta = 1.535;
    aPhi = 0.01088; bPhi = 0.0; cPhi = 2.223;
  } else if(1.930<=abs(eta) && abs(eta)<2.043) {
    aEt  = 0.0396; bEt  = 0.915; cEt  = 5.11;
    aEta = 0.00989; bEta = 0.0; cEta = 1.511;
    aPhi = 0.01146; bPhi = 0.0; cPhi = 2.071;
  } else if(2.043<=abs(eta) && abs(eta)<2.172) {
    aEt  = 0.032; bEt  = 0.907; cEt  = 4.44;
    aEta = 0.01029; bEta = 0.0; cEta = 1.495;
    aPhi = 0.01175; bPhi = 0.0; cPhi = 1.939;
  } else if(2.172<=abs(eta) && abs(eta)<2.322) {
    aEt  = 0.0347; bEt  = 0.875; cEt  = 3.96;
    aEta = 0.01098; bEta = 0.0; cEta = 1.428;
    aPhi = 0.01079; bPhi = 0.0; cPhi = 1.827;
  } else if(2.322<=abs(eta) && abs(eta)<2.500) {
    aEt  = 0.0199; bEt  = 0.851; cEt  = 3.36;
    aEta = 0.01314; bEta = 0.0; cEta = 1.43;
    aPhi = 0.01029; bPhi = 0.0; cPhi = 1.745;
  } else if(2.500<=abs(eta) && abs(eta)<3.000) {
    aEt  = 0.05; bEt  = 0.763; cEt  = 2.99;
    aEta = 0.02238; bEta = 0.0; cEta = 1.612;
    aPhi = 0.01396; bPhi = 0.0; cPhi = 1.5799;
  } else {
    return false;
  }
  etRes  = et * (sqrt(square(aEt) + square(bEt/sqrt(et)) + square(cEt/et)));
  etaRes =       sqrt(square(aEta)                       + square(cEta/et));
  phiRes =       sqrt(square(aPhi)                       + square(cPhi/et));
  return true;
}

bool
Resolution::udscPFJetResolution(const double et, const double eta,
                                double& etRes, double& etaRes, double& phiRes)
{
// Check that eta is in range

  if(abs(eta)>3.0)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.087) {
    aEt  = 0.0642; bEt  = 0.952; cEt  = 0.0;
    aEta = 0.00757; aEta = 0.0; cEta = 1.2578;
    aPhi = 0.01003; bPhi = 0.0; cPhi = 1.3972;
  } else if(0.087<=abs(eta) && abs(eta)<0.174) {
    aEt  = 0.069; bEt  = 0.9303; cEt  = 0.0;
    aEta = 0.0071; bEta = 0.0; cEta = 1.2661;
    aPhi = 0.01; bPhi = 0.0; cPhi = 1.3886;
  } else if(0.174<=abs(eta) && abs(eta)<0.261) {
    aEt  = 0.0675; bEt  = 0.938; cEt  = 0.8;
    aEta = 0.00795; bEta = 0.0; cEta = 1.2713;
    aPhi = 0.01017; aPhi = 0.0; aPhi = 1.4;
  } else if(0.261<=abs(eta) && abs(eta)<0.348) {
    aEt  = 0.0645; bEt  = 0.9409; cEt  = 0.0;
    aEta = 0.00729; bEta = 0.0; cEta = 1.2924;
    aPhi = 0.01004; bPhi = 0.0; cPhi = 1.39;
  } else if(0.348<=abs(eta) && abs(eta)<0.435) {
    aEt  = 0.0616; bEt  = 0.9614; cEt  = 0.0;
    aEta = 0.00689; bEta = 0.0; cEta = 1.3078;
    aPhi = 0.01024; bPhi = 0.0; cPhi = 1.4013;
  } else if(0.435<=abs(eta) && abs(eta)<0.522) {
    aEt  = 0.0708; bEt  = 0.896; cEt  = 1.34;
    aEta = 0.00716; bEta = 0.0; cEta = 1.3051;
    aPhi = 0.00976; bPhi = 0.0; cPhi = 1.4023;
  } else if(0.522<=abs(eta) && abs(eta)<0.609) {
    aEt  = 0.0647; bEt  = 0.9395; cEt  = 0.0;
    aEta = 0.00783; bEta = 0.0; cEta = 1.2687;
    aPhi = 0.00997; bPhi = 0.0; cPhi = 1.3834;
  } else if(0.609<=abs(eta) && abs(eta)<0.696) {
    aEt  = 0.0626; bEt  = 0.9445; cEt  = 0.0;
    aEta = 0.00782; bEta = 0.0; cEta = 1.2664;
    aPhi = 0.00952; bPhi = 0.0; cPhi = 1.4145;
  } else if(0.696<=abs(eta) && abs(eta)<0.783) {
    aEt  = 0.0642; bEt  = 0.9575; cEt  = 0.0;
    aEta = 0.00768; bEta = 0.0; cEta = 1.2863;
    aPhi = 0.0098; bPhi = 0.0; cPhi = 1.4062;
  } else if(0.783<=abs(eta) && abs(eta)<0.870) {
    aEt  = 0.0625; bEt  = 0.9851; cEt  = 0.0;
    aEta = 0.0071; bEta = 0.0; cEta = 1.3159;
    aPhi = 0.01023; bPhi = 0.0; cPhi = 1.4147;
  } else if(0.870<=abs(eta) && abs(eta)<0.957) {
    aEt  = 0.0617; bEt  = 1.0112; cEt  = 0.0;
    aEta = 0.00865; bEta = 0.0; cEta = 1.2837;
    aPhi = 0.01041; bPhi = 0.0; cPhi = 1.4286;
  } else if(0.957<=abs(eta) && abs(eta)<1.044) {
    aEt  = 0.0647; bEt  = 1.026; cEt  = 0.0;
    aEta = 0.0082; bEta = 0.0; cEta = 1.3122;
    aPhi = 0.01049; bPhi = 0.0; cPhi = 1.4245;
  } else if(1.044<=abs(eta) && abs(eta)<1.131) {
    aEt  = 0.0636; bEt  = 1.0591; cEt  = 0.0;
    aEta = 0.00828; bEta = 0.0; cEta = 1.3265;
    aPhi = 0.01083; bPhi = 0.0; cPhi = 1.4504;
  } else if(1.131<=abs(eta) && abs(eta)<1.218) {
    aEt  = 0.0661; bEt  = 1.0793; cEt  = 0.0;
    aEta = 0.00807; bEta = 0.0; cEta = 1.3559;
    aPhi = 0.01091; bPhi = 0.0; cPhi = 1.487;
  } else if(1.218<=abs(eta) && abs(eta)<1.305) {
    aEt  = 0.0614; bEt  = 1.1195; cEt  = 0.0;
    aEta = 0.01007; bEta = 0.0; cEta = 1.3581;
    aPhi = 0.01145; bPhi = 0.0; cPhi = 1.5019;
  } else if(1.305<=abs(eta) && abs(eta)<1.392) {
    aEt  = 0.0654; bEt  = 1.165; cEt  = 0.0;
    aEta = 0.014; bEta = 0.0; cEta = 1.327;
    aPhi = 0.01387; bPhi = 0.0; cPhi = 1.529;
  } else if(1.392<=abs(eta) && abs(eta)<1.479) {
    aEt  = 0.0575; bEt  = 1.205; cEt  = 0.0;
    aEta = 0.01072; bEta = 0.0; cEta = 1.348;
    aPhi = 0.01462; bPhi = 0.0; cPhi = 1.58;
  } else if(1.479<=abs(eta) && abs(eta)<1.566) {
    aEt  = 0.0469; bEt  = 1.19; cEt  = 0.0;
    aEta = 0.00992; bEta = 0.0; cEta = 1.395;
    aPhi = 0.01256; bPhi = 0.0; cPhi = 1.584;
  } else if(1.566<=abs(eta) && abs(eta)<1.653) {
    aEt  = 0.0; bEt  = 1.1632; cEt  = 0.0;
    aEta = 0.00975; bEta = 0.0; cEta = 1.396;
    aPhi = 0.01066; bPhi = 0.0; cPhi = 1.577;
  } else if(1.653<=abs(eta) && abs(eta)<1.740) {
    aEt  = 0.0; bEt  = 1.1109; cEt  = 0.0;
    aEta = 0.00967; bEta = 0.0; cEta = 1.365;
    aPhi = 0.01087; bPhi = 0.0; cPhi = 1.521;
  } else if(1.740<=abs(eta) && abs(eta)<1.830) {
    aEt  = 0.0; bEt  = 1.0841; cEt  = 0.0;
    aEta = 0.0093; bEta = 0.0; cEta = 1.405;
    aPhi = 0.01066; bPhi = 0.0; cPhi = 1.505;
  } else if(1.830<=abs(eta) && abs(eta)<1.930) {
    aEt  = 0.0; bEt  = 1.0288; cEt  = 0.0;
    aEta = 0.01057; bEta = 0.0; cEta = 1.365;
    aPhi = 0.01141; bPhi = 0.0; cPhi = 1.456;
  } else if(1.930<=abs(eta) && abs(eta)<2.043) {
    aEt  = 0.0; bEt  = 0.9821; cEt  = 0.0;
    aEta = 0.00992; bEta = 0.0; cEta = 1.329;
    aPhi = 0.01042; bPhi = 0.0; cPhi = 1.468;
  } else if(2.043<=abs(eta) && abs(eta)<2.172) {
    aEt  = 0.0; bEt  = 0.9441; cEt  = 0.0;
    aEta = 0.00938; bEta = 0.0; cEta = 1.327;
    aPhi = 0.01119; bPhi = 0.0; cPhi = 1.45;
  } else if(2.172<=abs(eta) && abs(eta)<2.322) {
    aEt  = 0.0; bEt  = 0.9134; cEt  = 0.0;
    aEta = 0.00973; bEta = 0.0; cEta = 1.312;
    aPhi = 0.01128; bPhi = 0.0; cPhi = 1.413;
  } else if(2.322<=abs(eta) && abs(eta)<2.500) {
    aEt  = 0.0; bEt  = 0.8322; cEt  = 2.0069;
    aEta = 0.01161; bEta = 0.0; cEta = 1.423;
    aPhi = 0.01256; bPhi = 0.0; cPhi = 1.471;
  } else if(2.500<=abs(eta) && abs(eta)<3.000) {
    aEt  = 0.0526; bEt  = 0.774; cEt  = 2.39;
    aEta = 0.0; bEta = 0.0; cEta = 1.4;
    aPhi = 0.02829; bPhi = 0.0; cPhi = 1.498;
  } else {
    return false;
  }
  etRes  = et * (sqrt(square(aEt) + square(bEt/sqrt(et)) + square(cEt/et)));
  etaRes =       sqrt(square(aEta)                       + square(cEta/et));
  phiRes =       sqrt(square(aPhi)                       + square(cPhi/et));
  return true;
}

bool
Resolution::bPFJetResolution(const double et, const double eta,
                             double& etRes, double& etaRes, double& phiRes)
{ 
// Check that eta is in range

  if(abs(eta)>3.0)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.087) {
    aEt  = 0.0876; bEt  = 0.93; cEt  = 0.0;
    aEta = 0.00658; bEta = 0.0; cEta = 1.3618;
    aPhi = 0.00914; bPhi = 0.0; cPhi = 1.5326;
  } else if(0.087<=abs(eta) && abs(eta)<0.174) {
    aEt  = 0.0892; bEt  = 0.905; cEt  = 1.6;
    aEta = 0.00578; bEta = 0.0; cEta = 1.3927;
    aPhi = 0.0091; bPhi = 0.0; cPhi = 1.5446;
  } else if(0.174<=abs(eta) && abs(eta)<0.261) {
    aEt  = 0.0856; bEt  = 0.946; cEt  = 0.2;
    aEta = 0.0063; bEta = 0.0; cEta = 1.3873;
    aPhi = 0.00892; bPhi = 0.0; cPhi = 1.5446;
  } else if(0.261<=abs(eta) && abs(eta)<0.348) {
    aEt  = 0.0838; bEt  = 0.911; cEt  = 1.76;
    aEta = 0.00587; bEta = 0.0; cEta = 1.4045;
    aPhi = 0.00889; bPhi = 0.0; cPhi = 1.5435;
  } else if(0.348<=abs(eta) && abs(eta)<0.435) {
    aEt  = 0.0792; bEt  = 0.961; cEt  = 0.5;
    aEta = 0.00562; bEta = 0.0; cEta = 1.4079;
    aPhi = 0.00883; bPhi = 0.0; cPhi = 1.54;
  } else if(0.435<=abs(eta) && abs(eta)<0.522) {
    aEt  = 0.0791; bEt  = 0.0; cEt  = 0.9;
    aEta = 0.00602; bEta = 0.0; cEta = 1.4112;
    aPhi = 0.00846; bPhi = 0.0; cPhi = 1.5708;
  } else if(0.522<=abs(eta) && abs(eta)<0.609) {
    aEt  = 0.0748; bEt  = 0.98; cEt  = 0.4;
    aEta = 0.00616; bEta = 0.0; cEta = 1.4132;
    aPhi = 0.00836; bPhi = 0.0; cPhi = 1.5673;
  } else if(0.609<=abs(eta) && abs(eta)<0.696) {
    aEt  = 0.0753; bEt  = 0.969; cEt  = 0.0;
    aEta = 0.00664; bEta = 0.0; cEta = 1.3955;
    aPhi = 0.00826; bPhi = 0.0; cPhi = 1.588;
  } else if(0.696<=abs(eta) && abs(eta)<0.783) {
    aEt  = 0.0831; bEt  = 0.947; cEt  = 0.0;
    aEta = 0.00591; bEta = 0.0; cEta = 1.4045;
    aPhi = 0.00886; bPhi = 0.0; cPhi = 1.561;
  } else if(0.783<=abs(eta) && abs(eta)<0.870) {
    aEt  = 0.0781; bEt  = 0.961; cEt  = 1.16;
    aEta = 0.00683; bEta = 0.0; cEta = 1.3992;
    aPhi = 0.00811; bPhi = 0.0; cPhi = 1.583;
  } else if(0.870<=abs(eta) && abs(eta)<0.957) {
    aEt  = 0.078; bEt  = 1.004; cEt  = 0.7;
    aEta = 0.00695; bEta = 0.0; cEta = 1.425;
    aPhi = 0.00865; bPhi = 0.0; cPhi = 1.582;
  } else if(0.957<=abs(eta) && abs(eta)<1.044) {
    aEt  = 0.0787; bEt  = 1.025; cEt  = 0.0;
    aEta = 0.00618; bEta = 0.0; cEta = 1.452;
    aPhi = 0.00866; bPhi = 0.0; cPhi = 1.619;
  } else if(1.044<=abs(eta) && abs(eta)<1.131) {
    aEt  = 0.081; bEt  = 1.035; cEt  = 0.0;
    aEta = 0.00675; bEta = 0.0; cEta = 1.459;
    aPhi = 0.0087; bPhi = 0.0; cPhi = 1.613;
  } else if(1.131<=abs(eta) && abs(eta)<1.218) {
    aEt  = 0.0853; bEt  = 1.048; cEt  = 0.0;
    aEta = 0.00738; bEta = 0.0; cEta = 1.489;
    aPhi = 0.00942; bPhi = 0.0; cPhi = 1.644;
  } else if(1.218<=abs(eta) && abs(eta)<1.305) {
    aEt  = 0.0875; bEt  = 1.04; cEt  = 0.0;
    aEta = 0.00873; bEta = 0.0; cEta = 1.49;
    aPhi = 0.0094; bPhi = 0.0; cPhi = 1.68;
  } else if(1.305<=abs(eta) && abs(eta)<1.392) {
    aEt  = 0.0906; bEt  = 1.081; cEt  = 0.0;
    aEta = 0.01038; bEta = 0.0; cEta = 1.495;
    aPhi = 0.01143; bPhi = 0.0; cPhi = 1.701;
  } else if(1.392<=abs(eta) && abs(eta)<1.479) {
    aEt  = 0.0919; bEt  = 1.096; cEt  = 0.0;
    aEta = 0.00822; bEta = 0.0; cEta = 1.537;
    aPhi = 0.011; bPhi = 0.0; cPhi = 1.785;
  } else if(1.479<=abs(eta) && abs(eta)<1.566) {
    aEt  = 0.0825; bEt  = 1.124; cEt  = 0.0;
    aEta = 0.00871; bEta = 0.0; cEta = 1.537;
    aPhi = 0.01065; bPhi = 0.0; cPhi = 1.786;
  } else if(1.566<=abs(eta) && abs(eta)<1.653) {
    aEt  = 0.0504; bEt  = 1.174; cEt  = 0.0;
    aEta = 0.00644; bEta = 0.0; cEta = 1.575;
    aPhi = 0.00833; bPhi = 0.00; cPhi = 1.77;
  } else if(1.653<=abs(eta) && abs(eta)<1.740) {
    aEt  = 0.0432; bEt  = 1.122; cEt  = 0.0;
    aEta = 0.00791; bEta = 0.0; cEta = 1.545;
    aPhi = 0.00841; bPhi = 0.0; cPhi = 1.712;
  } else if(1.740<=abs(eta) && abs(eta)<1.830) {
    aEt  = 0.0244; bEt  = 1.113; cEt  = 0.0;
    aEta = 0.00574; bEta = 0.0; cEta = 1.578;
    aPhi = 0.00697; bPhi = 0.0; cPhi = 1.702;
  } else if(1.830<=abs(eta) && abs(eta)<1.930) {
    aEt  = 0.0303; bEt  = 1.067; cEt  = 0.0;
    aEta = 0.00727; bEta = 0.0; cEta = 1.552;
    aPhi = 0.00675; bPhi = 0.0; cPhi = 1.672;
  } else if(1.930<=abs(eta) && abs(eta)<2.043) {
    aEt  = 0.0193; bEt  = 1.052; cEt  = 0.0;
    aEta = 0.00823; bEta = 0.0; cEta = 1.494;
    aPhi = 0.00676; bPhi = 0.0; cPhi = 1.609;
  } else if(2.043<=abs(eta) && abs(eta)<2.172) {
    aEt  = 0.0372; bEt  = 0.985; cEt  = 0.0;
    aEta = 0.0075; bEta = 0.0; cEta = 1.484;
    aPhi = 0.00773; bPhi = 0.0; cPhi = 1.586;
  } else if(2.172<=abs(eta) && abs(eta)<2.322) {
    aEt  = 0.0292; bEt  = 0.967; cEt  = 0.0;
    aEta = 0.00629; bEta = 0.0; cEta = 1.484;
    aPhi = 0.00676; bPhi = 0.0; cPhi = 1.631;
  } else if(2.322<=abs(eta) && abs(eta)<2.500) {
    aEt  = 0.014; bEt  = 0.963; cEt  = 1.24;
    aEta = 0.0; bEta = 0.0; cEta = 1.775;
    aPhi = 0.00652; bPhi = 0.0; cPhi = 1.697;
  } else if(2.500<=abs(eta) && abs(eta)<3.000) {
    aEt  = 0.0653; bEt  = 0.889; cEt  = 2.05;
    aEta = 0.01595; bEta = 0.0; cEta = 2.003;
    aPhi = 0.01746; bPhi = 0.0; cPhi = 1.9;
  } else {
    return false;
  }
  etRes  = et*(sqrt(square(aEt)  + square(bEt/sqrt(et)) + square(cEt/et)));
  etaRes =     sqrt(square(aEta)                        + square(cEta/et));
  phiRes =     sqrt(square(aPhi)                        + square(cPhi/et));
  return true;
}

bool
Resolution::bCaloJetResolution(const double et, const double eta,
                               double& etRes, double& etaRes, double& phiRes) 
{
// Check that eta is in range

  if(abs(eta)>3.0)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

// Set the coefficients according to the eta interval
// If no eta interval qualifies, return false to signal failure.
  if(0.000<=abs(eta) && abs(eta)<0.087) {
    aEt  = 0.0901; bEt  = 1.035; cEt  = 6.2;
    aEta = 0.00516; bEta = 0.0; cEta = 1.683;
    aPhi = 0.0024; bPhi = 0.0; cPhi = 3.159;
  } else if(0.087<=abs(eta) && abs(eta)<0.174) {
    aEt  = 0.0715; bEt  = 1.277; cEt  = 4.77;
    aEta = 0.00438; bEta = 0.0; cEta = 1.72;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.179;
  } else if(0.174<=abs(eta) && abs(eta)<0.261) {
    aEt  = 0.0812; bEt  = 1.192; cEt  = 5.35;
    aEta = 0.00517; bEta = 0.0; cEta = 1.71;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.136;
  } else if(0.261<=abs(eta) && abs(eta)<0.348) {
    aEt  = 0.0713; bEt  = 1.257; cEt  = 4.75;
    aEta = 0.00474; bEta = 0.0; cEta = 1.732;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.166;
  } else if(0.348<=abs(eta) && abs(eta)<0.435) {
    aEt  = 0.0835; bEt  = 1.158; cEt  = 5.08;
    aEta = 0.0047; bEta = 0.0; cEta = 1.744;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.15;
  } else if(0.435<=abs(eta) && abs(eta)<0.522) {
    aEt  = 0.0638; bEt  = 1.298; cEt  = 4.24;
    aEta = 0.00404; bEta = 0.0; cEta = 1.793;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.152;
  } else if(0.522<=abs(eta) && abs(eta)<0.609) {
    aEt  = 0.0676; bEt  = 1.257; cEt  = 4.48;
    aEta = 0.00533; bEta = 0.0; cEta = 1.747;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.112;
  } else if(0.609<=abs(eta) && abs(eta)<0.696) {
    aEt  = 0.0723; bEt  = 1.185; cEt  = 5.28;
    aEta = 0.00511; bEta = 0.0; cEta = 1.745;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.173;
  } else if(0.696<=abs(eta) && abs(eta)<0.783) {
    aEt  = 0.0661; bEt  = 1.292; cEt  = 4.02;
    aEta = 0.00623; bEta = 0.0; cEta = 1.724;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.127;
  } else if(0.783<=abs(eta) && abs(eta)<0.870) {
    aEt  = 0.0773; bEt  = 1.249; cEt  = 4.12;
    aEta = 0.00522; bEta = 0.0; cEta = 1.796;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.123;
  } else if(0.870<=abs(eta) && abs(eta)<0.957) {
    aEt  = 0.082; bEt  = 1.18; cEt  = 5.24;
    aEta = 0.00564; bEta = 0.0; cEta = 1.772;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.125;
  } else if(0.957<=abs(eta) && abs(eta)<1.044) {
    aEt  = 0.0703; bEt  = 1.322; cEt  = 3.81;
    aEta = 0.00337; bEta = 0.0; cEta = 1.832;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.143;
  } else if(1.044<=abs(eta) && abs(eta)<1.131) {
    aEt  = 0.0578; bEt  = 1.39; cEt  = 3.69;
    aEta = 0.00323; bEta = 0.0; cEta = 1.85;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.175;
  } else if(1.131<=abs(eta) && abs(eta)<1.218) {
    aEt  = 0.039; bEt  = 1.508; cEt  = 1.3;
    aEta = 0.00309; bEta = 0.0; cEta = 1.916;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.182;
  } else if(1.218<=abs(eta) && abs(eta)<1.305) {
    aEt  = 0.0722; bEt  = 1.347; cEt  = 4.38;
    aEta = 0.00618; bEta = 0.0; cEta = 1.933;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.223;
  } else if(1.305<=abs(eta) && abs(eta)<1.392) {
    aEt  = 0.0807; bEt  = 1.35; cEt  = 4.38;
    aEta = 0.00889; bEta = 0.0; cEta = 1.961;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.331;
  } else if(1.392<=abs(eta) && abs(eta)<1.479) {
    aEt  = 0.066; bEt  = 1.457; cEt  = 3.54;
    aEta = 0.00747; bEta = 0.0; cEta = 2.079;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.484;
  } else if(1.479<=abs(eta) && abs(eta)<1.566) {
    aEt  = 0.0685; bEt  = 1.42; cEt  = 3.67;
    aEta = 0.01005; bEta = 0.0; cEta = 2.045;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.583;
  } else if(1.566<=abs(eta) && abs(eta)<1.653) {
    aEt  = 0.0; bEt  = 1.561; cEt  = 1.59;
    aEta = 0.0036; bEta = 0.0; cEta = 2.024;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.339;
  } else if(1.653<=abs(eta) && abs(eta)<1.740) {
    aEt  = 0.0736; bEt  = 1.264; cEt  = 4.34;
    aEta = 0.0038; bEta = 0.0; cEta = 2.042;
    aPhi = 0.0; bPhi = 0.0; cPhi = 3.11;
  } else if(1.740<=abs(eta) && abs(eta)<1.830) {
    aEt  = 0.0648; bEt  = 1.234; cEt  = 4.5;
    aEta = 0.0037; bEta = 0.0; cEta = 2.109;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.923;
  } else if(1.830<=abs(eta) && abs(eta)<1.930) {
    aEt  = 0.049; bEt  = 1.243; cEt  = 3.83;
    aEta = 0.0054; aEta = 0.0; aEta = 1.944;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.716;
  } else if(1.930<=abs(eta) && abs(eta)<2.043) {
    aEt  = 0.0661; bEt  = 1.081; cEt  = 4.16;
    aEta = 0.0033; aEta = 0.0; aEta = 1.871;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.548;
  } else if(2.043<=abs(eta) && abs(eta)<2.172) {
    aEt  = 0.0644; bEt  = 1.02; cEt  = 3.89;
    aEta = 0.0; bEta = 0.0; cEta = 1.803;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.365;
  } else if(2.172<=abs(eta) && abs(eta)<2.322) {
    aEt  = 0.0892; bEt  = 0.779; cEt  = 4.28;
    aEta = 0.0; bEta = 0.0; cEta = 1.682;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.148;
  } else if(2.322<=abs(eta) && abs(eta)<2.500) {
    aEt  = 0.0498; bEt  = 0.912; cEt  = 3.53;
    aEta = 0.0; bEta = 0.0; cEta = 1.732;
    aPhi = 0.0; bPhi = 0.0; cPhi = 2.019;
  } else if(2.500<=abs(eta) && abs(eta)<3.000) {
    aEt  = 0.0605; bEt  = 0.861; cEt  = 3.08;
    aEta = 0.0; bEta = 0.0; cEta = 2.032;
    aPhi = 0.0; bPhi = 0.0; cPhi = 1.805;
  } else {
    return false;
  } 
  etRes  = et*(sqrt(square(aEt)  + square(bEt/sqrt(et)) + square(cEt/et)));
  etaRes =     sqrt(square(aEta)                        + square(cEta/et));
  phiRes =     sqrt(square(aPhi)                        + square(cPhi/et));
  return true;
}

