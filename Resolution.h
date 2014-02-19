#ifndef Resolution_h
#define Resolution_h

class Resolution {

public:

  Resolution();

  virtual ~Resolution();

  bool electronResolution(const double et, const double eta,
                          double& etRes, double& etaRes, double& phiRes);

  bool muonResolution(const double et, const double eta,
                      double& etRes, double& etaRes, double& phiRes);

  bool caloMETResolution(const double et, double& etRes, double& etaRes, double& phiRes);

  bool PFMETResolution(const double et, double& etRes, double& etaRes, double& phiRes);

  bool udscCaloJetResolution(const double et, const double eta,
                             double& etRes, double& etaRes, double& phiRes);

  bool udscPFJetResolution(const double et, const double eta,
                           double& etRes, double& etaRes, double& phiRes);

  bool bCaloJetResolution(const double et, const double eta,
                          double& etRes, double& etaRes, double& phiRes);

  bool bPFJetResolution(const double et, const double eta,
                        double& etRes, double& etaRes, double& phiRes);

  inline double square(const double x) {return x*x;}

};  

#endif
