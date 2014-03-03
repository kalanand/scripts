
#include "BToDKTo3piK/ErrMatStatGauss.hh"

ClassImp(ErrMatStatGauss)

ErrMatStatGauss::ErrMatStatGauss() {
  for (int c = 0; c < NCHARGES; ++c) {
    _rhoErrNSpreads[c] = 0.;
    _thetaErrNSpreads[c] = 0;
  }
}

ErrMatStatGauss::ErrMatStatGauss(const Poly3 & rhoErr,
				 const Poly3 & thetaErr,
				 const Poly3 & rhoErrSpread,
				 const Poly3 & thetaErrSpread,
				 const Poly3 & rhoPullSigma,
				 const Poly3 & thetaPullSigma,
				 double rhoErrNSpreads[NCHARGES],
				 double thetaErrNSpreads[NCHARGES]) {
init(rhoErr,
     thetaErr,
     rhoErrSpread,
     thetaErrSpread,
     rhoPullSigma,
     thetaPullSigma,
     rhoErrNSpreads, 
     thetaErrNSpreads);
}

void ErrMatStatGauss::init(const Poly3 & rhoErr,
			   const Poly3 & thetaErr,
			   const Poly3 & rhoErrSpread,
			   const Poly3 & thetaErrSpread,
			   const Poly3 & rhoPullSigma,
			   const Poly3 & thetaPullSigma,
			   double rhoErrNSpreads[NCHARGES],
			   double thetaErrNSpreads[NCHARGES]) {
  _rhoErr = rhoErr;
  _thetaErr = thetaErr;
  _rhoErrSpread = rhoErrSpread;
  _thetaErrSpread = thetaErrSpread;
  _rhoPullSigma = rhoPullSigma;
  _thetaPullSigma = thetaPullSigma;

  for (int c = 0; c < NCHARGES; ++c) {
    _rhoErrNSpreads[c] = rhoErrNSpreads[c];
    _thetaErrNSpreads[c] = thetaErrNSpreads[c];
  }
}

ErrMatStatGauss::~ErrMatStatGauss() {}
  
TMatrixD ErrMatStatGauss::matrix(double rM, double tM, double rP, double tP) 
  const {
  // The algorithm is as follows: For a given true rho, find the average
  // rho and theta errors and the spreads in the rho and theta error 
  // distributions. The final error is the average error plus the spread
  // times nSpreads:

  TMatrixD result(4,4);

  double errRhoN = _rhoPullSigma.getVal(rM) * measuredErrRho(rM,MINUS);
  double thetaRhoN = _thetaPullSigma.getVal(rM) * measuredErrTheta(rM,MINUS);  
  double errRhoP = _rhoPullSigma.getVal(rP) * measuredErrRho(rP,PLUS);
  double thetaRhoP = _thetaPullSigma.getVal(rP) * measuredErrTheta(rP,PLUS);

  result[0][0] = errRhoN * errRhoN;
  result[1][1] = thetaRhoN * thetaRhoN;
  result[2][2] = errRhoP * errRhoP;
  result[3][3] = thetaRhoP * thetaRhoP;

  return result;
}


double ErrMatStatGauss::measuredErrRho(double rho, Charge c) const {
  return 
    _rhoErr.getVal(rho) + _rhoErrSpread.getVal(rho) *_rhoErrNSpreads[c];
}


double ErrMatStatGauss::measuredErrTheta(double rho, Charge c) const {
  return 
    _thetaErr.getVal(rho) + _thetaErrSpread.getVal(rho) *_thetaErrNSpreads[c];
}
