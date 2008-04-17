
#include <iostream>
using namespace std;
#include "BToDKTo3piK/ErrMatSyst.hh"

ClassImp(ErrMatSyst)

ErrMatSyst::ErrMatSyst(Bool_t includeModel, Bool_t modelOnly) :
  _includeModel(includeModel),
  _modelOnly(modelOnly),
  // For the rho Err objects, the total error is from tables 43+44 in the bad:
  _errRm(sqrt(square(0.0462) + square(0.0375)), 
	 // The proportional error is from table 44 plus the first line
	 // of table 43:
	 sqrt(square(0.0230)+square(0.0375)),
	 // The model error:
	 0.0261),
  // Same for rho+:
  _errRp(sqrt(square(0.0464) + square(0.0360)), 
	 sqrt(square(0.0203)+square(0.0360)),
	 0.0300),
  // For the theta Err objects, use only the total error and model
  // error from table 43:
  _errTm(19.06, 0, 15.42),
  _errTp(12.69, 0, 12.16)
{
}
    
TMatrixD ErrMatSyst::matrix(double rM, double tM, double rP, double tP) 
  const {
  // start with the constant errors, those that don't directly affect rho.
  // they are the ones that don't have an entry in the A_DKsig or N_DKsig 
  // columns of the systematics table (tab 43) in the BAD.

  TMatrixD result(4,4);
  
  if (kTRUE == _modelOnly) {
    // Only model errors:
    result[0][0] = square(_errRm._errModel);
    result[1][1] = square(_errTm._errModel);
    result[2][2] = square(_errRp._errModel);
    result[3][3] = square(_errTp._errModel);

    return result;
  }

  // If we are here, then doing non-model errors and (if includeModel==kTRUE)
  // model errors:
  
  int subtractModel = 0;
  if (kFALSE == _includeModel) {
    subtractModel = 1;
  }
  
  result[0][0] = square(_errRm._errConst) - subtractModel * square(_errRm._errModel);
  result[1][1] = square(_errTm._errConst) - subtractModel * square(_errTm._errModel);
  result[2][2] = square(_errRp._errConst) - subtractModel * square(_errRp._errModel);
  result[3][3] = square(_errTp._errConst) - subtractModel * square(_errTp._errModel);

  /* Add the part of the rho error that is proportional to the yield.
     this error equals to the proportional error measured on the data,
     times the ratio 
     (1 + rho^2 - x_0^2)/(2 rho) 
     -------------------------------------
     (1 + rho_data^2 - x_0^2)/(2 rho_data) 
     
     where rho in the numerator is the true rho used to generate the point,
     and rho_data in the denominator is the value of rho found on the data.
   */

  const double x0Sq = square(0.850);
  const double rMData = 0.715;
  const double rPData = 0.748;

  double numerM = 1 + square(rM) - x0Sq;
  double denomM = 1 + square(rMData) - x0Sq;
  double propErrM = numerM / denomM * _errRm._errProp;
  
  double numerP = 1 + square(rP) - x0Sq;
  double denomP = 1 + square(rPData) - x0Sq;
  double propErrP = numerP / denomP * _errRp._errProp;

  // add this to the diagonal elements:
  result[0][0] += square(propErrM);
  result[2][2] += square(propErrP);

  // and the off-diagonal elements with 100% correlation:
  result[0][2] += propErrM * propErrP;
  result[2][0] += propErrM * propErrP;
  
  return result;
}


void ErrMatSyst::testClass() {
  ErrMatSyst wModel;
  cout << "Error matrix with model:" << endl;
  wModel.matrix(0.715, 172.7, 0.748, 146.9).Print();

  ErrMatSyst woModel(kFALSE);
  cout << "Error matrix without model:" << endl;
  woModel.matrix(0.715, 172.7, 0.748, 146.9).Print();
}
