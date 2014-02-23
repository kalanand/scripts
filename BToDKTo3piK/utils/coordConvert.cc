// Handles coordinate conversions

#include <math>


double xFromCP(double r, double phase) {
  return r * cos(phase);
}


double yFromCP(double r, double phase) {
  return r * sin(phase);
}


double rhoFromCart(double x, double y) {
  double X0 = dalitzHolderN.sigGoodD0Type().x0()->getVal();
  return sqrt((x - X0)*(x - X0) + y * y);
}


double thetaFromCart(double x, double y) {
  double X0 = dalitzHolderN.sigGoodD0Type().x0()->getVal();
  double result = atan2(y, (x - X0)) / BdkDDalitzAmp::DEGTORAD;
  if (result < 0) result += 360;
  if (result > 360) result -= 360;
  return result;
}


double rhoFromCP(double r, double phase) {
  return rhoFromCart(xFromCP(r, phase), yFromCP(r, phase));
}


double thetaFromCP(double r, double phase) {
  return thetaFromCart(xFromCP(r, phase), yFromCP(r, phase));
}
