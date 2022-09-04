#include "EX1-funciones.h"

int main(int argc, char **argv){
  double gamma = std::atof(argv[1]);
  double t = std::atof(argv[2]);
  
  double S0 = 0.999;   double I0 = 0.001;   double R0 = S0 + I0 - 1;
  vector3D R; R.load(S0, I0, R0);
  
  
  double dt = 1e-2;
  bool print = 0;

  double deltaRatio = 1e-2;

  for(double Ratio = 0; Ratio < 5.0; Ratio += deltaRatio)
    std::cout<<Ratio<< "\t"<< RungeKutaIterado(t,R,dt,Ratio*gamma,gamma, print).x() << "\n";
 
  return 0;
}
