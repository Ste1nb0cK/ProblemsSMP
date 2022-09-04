#include "EX1-funciones.h"

int main(int argc, char **argv){
  double beta = std::atof(argv[1]);
  double gamma = std::atof(argv[2]);
  double t = std::atof(argv[3]);
  
  double S0 = 0.999; double I0 = 0.001;   double R0 = S0 + I0 - 1;
  vector3D inicial; inicial.load(S0,I0,R0);
  
  double dt = 1e-2;
  bool print = 1;
  
  vector3D R = RungeKutaIterado(t, inicial,dt,  beta, gamma, print);
  return 0;
}

 


