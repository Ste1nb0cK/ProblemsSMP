#ifndef _ITERATE_H_
#define _ITERATE_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "vector.h"

vector3D derivada(double t, vector3D R, double beta, double gamma);
void UnPasoRungeKutta(double & t, vector3D & R, double dt, double beta, double gamma);
vector3D Asintota(vector3D R0, double dt, double beta, double gamma);
vector3D RungeKutaIterado(double  t, vector3D  R, double dt, double beta, double gamma,bool print);


vector3D derivada(double t, vector3D R, double beta, double gamma){
  vector3D derivada;
  derivada.load(-beta*R.x()*R.y(), beta*R.x()*R.y() - gamma*R.y(), gamma*R.y());
  return derivada;
}


// las variables a actualizar se pasan por referencia
void UnPasoRungeKutta(double & t, vector3D & R, double dt, double beta, double gamma){
  vector3D dR, dR1, dR2, dR3, dR4;
  
  dR1 = dt*derivada(t, R, beta, gamma);
  dR2 = dt*derivada(t + dt/2, R + dR1/2, beta, gamma);
  dR3 = dt*derivada(t + dt/2, R + dR2/2, beta, gamma);
  dR4 = dt*derivada(t + dt, R + dR3, beta, gamma);
  

  dR = (dR1 + 2*dR2 + 2*dR3 + dR4)/6;
  t += dt;
  R += dR;
}

vector3D RungeKutaIterado(double  t, vector3D  R0, double dt, double beta, double gamma,bool print){
  int N = std::ceil(t/dt);
  dt = t/N;
  double t0 = 0;
  
  for( int ii = 0; ii < N ;  ){
    if( print == 1)
      std::cout<< t0 << "\t" << R0.x()  << "\t" << R0.y() << "\t" << R0.z() << std::endl;
    UnPasoRungeKutta(t0,R0,dt,beta, gamma);
    ++ii;
  }
  return R0;
}

#endif
