#include "RK4CoupledStep.hpp"

void RK4CoupledStep( std::function <long double (long double, long  double, long double)> f1,
                     std::function <long double (long double, long double, long double)> f2 , long double &t, long double &x1,
                   long double &x2, long double dt){
   long double dx11, dx21, dx31, dx41; long double dx12, dx22, dx32, dx42;

    dx11 = dt*f1(t,x1, x2); dx12 = dt*f2(t,x1, x2);

    dx21 = dt *  f1(t + dt/2, x1 + dx11/2, x2 + dx12/2);
    dx22 = dt *  f2(t + dt/2, x1 + dx11/2, x2 + dx12/2);

    dx31 = dt * f1(t+dt/2, x1 + dx21/2, x2 + dx22/2);
    dx32 = dt * f2(t+dt/2, x1 + dx21/2, x2 + dx22/2);

    dx41 = dt * f1(t + dt, x1 + dx31, x2 + dx32);
    dx42 = dt * f2(t + dt, x1 + dx31, x2 + dx32);

    long double dx1 = (dx11 + 2*(dx21 +dx31) + dx41)/6;
    long double dx2 = (dx12 + 2*(dx22 +dx32) + dx42)/6;

    x1 +=dx1;
    x2 +=dx2;
    t += dt;


}
