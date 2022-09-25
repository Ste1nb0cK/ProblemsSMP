#include<iostream>
#include<cmath>
double f1(double t, double x1, double x2); 
double f2(double t, double x1, double x2, double lambda);
void Runge_one_step(double &t0, double &x10, double &x20, double dt, double lambda);

int main(void){

  double dt = 0.01;
  double x1 = 1;
  double x2 = 0;

  for(double t=0.01; t<10; ){
    Runge_one_step(t, x1, x2, dt, 14.9323);
    std::cout<<t<<"\t"<<x1<<std::endl;
  }

  
  return 0;
}


double f2(double t, double x1, double x2, double lambda){

  
  
  return -lambda*lambda*x1-(x2/t);
}

double f1(double t, double x1, double x2){

  return x2;
}

void Runge_one_step(double &t0, double &x10, double &x20, double dt, double lambda){

  double dx11, dx21, dx31, dx41;
  double dx12, dx22, dx32, dx42;


  dx11 = f1(t0,x10,x20)*dt;                                  dx12 = f2(t0,x10,x20, lambda)*dt;

  dx21 = dt*f1(t0+ 0.5*dt, x10 + 0.5*dx11, x20 + 0.5*dx12);  dx22 = dt*f2(t0+ 0.5*dt, x10 + 0.5*dx11, x20 + 0.5*dx12, lambda);

  dx31 = dt*f1(t0+ 0.5*dt, x10 + 0.5*dx21, x20 + 0.5*dx22);  dx32 = dt*f2(t0+ 0.5*dt, x10 + 0.5*dx21, x20 + 0.5*dx22, lambda);

  dx41 = dt*f1(t0+ 0.5*dt, x10 + dx31, x20 + dx32);          dx42 = dt*f2(t0+ 0.5*dt, x10 + dx31, x20 + dx32, lambda);  

  
  x10+=(dx11 + 2*dx21 + 2*dx31 + dx41)/6;
  x20+=(dx12 + 2*dx22 + 2*dx32 + dx42)/6;
  t0+=dt;

}


