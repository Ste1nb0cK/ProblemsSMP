#include<iostream>
#include<cmath>
#include<vector>
double f1(double t, double x1, double x2);
double f2(double t, double x1, double x2, double lambda);
void Runge_one_step(double &t0, double &x10, double &x20, double dt, double lambda);
double RGK_bisection(double a, double b, double x01, double x02);
double f3(double lambda, double x10, double x20);

int main(){

  std::cout<<"Los valores de lambda que para los que se tienen modos normales del oscilación son "<<std::endl;
  std::cout<<RGK_bisection(2, 5, 0.17, -0.56)<<std::endl;
  std::cout<<RGK_bisection( 5, 7, 0.007,-0.52)<<std::endl;
  std::cout<<RGK_bisection( 8 ,9, 0.15, -0.24)<<std::endl;
  std::cout<<RGK_bisection(10, 12, 0.01, -0.27)<<std::endl;
  std::cout<<RGK_bisection(14 ,15, 0.04, -0.2)<<std::endl;
    

  
  return 0;
}

double f3(double lambda, double x10, double x20){

  double dt = 0.001;
  double t;
  double x1, x2;
  for(t=0.001, x1=x10, x2=x20; t<=1; ){

    Runge_one_step(t,x1,x2,dt,lambda);
  }

  return x1;

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

double RGK_bisection(double a, double b, double x01, double x02){
  double x, m, fa, fm;
  double Errmax= 1e-7;

  fa=f3(a,x01,x02);

  while(b-a>Errmax){
    m = (b+a)/2; fm = f3(m, 1, 0);
    if(fa*fm > 0){
      a = m;
      fa = fm;
      
    }
    else{ b = m;}

   }
  return (a+b)/2;
}
