#include<iostream>
#include<cmath>
#include<fstream>

double f(double x);
double Simpson(double a, double b, int n);
double Bessel(double alpha, double x);
double Bisection(double a, double b, double alpha);


int main(int argc, char *argv[]){
  double alpha = 0; double x = 0;;
  double lambda = 0;
  //El usuario debe especificar a la hora de ejecutar para qué modo normal desea graficar obtener los datos de R(r).
  if (std::atof(argv[1])==1) lambda = 2.43797;
  if (std::atof(argv[1])==2) lambda = 5.4858;
  if (std::atof(argv[1])==3) lambda = 8.66974;
  if (std::atof(argv[1])==4) lambda = 11.7982;
  if (std::atof(argv[1])==5) lambda = 14.9112;


  for (x=0; x<=10; x+=0.01){

     std::cout<<x<<"\t"<< Bessel(alpha, lambda*x)<<std::endl;

     }

  //Descomentar si se desea acceder a los ceros teóricos calculados con la bisección.

  /* std::cout<<"Los ceros de la solución teórica están en"<<std::endl;
  std::cout<<Bisection(2, 3, alpha)<<std::endl;
  std::cout<<Bisection(5, 6, alpha)<<std::endl;
  std::cout<<Bisection(6, 9, alpha)<<std::endl;
  std::cout<<Bisection(10, 12, alpha)<<std::endl;
  std::cout<<Bisection(13, 16, alpha)<<std::endl;*/
   
  
  

  return 0; 
}
double f(double t, double x, double alpha){

  return cos(alpha*t-x*sin(t));
}

double Simpson(double a, double b, int n, double x, double alpha){
  double t = 0;
  double h = 0;
  double sum = 0;
  n*=2; h=(b-a)/n;
  for(int ii = 0; ii<=n; ii++){
    t = a + ii*h;
    if(ii==0 || ii == n){
      sum += f(t, x, alpha);
    }
    else if(ii%2==0){
      sum += 2*f(t, x, alpha);
      }
    else {
      sum += 4*f(t, x, alpha); 
    }
    
  
  }
   return sum*(h/3);

}

double Bessel(double alpha, double x){

  return (1/(M_PI))*Simpson(0, M_PI, 1000, x, alpha);

}

double Bisection(double a, double b, double alpha){
  double m, Bessela, Besselm;
  double Errmax= 1e-7;
  Bessela=Bessel(alpha, a);
  while(b-a>Errmax){
    m = (b+a)/2; Besselm = Bessel(alpha, m);
    if(Bessela*Besselm > 0){
      a = m;
      Bessela = Besselm;
      
    }
    else{ b = m;}

   }
  return (a+b)/2;
}

