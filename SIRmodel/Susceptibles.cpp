#include <iostream>
#include <cmath>

//Declaracion de funciones
double df1(double t, double s, double i, double Beta);
double df2(double t, double s, double i, double Beta);
void UnPasoRungeKutta(double &t0, double &s0,
                      double &i0, double dt, double Beta);

//Constantes modelo SIR
const double Gamma = 0.08;
int main(){
    
    double t,s,i; //Condiciones iniciales
    double dt = 0.5; //Paso Delta t
    for(double beta=0.01; beta <0.5 ; beta+=0.01){
    //RungeKutta
        for(t=0,s=0.999,i=0.001; t<120;){
            UnPasoRungeKutta(t,s,i,dt,beta);
        }
        std::cout << beta/Gamma << "\t" << s << std::endl;
    }
        return 0;
}

//Implementacion

double df1(double t, double s, double i, double Beta){
    
    return -Beta*s*i;
}
double df2(double t, double s, double i, double Beta){
    return Beta*s*i - Gamma*i;
}

void UnPasoRungeKutta(double &t0, double &s0,
                      double &i0, double dt, double Beta){
    
    double ds1, ds2, ds3, ds4;
    double di1, di2, di3, di4;
    //Primer delta
    ds1 = dt*df1(t0,s0,i0, Beta);
    di2 = dt*df2(t0,s0,i0, Beta);
    //Segundo delta
    ds2 = dt*df1(t0+dt/2,s0+ds1/2,i0+di1/2,Beta);
    di2 = dt*df2(t0+dt/2,s0+ds1/2,i0+di1/2,Beta);
    //Tercer delta
    ds3 = dt*df1(t0+dt/2,s0+ds2/2,i0+di2/2,Beta);
    di3 = dt*df2(t0+dt/2,s0+ds2/2,i0+di2/2,Beta);
    //Cuarto delta
    ds4 = dt*df1(t0+dt,s0+ds3,i0+di3,Beta);
    di4 = dt*df2(t0+dt,s0+ds3,i0+di3,Beta);
    //Deltas totales
    s0+=(ds1+2*(ds2+ds3)+ds4)/6;
    i0+=(di1+2*(di2+di3)+di4)/6;
    
    t0+=dt;
}
