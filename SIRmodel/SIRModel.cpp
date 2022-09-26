#include <iostream>
#include <cmath>

//Declaracion de funciones
double df1(double t, double s, double i);
double df2(double t, double s, double i);
void UnPasoRungeKutta(double &t0, double &s0, double &i0, double dt);

//Constantes modelo SIR
const double Beta = 0.35;
const double Gamma = 0.08;
int main(){
    
    double t,s,i; //Condiciones iniciales
    double dt = 1; //Paso Delta t
    
    for(t=0,s=0.999,i=0.001; t<120; ){
        std::cout << t << "\t" << s << "\t" << i  << "\t" << 1-s-i << std::endl;
        UnPasoRungeKutta(t,s,i,dt);
    }
    return 0;
}

//Implementacion

double df1(double t, double s, double i){
    
    return -Beta*s*i;
}
double df2(double t, double s, double i){
    return Beta*s*i - Gamma*i;
}

void UnPasoRungeKutta(double &t0, double &s0,
                      double &i0, double dt){
    
    double ds1, ds2, ds3, ds4;
    double di1, di2, di3, di4;
    //Primer delta
    ds1 = dt*df1(t0,s0,i0);
    di2 = dt*df2(t0,s0,i0);
    //Segundo delta
    ds2 = dt*df1(t0+dt/2,s0+ds1/2,i0+di1/2);
    di2 = dt*df2(t0+dt/2,s0+ds1/2,i0+di1/2);
    //Tercer delta
    ds3 = dt*df1(t0+dt/2,s0+ds2/2,i0+di2/2);
    di3 = dt*df2(t0+dt/2,s0+ds2/2,i0+di2/2);
    //Cuarto delta
    ds4 = dt*df1(t0+dt,s0+ds3,i0+di3);
    di4 = dt*df2(t0+dt,s0+ds3,i0+di3);
    //Deltas totales
    s0+=(ds1+2*(ds2+ds3)+ds4)/6;
    i0+=(di1+2*(di2+di3)+di4)/6;
    
    t0+=dt;
}
