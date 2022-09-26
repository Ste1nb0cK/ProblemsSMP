// Simular el movimiento de N particulas Lennard Jones
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double Epsilon = 1.0;
const double r_0=10.0; //Distancia de equilibrio
const double K=1.0e4, Gamma=50, Kcundall=10, MU=0.4;
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau; double I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
	      double theta0,double omega0);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void AdicioneTorque(double tau0){tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double Gettheta(void){return theta;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
		    double theta0,double omega0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);  omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
    void CalculeFuerzas(Cuerpo * Grano);
    void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
    void FuerzaParedes(Cuerpo & Grano);
};
void Colisionador::FuerzaParedes(Cuerpo & Grano){
    //Pared izquierda
    double si;
        si = Grano.Getx()-Grano.R;
        if(si<0.0){
            vector3D Fi;
            Fi.load(K*pow(-si,1.5),0.0,0.0);
            Grano.AdicioneFuerza(Fi);
    }
            //Pared abajo
    double sd;
        sd = Grano.Gety()-Grano.R;
        if(sd<0.0){
            vector3D Fd;
            Fd.load(0.0,K*pow(-sd,1.5),0.0);
            Grano.AdicioneFuerza(Fd);
    }
        //Pared Arriba
        double su;
        su = -Grano.Gety()+Ly-Grano.R;
        if(su<0.0){
            vector3D Fu;
            Fu.load(0.0,-K*pow(-su,1.5),0.0);
            Grano.AdicioneFuerza(Fu);
    }
        //Pared derecha
            double sr;
            sr = -Grano.Getx()+Lx-Grano.R;
            if(sr<0.0){
                vector3D Fr;
                Fr.load(-K*pow(-sr,1.5),0.0,0.0);
                Grano.AdicioneFuerza(Fr);
    }
}
void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++)
    Grano[i].BorreFuerza();
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++){
      for(j=i+1;j<N;j++) CalculeFuerzaEntre(Grano[i], Grano[j]);
      FuerzaParedes(Grano[i]);
  }
  
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21=Grano2.r-Grano1.r;
  double d = r21.norm();
  double s = Grano1.R+Grano2.R-d;
    vector3D n=r21*(1.0/d);
    double pauli=pow(r_0/d,12.0);
    double waals=pow(r_0/d,6.0);
    double aux = 12.0*Epsilon/d*(pauli-waals);
    vector3D F2=n*aux;
    Grano2.AdicioneFuerza(F2);   Grano1.AdicioneFuerza(F2*(-1));   
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
    // cout<<"set terminal gif animate"<<endl; 
    //cout<<"set output 'Lennard.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[N];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5, kT=10.0, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,tmax=200,tcuadro=tmax/1000,dt=0.001;
  double dx=10.0, dy=10.0;
  double Theta;
  
  InicieAnimacion(); //Dibujar
  //Inicializar las moléculas
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
      Grano[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy,V0*cos(Theta),V0*sin(Theta), m0,R0,0,1);//OJO
    }
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }

    //--- Muevase por PEFRL ---
    
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);  

  }   

  
  return 0;
}

  
