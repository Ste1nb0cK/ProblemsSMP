// Simular el movimiento de N particulas Lennard Jones
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double Epsilon = 1.0;
const double r_0=10.0; //Distancia de equilibrio
const double K=1.0e4;
const double Lx=60, Ly=60;
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
  vector3D r,V,F,V_A; double m,R, theta,P; int idy,idx;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,double theta0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Mueva_P(void);
  void Sumeinteraccionx();
  void Sumeinteracciony();
  void CambieVAntigua(vector3D VV){V_A=VV;};
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inlin
  double GetP(void){return P;};
  double Gettheta(void){return theta;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,double theta0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0;
  P=0; idx=0;idy=0; V_A.load(0.0,0.0,0.0);
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);
}
void Cuerpo::Mueva_P(void){ 
  if(idx>0){
    P+=abs(V.x()-V_A.x());
    idx=0;
  }
  if(idy>0){
    P+=abs(V.y()-V_A.y());
    idy=0;
  }
}
void Cuerpo::Sumeinteraccionx(void){
  idx+=1;
}
void Cuerpo::Sumeinteracciony(void){
  idy+=1;
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
  
  vector3D V_ant,Fi,Fd,Fu,Fr;

  double si;
  si = Grano.Getx()-Grano.R;
  if(si<0.0){
    Grano.Sumeinteraccionx();       
      if (Grano.idx==1){
	V_ant=Grano.V;
	Grano.CambieVAntigua(V_ant);
      }
      
      Fi.load(K*pow(-si,1.5),0.0,0.0);
      Grano.AdicioneFuerza(Fi);
      
  }
  //Pared abajo
  double sd;
  sd = Grano.Gety()-Grano.R;
  if(sd<0.0){
    Grano.Sumeinteracciony();       
      if (Grano.idx==1){
	V_ant=Grano.V;
	Grano.CambieVAntigua(V_ant);
      }
      
      Fd.load(0.0,K*pow(-sd,1.5),0.0);
      Grano.AdicioneFuerza(Fd);
      
  }
  //Pared Arriba
  double su;
  su = -Grano.Gety()+Ly-Grano.R;
  if(su<0.0){
    
    Grano.Sumeinteracciony();       
    if (Grano.idx==1){
      V_ant=Grano.V;
      Grano.CambieVAntigua(V_ant);
    }
    
      
      Fu.load(0.0,-K*pow(-su,1.5),0.0);
      Grano.AdicioneFuerza(Fu);
  }
  //Pared derecha
  double sr;
  sr = -Grano.Getx()+Lx-Grano.R;
  if(sr<0.0){
    Grano.Sumeinteraccionx();       
    if (Grano.idx==1){
      V_ant=Grano.V;
	Grano.CambieVAntigua(V_ant);
    }
    
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
    FuerzaParedes(Grano[i]);//orden
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
  double m0=1.0, R0=2.5;
  int i,ix,iy;
  double t,tdibujo,tmax=260,tcuadro=tmax/1000,dt=0.001;
  double dx=10.0, dy=10.0;
  double Theta;

  for (double kT=2.0; kT<21 ;kT++){ 
    //InicieAnimacion(); //Dibujar
    double V0=sqrt(2*kT/m0);
    //Inicializar las moléculas
    for(ix=0;ix<Nx;ix++)
      for(iy=0;iy<Ny;iy++){
	Theta=2*M_PI*ran64.r();
	//--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
	Grano[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy,V0*cos(Theta),V0*sin(Theta), m0,R0,0);//OJO
      }
    for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    /*
      if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    */
      
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
      for(i=0;i<N;i++){
	Grano[i].Mueva_r(dt,epsilon);
	if ( t>60 && t<260){
	  Grano[i].Mueva_P();
	}
      }
    }
    double Presion=0;
    for(int j=0;j<N;j++)Presion+=Grano[j].GetP();
    
    cout<<kT<<" " <<Presion/(200*Ly)<<endl;
  }
  return 0;
}

  
