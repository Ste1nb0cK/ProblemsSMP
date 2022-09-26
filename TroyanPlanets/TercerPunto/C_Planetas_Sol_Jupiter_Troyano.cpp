#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes globales

const int N=3;
const double G=1.0;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------------- Clase Cuerpo --------------
class Cuerpo{
private:
  vector3D r,V,F;  double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};
  void SumeFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  void PrintOneRotado(Cuerpo & Referencia);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(dt*coeficiente);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(dt*coeficiente/m);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//--------------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Planeta);
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);    
};
void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  //Borrar fuerzas
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D r21,n,F1; double d21,F;
  r21=Planeta2.r-Planeta1.r; d21=r21.norm(); n=r21/d21;
  F=G*Planeta1.m*Planeta2.m*std::pow(d21,-2.0);
  F1=F*n; Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza(F1*(-1));
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<endl; 
  std::cout<<"set output 'DosPlanetas.gif'"<<endl;
  std::cout<<"unset key"<<endl;
  std::cout<<"set xrange[-12:12]"<<endl;
  std::cout<<"set yrange[-12:12]"<<endl;
  std::cout<<"set size ratio -1"<<endl;
  std::cout<<"set parametric"<<endl;
  std::cout<<"set trange [0:7]"<<endl;
  std::cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    std::cout<<endl;
}
/// --------------- Rotación del nuevo eje ----------------
void PrintAllRotado( Cuerpo * Planeta,Cuerpo & Jup, double t){
  std::cout<< t;
  for (int i = 0; i < N ; i++) Planeta[i].PrintOneRotado(Jup);
  std::cout <<std::endl;
}
void Cuerpo::PrintOneRotado(Cuerpo & Referencia){
  vector3D Proy = (r*Referencia.r/(Referencia.r.norm2()))*Referencia.r;
  double x = Proy.norm();
  double y = (r-Proy).norm();
  std::cout<<"\t"<<x<<"\t"<<y;
}

 
int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m0 =1047.346599, m1=1.0, r=1000.0;
  //  double m0  = 1.0e5, m1=1.0, r=1000.0;
  double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  double omega= std::sqrt(G*M/(r*r*r)), T=2*M_PI/omega, V0=omega*x0, V1=omega*x1;
  double t,tmax=80.1*T,dt=25.1;
  double tdibujo,tcuadro=T/100.0;
  int i;
  // T es mas o menos 7000
  double x2 = r*std::cos(M_PI/3); double y2 = r*std::sin(M_PI/3);
  double Vx2 = V1*(-std::sin(M_PI/3)); double Vy2 = V1*std::cos(M_PI/3);
  double  m2 = 5e-3;
    //---------------(x0,y0,z0,Vx0,Vy0,Vz,m0,R0)
  Planeta[0].Inicie(x0, 0.0, 0.0,  0.0, V0, 0, m0, 1.0); // Sol
  Planeta[1].Inicie(x1, 0.0, 0.0,  0, V1, 0.0, m1, 0.5); // Jupiter
  Planeta[2].Inicie(x2, y2, 0.0, Vx2, Vy2, 0.0, m2, 0.1);// troyanos
  //  InicieAnimacion();
  
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){

    if(tdibujo > tcuadro){
      PrintAllRotado( Planeta, Planeta[1], t);
      
      /*
      std::cout<<Planeta[0].Getx()<<"\t"<<Planeta[0].Gety()<<
      "\t"<<Planeta[1].Getx()<<"\t"<<Planeta[1].Gety()<<endl;
      /*
      InicieCuadro();
      for(i=0;i<N;i++) Planeta[i].Dibujese();
      TermineCuadro();
      */
      tdibujo=0;
    }
        


    // Mover por PEFRL
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zeta);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zeta);   
  }
  
return 0;
}
