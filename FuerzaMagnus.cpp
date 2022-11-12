#include <iostream>
#include <fstream>
#include <cmath>
//#include <math.h>
using namespace std;

const int Lx=512;
const int Ly=64;

const int Q=9;

//--------------------- Clase LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  double w[Q];      //Weights 
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Start(double rho0,double Ux0,double Uy0);
  void Collision(double tau);
  void ImposeFields(double Ufan,double omega,double R, int ixc, int iyc);
  void Advection(void);
  void Print(const char * NameFile,double Ufan);
  void Sigma(double ix, double iy, double tau,double & sigmaxx, double & sigmaxy, double & sigmayy);
  void dF_drag(double x, double y, double dAx, double dAy,double tau, double & dFx, double & dFy);
  void FuerzaTotalDrag(int NumPoints, double *Points, double *VecArea,double tau,double & Fx, double & Fy);
};
LatticeBoltzmann::LatticeBoltzmann(void){
   //Set the weights
  w[0]=4.0/9;  w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

            Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
            Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
 //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}  
double  LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}  
void LatticeBoltzmann::Start(double rho0,double Ux0,double Uy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i);
	f[n0]=feq(rho0,Ux0,Uy0,i);
      }
}  
void LatticeBoltzmann::Collision(double tau){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  double Utau=1.0/tau;
  double UmUtau=1-Utau;

  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }  
}
void LatticeBoltzmann::ImposeFields(double Ufan,double omega,double R,int ixc, int iyc){
  int i,ix,iy,n0; double rho0;
   double R2=R*R;
  //go through all cells, looking if they are fan or obstacle
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //fan
      if(ix==0)
	for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
      //obstacle
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2) 
	for(i=0;i<Q;i++) {n0=n(ix,iy,i); fnew[n0]=feq(rho0,-omega*(iy-iyc),omega*(ix-ixc),i);}
      //An extra point at one side to break the symmetry
      //else if(ix==ixc && iy==iyc+R+1)
      //for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}	
    }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(const char * NameFile,double Ufan){
  ofstream MyFile(NameFile); double rho0,Ux0,Uy0; int ix,iy;
  for(ix=0;ix<Lx;ix+=4){
    for(iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

void LatticeBoltzmann::Sigma(double ix, double iy,double tau, double & sigmaxx, double & sigmaxy, double & sigmayy){
  double eta = (tau - 0.5)/3.0;
  double Ux[Q]{0};
  double Uy[Q]{0};
  double rho0[Q]{0};
  double dUxdX = 0, dUydY = 0, dUxdY = 0, dUydX = 0;
  for(int q = 0; q < Q; q++){
    rho0[q] = rho(ix + Vx[q],iy + Vy[q],true);
    Ux[q] = Jx(ix + Vx[q],iy + Vy[q],true)/rho0[q];
    Uy[q] = Jy(ix + Vx[q],iy + Vy[q],true)/rho0[q];
  }
  
  for(int q = 0; q < Q; q++){
    dUxdX += w[q]*Vx[q]*Ux[q];
  }
  dUxdX *= 3;

  for(int q = 0; q < Q; q++){
    dUydX += w[q]*Vx[q]*Uy[q];
  }
  dUydX *= 3;
  
  for(int q = 0; q < Q; q++){
    dUydY += w[q]*Vy[q]*Uy[q];
  }
  dUydY *= 3;

  for(int q = 0; q < Q; q++){
    dUydY += w[q]*Vy[q]*Ux[q];
  }
  dUxdY *= 3;

  double p = rho0[0]/3;
  sigmaxx = -p + eta*2*dUxdX;
  sigmayy = -p + eta*2*dUydY;
  sigmaxy = eta*(dUxdY + dUydY);
}

void LatticeBoltzmann::dF_drag(double x, double y, double dAx, double dAy,double tau, double & dFx, double & dFy){
  int ix = std::floor(x);
  int iy = std::floor(y);
  
  double AUXsigmaxx =0, sigmaxx = 0;
  double AUXsigmaxy = 0, sigmaxy = 0;
  double AUXsigmayy = 0, sigmayy = 0;
  
  double u = x - ix;
  double v = y - iy;
  
  Sigma(ix, iy,tau, AUXsigmaxx, AUXsigmaxy, AUXsigmayy);
  sigmaxx += (1-u)*(1-v)*AUXsigmaxx;
  sigmaxy += (1-u)*(1-v)*AUXsigmaxy;
  sigmayy += (1-u)*(1-v)*AUXsigmayy;

  Sigma(ix + 1, iy,tau, AUXsigmaxx, AUXsigmaxy, AUXsigmayy);
  sigmaxx += u*(1-v)*AUXsigmaxx;
  sigmaxy += u*(1-v)*AUXsigmaxy;
  sigmayy += u*(1-v)*AUXsigmayy;

  Sigma(ix, iy + 1, tau,AUXsigmaxx, AUXsigmaxy, AUXsigmayy);
  sigmaxx += v*(1-u)*AUXsigmaxx;
  sigmaxy += v*(1-u)*AUXsigmaxy;
  sigmayy += v*(1-u)*AUXsigmayy;
  
  Sigma(ix + 1, iy + 1, tau,AUXsigmaxx, AUXsigmaxy, AUXsigmayy);
  sigmaxx += v*u*AUXsigmaxx;
  sigmaxy += v*u*AUXsigmaxy;
  sigmayy += v*u*AUXsigmayy;

  dFx = sigmaxx*dAx + sigmaxy*dAy;
  dFy = sigmaxy*dAx + sigmayy*dAy;
}

void LatticeBoltzmann::FuerzaTotalDrag(int NumPoints, double *Points, double *VecArea,double tau,double & Fx, double & Fy){
  double AUXFx = 0; double AUXFy = 0;
  Fx = 0; Fy =0;
  for(int i = 0; i < NumPoints; i++){
    dF_drag(Points[2*i], Points[2*i + 1], VecArea[2*i], VecArea[2*i + 1],tau, AUXFx, AUXFy);
    Fx += AUXFx;
    Fy += AUXFy;
  }
}

void PartiPresion(int NumPoints, double* Points,  double* VecArea, double R, int ixc, int iyc){
  double VecNorma = 2*R*std::tan(M_PI/NumPoints);
  double x = 0; double y = 0;
  
  for(int i = 0; i < NumPoints; i++){  
    x = std::cos(M_PI*2/NumPoints*i); //componentes del vector unitario desde el
    y = std::sin(M_PI*2/NumPoints*i); //centro del cilindro hacia el punto
    Points[i*2] = ixc + R*x;
    Points[i*2 + 1] = iyc + R*y;
    VecArea[i*2] = VecNorma*x;
    VecArea[i*2 + 1] = VecNorma*y;
  }  
}

//------------------- Funciones Globales ------------

int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=3000;//=10000;
  double rho0=1.0,Ufan0=0.1;

  double tau = 1.5;
  int ixc=128, iyc=32;
  double R = 8.0;
  double omega = 2*M_PI/1000;
  int NumPoints = 24;
  double Points [NumPoints*2] {0};
  double VecArea [NumPoints*2] {0};
  PartiPresion(NumPoints, Points, VecArea,R,ixc,iyc);

  double FXdrag = 0;
  double FYdrag = 0;

  ofstream Fuerza("FuerzasDragTiempo-2.txt");

  
  //Start
  Aire.Start(rho0,Ufan0,0);
  //Run
  for(t=0;t<tmax;t++){
    Aire.Collision(tau);
    Aire.ImposeFields(Ufan0,omega,R, ixc,iyc);
    Aire.FuerzaTotalDrag(NumPoints, Points, VecArea, tau,FXdrag, FYdrag);
    Fuerza <<t << "\t" << FXdrag<<"\t"<< FYdrag<<std::endl;
    Aire.Advection();
  }

  return 0;
}  
   //  std::cout<<"La Fuerza de arrastre es F = ("<<FXdrag<<" , "<<FYdrag<<")"<<std::endl;
