#include <iostream>
#include "RK4CoupledStep.hpp"

//declaration of the functions defining the system
double f1(long double x, long double x1, long double x2);
double f2(long double x, long double x1, long double x2);

int main(void){
    //print with scientific long double precision
    std::cout.precision(14); std::cout.setf(std::ios::scientific);
    long double x, x1, x2, x0, x10, x20;
    //Initial conditions
    x0 = 1e-15;
    x10 = 1; //function
    x20 = 0; //derivative
    //Step and final x
    long double dx = 0.01;
    long double xf = 15;
    //Set initial conditions
    x=x0; x1=x10; x2=x20;
    //run Simulation
    while( x<xf+dx){
        RK4CoupledStep(f1, f2, x, x1, x2, dx);
        std::cout << x << "\t" << x1 << "\t" << x1 << std::endl;
    }
    return 0;
}


double f1(long double x, long double x1, long double x2){
    return x2/x;
}

double f2(long double x, long double x1, long double x2){
    return -x1*x;
}
