#pragma once
#include <functional>
#include <cmath>
void RK4CoupledStep( std::function <long double (long double, long double, long double)> f1,
                     std::function <long double (long double, long double, long double)> f2 , long double &t, long double &x,
                     long double &y, long double dt);
