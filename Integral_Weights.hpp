#ifndef __INTEGRAL_WEIGHTS__
#define __INTEGRAL_WEIGHTS__

#include <iostream>
#include <vector>
// uncomment the typdefs when testing the integral weights separately
typedef double decimal;
typedef int integer;
//typedef Number decimal;
//typedef Index integer;
template<class decimal, class integer>
std::vector<decimal > compute_integral_weights
(
    integer N_
);

#endif