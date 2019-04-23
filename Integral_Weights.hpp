#ifndef __INTEGRAL_WEIGHTS__
#define __INTEGRAL_WEIGHTS__

#include <iostream>
#include <vector>

typedef double decimal;
typedef int integer;

template<class decimal, class integer>
std::vector<decimal > compute_integral_weights
(
    integer N_
);

#endif