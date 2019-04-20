#ifndef __Differentiation_Matrix__
#define __Differentitaion_Matrix__

//#include "IpTNLP.hpp"
#include <iostream>
#include <vector>
//using namespace Ipopt;
typedef double decimal;
typedef int integer;
template<class decimal, class integer>
std::vector<decimal > compute_c
(
    integer N
);

template<class decimal, class integer>
std::vector<decimal > define_time_stamps
(
    integer N
);

template<class decimal, class integer>
std::vector<decimal > multiply_D_X
(
    std::vector<std::vector<decimal > > 	D,	// matrix D
    std::vector<decimal > 								X,	// vector X
    integer 																N 	// length of vector X
);

template<class decimal, class integer>
std::vector<std::vector<decimal > > formulate_differentiation_matrix
(
    std::vector<decimal > c, //
    std::vector<decimal > t,
    integer 				 N //
);

#endif