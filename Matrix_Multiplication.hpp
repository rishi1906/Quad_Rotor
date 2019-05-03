#ifndef __Matrix_Multiplication__
#define __Matrix_Multiplication__

#include <iostream>
#include <vector>
// uncomment the typdefs when testing the integral weights separately
typedef double decimal;
typedef int integer;
//typedef Number decimal;
//typedef Index integer;
template<class decimal, class integer>
std::vector<decimal> multiply_M_V
(
    const std::vector<std::vector<decimal > >& A,
    const std::vector<decimal >& B,
    integer sz
);
template<class decimal, class integer>
std::vector<decimal> multiply_M_V
(
    const std::vector<std::vector<decimal > >& A,
    const std::vector<decimal >& B,
    integer sz1,
    integer sz2
);
template<class decimal, class integer>
decimal multiply_V_V
(
    const std::vector<decimal>& A,
    const std::vector<decimal>& B,
    integer sz
);

template<class decimal, class integer>
std::vector<std::vector<decimal > > multiply_D_X
(
    const std::vector<std::vector<decimal > >& D,  // matrix D
    const std::vector<decimal >& X,  // vector X
    integer                                 N,  // length of vector X
    integer                                 M
);


#endif


