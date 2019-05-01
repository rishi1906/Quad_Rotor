#ifndef __Matrix_Multiplication__
#define __Matrix_Multiplication__

#include <iostream>
#include <vector>

typedef double decimal;
typedef int integer;
template<class decimal, class integer>
vector<decimal> multiply_M_V
(
    vector<vector<decimal > >A,
    vector<decimal > B,
    integer sz
);
template<class decimal, class integer>
vector<decimal> multiply_M_V
(
    vector<vector<decimal > > A,
    vector<decimal >          B,
    integer sz1,
    integer sz2
);
template<class decimal, class integer>
decimal multiply_V_V
(
    vector<decimal> A,
    vector<decimal> B,
    integer sz
);

template<class decimal, class integer>
std::vector<std::vector<decimal > > multiply_D_X
(
    std::vector<std::vector<decimal > > 	D,	// matrix D
    std::vector<decimal > 					X,	// vector X
    integer 								N, 	// length of vector X
    integer 								M
);


#endif


