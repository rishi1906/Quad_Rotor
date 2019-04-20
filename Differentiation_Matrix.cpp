#include <iostream>
#include "Differentiation_Matrix.hpp"
#include <cmath>
#include <vector>
#define PI acos(-1)
typedef double decimal;
typedef int integer;
//decimal DEG_to_RAD(decimal d) { return d * PI / 180.0; }
//decimal RAD_to_DEG(decimal r) { return r * 180.0 / PI; }

template<class decimal, class integer>
std::vector<decimal > compute_c
(
    integer N
)
{
	std::vector<decimal > c(N);
	c[0] = c[N - 1] = 2;
	for (integer i = 1 ; i < N - 1 ; i++) {
		c[i] = 1;
	}
	return c;
}

template<class decimal, class integer>
std::vector<decimal > define_time_stamps
(
    integer N
)
{
	std::vector<decimal > T(N);
	for (integer k = 0 ; k < N ; k++) {
		T[k] = cos(((PI * k ) / (N - 1) ));
		//T[k] = -(1.0) *T[k];
	}
	/*
	decimal tf = 2.00;
	for (integer i = 0; i < N; i++)
	{
		t[i] = (tf / 2.0) * (t[i] + 1);
		//std::cout << "t[" << i << "]" << " : " << t[i] << std::endl;
	}*/
	return T;
}

template<class decimal, class integer>
std::vector<decimal > multiply_D_X
(
    std::vector<std::vector<decimal > > 	D,	// matrix D
    std::vector<decimal > 					X,	// vector X
    integer 								N 	// length of vector X
)
{
	std::vector<decimal > prod(N);

	for ( integer i = 0; i < N; i++ )
	{
		prod[i] = 0.0;                        // <===== Needs initialising
		for ( integer j = 0; j < N; j++ )
		{
			prod[i] = prod[i] + (D[i][j] * X[j]);     // <===== Add terms to sum for ith element
		}
	}
	//D.clear();
	//X.clear();
	return prod;
}

template<class decimal, class integer>
std::vector<std::vector<decimal > > formulate_differentiation_matrix
(
    std::vector<decimal > c, //
    std::vector<decimal > t, //
    integer 			  N  // N = N_ + 1
)
{
	std::vector<std::vector<decimal > > D (N , std::vector <decimal > (N));
	for (integer k = 0 ; k < N; k++) {
		for (integer j = 0 ; j < N ; j++)
		{
			if (j == k)
			{
				if (j == 0)
				{
					D[k][j] = ((2.0 * (N - 1) * (N - 1)) + 1) / 6.0;
				} else if (j == (N - 1) )
				{
					D[k][j] = -1.0 * (((2.0 * (N - 1) * (N - 1)) + 1) / 6.0);
				} else
				{
					D[k][j] = (-1.0 * t[k]) / (2.0 * (1 - (t[k] * t[k])));
				}
			} else if (j != k)
			{
				D[k][j] = (c[k] / c[j]) * ( (pow(-1, (j + k))) / (t[k] - t[j]) );
			}
			//D[j][k] = -D[j][k];
		}
	}
	/*
	for (integer k = 0 ; k < N; k++) {
		for (integer j = 0 ; j < N ; j++)
		{
			D[k][j] = (-1.0)*(D[k][j]);
		}
	}*/
	return D;
}
