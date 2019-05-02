#include <iostream>
#include <cmath>
#include <vector>
#include "Integral_Weights.hpp"
#define PI acos(-1)
using namespace std;
// uncomment the typdefs when testing the integral weights separately
//typedef double decimal;
//typedef int integer;


template<class decimal, class integer>
std::vector<decimal > compute_integral_weights(integer N_)
{
	//cout << "Entered\n";
	integer N = N_ - 1;
	std::vector<decimal > w(N);

	if (N % 2 == 0)
	{
		w[0] = 1.0 / ((N * N) - 1.0);

		w[N] = 1.0 / ((N * N) - 1.0);
		for (integer s = 1 ; s <= N / 2 ; s++)
		{
			decimal sum_trm = 0.0;
			for (integer j = 0 ; j <= N / 2 ; j++)
			{
				decimal val = ( 1.0 / (1.0 - (4.0 * j * j)) ) * cos((2.0 * PI * j * s) / N );
				if (j == 0 || j == N / 2)
				{
					sum_trm += val / 2.0;
				}
				else
				{
					sum_trm += val;
				}
			}
			w[s] = w[N - s] = (4.0 / N) * sum_trm;
		}
	}

	else if (N  % 2 != 0)
	{
		w[0] = w[N] = 1.0 / (N * N);
		for (integer s = 1 ; s <= (N - 1) / 2 ; s++)
		{
			decimal sum_trm = 0.0;
			for (integer j = 0 ; j <= (N - 1) / 2 ; j++)
			{
				decimal val = ( 1.0 / (1.0 - (4.0 * j * j)) ) * cos((2.0 * PI * j * s) / N );
				if (j == 0 || j == ((N - 1) / 2) )
				{
					sum_trm += val / 2.0;
				}
				else
				{
					sum_trm += val;
				}
			}
			w[s] = w[N - s] = (4.0 / N) * sum_trm;
		}
	}
	return w;
}
