#include <iostream>
#include "Matrix_Multiplication.hpp"
#include <cmath>
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
)
{
  std::vector<decimal> P(sz);
  for (integer i = 0 ; i < sz ; i++)
  {
    decimal prod = 0.0;
    for (integer j = 0 ; j < sz ; j++)
    {
      prod += (A[i][j] * B[j]);
    }
    P[i] = prod;
  }

  return P;
}

template<class decimal, class integer>
std::vector<decimal> multiply_M_V
(
  const std::vector<std::vector<decimal > >& A,
  const std::vector<decimal >& B,
  integer sz1,
  integer sz2
)
{
  std::vector<decimal> P(sz1);
  for (integer i = 0 ; i < sz1 ; i++)
  {
    decimal prod = 0.0;
    for (integer j = 0 ; j < sz2 ; j++)
    {
      prod += (A[i][j] * B[j]);
    }
    P[i] = prod;
  }

  return P;
}

template<class decimal, class integer>
decimal multiply_V_V
(
  const std::vector<decimal>& A,
  const std::vector<decimal>& B,
  integer sz
)
{
  decimal value = 0.0;
  for (integer i = 0 ; i < sz ; i++)
  {
    value += A[i] * B[i];
  }
  return value;
}

template<class decimal, class integer>
std::vector<std::vector<decimal > > multiply_M_M
(
  const std::vector<std::vector<decimal > >& D,  // matrix D
  const std::vector<std::vector<decimal > >& X,  // vector X
  integer                               N,  // length of vector X
  integer                               M
)
{
  std::vector<std::vector<decimal > > prod(N, std::vector <decimal > (M));

  for ( integer i = 0; i < N; i++ )
  {
    for ( integer j = 0; j < M; j++ )
    {
      prod[i][j] = 0.0;
      for (integer k = 0 ; k < N ; k++)
      {
        prod[i][j] = prod[i][j] + (D[i][k] * X[k][j]);
      }
    }
    //D.clear();
    //X.clear();
  }
  return prod;
}

