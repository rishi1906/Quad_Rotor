#include <iostream>
#include "Matrix_Multiplication.hpp"
#include <cmath>
#include <vector>

typedef double decimal;
typedef int integer;

template<class decimal, class integer>
vector<decimal> multiply_M_V
(
  vector<vector<decimal > > A,
  vector<decimal >          B,
  integer sz
)
{
  vector<decimal> P;
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
vector<decimal> multiply_M_V
(
  vector<vector<decimal > > A,
  vector<decimal >          B,
  integer sz1,
  integer sz2
)
{
  vector<decimal> P;
  for (integer i = 0 ; i < sz1 ; i++)
  {
    integer prod = 0.0;
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
  vector<decimal> A,
  vector<decimal> B,
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
  std::vector<std::vector<decimal > >   D,  // matrix D
  std::vector<std::vector<decimal > >   X,  // vector X
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
        prod[i][j] = prod[i][j] + (D[i][k] * X[j][k]);
      }
    }
    //D.clear();
    //X.clear();
    return prod;
  }
}

