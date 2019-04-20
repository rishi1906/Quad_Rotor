vector<Number> multiply_M_V
(
  vector<vector<Number > >Q,
  vector<Number > X,
  Index sz
)
{
  vector<Number> P;
  for (Index i = 0 ; i < sz ; i++)
  {
    Index prod = 0.0;
    for (Index j = 0 ; j < sz ; j++)
    {
      prod += (Q[i][j] * X[j]);
    }
    P[i] = prod;
  }

  return P;
}

Number multiply_V_M
(
  vector<Number> A,
  vector<Number> B,
  Index sz
)
{
  Number value = 0.0;
  for (Index i = 0 ; i < sz ; i++)
  {
    value += A[i] * B[i];
  }
  return value;
}