vector<Number> multiply_M_V
(
  vector<vector<Number > >A,
  vector<Number > B,
  Index sz
)
{
  vector<Number> P;
  for (Index i = 0 ; i < sz ; i++)
  {
    Index prod = 0.0;
    for (Index j = 0 ; j < sz ; j++)
    {
      prod += (A[i][j] * B[j]);
    }
    P[i] = prod;
  }

  return P;
}

Number multiply_V_V
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