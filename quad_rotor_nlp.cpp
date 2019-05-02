// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "quad_rotor_nlp.hpp"
#include "Differentiation_Matrix.cpp"
#include "_Index_.cpp"
#include "Matrix_Multiplication.cpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include <map>
#include <string>
#include <set>
using namespace Ipopt;
#define PI acos(-1)
#define nxt_indx (((++idx) * N) + idx)
#define db_min std::numeric_limits<double>::min()
#define db_max std::numeric_limits<double>::max()
const Number t_0 = 0.00;
const Number t_f = 60.0;
const no_of_stt_var = 12 ; // define the length of X
const no_of_ctrl_var = 4 ; // define the length of U
// *  X = [p q r phi theta psi z Vz y Vy x Vx]'
// *  U = [netT Mx My Mz]'
// *  J = 1/2 integral(t_0,t_f,(X'.Q.X+U'.R.U))

// define no of constraints and size of decision vector

//#define N_ 10 // no of grid points

// define size of stepsize for finite difference scheme to find gradient
//const Number step_size = 1e-8;

#define PI acos(-1)
//Number RAD_to_DEG(Number r) { return r * 180.0 / PI; }

/*
// change me : define nxtutation of objective function to be used in finite differenc scheme
inline Number Obj_func(Number* X, Index n)
{

  Number value = X[n - 1];

  return value;
}*/

//Default Constructor
QUAD_ROTOR_NLP::QUAD_ROTOR_NLP
(
  Index N
)
{
  // define time stemps
  T = define_time_stamps<Number, Index>(N + 1);
  // test if time stamps are generated correctly
  /*
   for (Index i = 0; i <= N_; i++)
   {
     std::cout << "T[" << i << "]" << " : " << T[i] << "\n";
   }
  */

  // get initial indices of state variables
  INDX = get_indices(N);
  _Index_.insert(INDX.begin(), INDX.end());

  // test if _Index_ has the right vlaues
  /*
  for (auto const &pair : set)
  {
    std::cout << '{' << pair.first << "->" << pair.second << '}' << '\n';
  }
  */
}

// constructor
QUAD_ROTOR_NLP::QUAD_ROTOR_NLP() {}

// destructor
QUAD_ROTOR_NLP::~QUAD_ROTOR_NLP() {}


// returns the size of the problem
bool QUAD_ROTOR_NLP::get_nlp_info
(
  Index & n,         // size of problem
  Index & m,         // no of constraints
  Index & nnz_jac_g, // no of non zero elements in jacobain
  Index & nnz_h_lag, // no of non zero elements in hessian
  IndexStyleEnum & index_style
)
{
  // size of problem
  n = ((no_of_stt_var + no_of_ctrl_var) * (N + 1));

  // size of constraints
  m = ((N + 1) * no_of_stt_var);

  // size of jacobian matrix
  nnz_jac_g = m * n;

  // size of hessian of lagrangian
  // nnz_h_lag = 10;

  // style of indexing
  index_style = TNLP::C_STYLE; // use the C style indexing (0-based)

  return true;
}

// returns the variable bounds
bool QUAD_ROTOR_NLP::get_bounds_info
(
  Index n,      // size of problem
  Number * x_l, // lower limits for decision variables
  Number * x_u, // uppe limits for decision variables
  Index m,      // no of constraints
  Number * g_l, // lower limits for constraints
  Number * g_u  // upper limits for constraints
) // upper limits for constraints
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.

  //assertain the values of m and n
  assert(n == ((no_of_stt_var + no_of_ctrl_var) * (N + 1)));
  assert(m == ((N + 1) * no_of_stt_var));

  // Lower bounds
  for (int i = 0 ; i <= N ; i++)
  {
    x_l[_Index_["p"     ]  + i] = db_min;
    x_l[_Index_["q"     ]  + i] = db_min;
    x_l[_Index_["r"     ]  + i] = db_min;
    x_l[_Index_["phi"   ]  + i] = -2.0 * PI;
    x_l[_Index_["theta" ]  + i] = -2.0 * PI;
    x_l[_Index_["psi"   ]  + i] = -2.0 * PI;
    x_l[_Index_["z"     ]  + i] = db_min;
    x_l[_Index_["Vz"    ]  + i] = db_min;
    x_l[_Index_["y"     ]  + i] = db_min;
    x_l[_Index_["Vy"    ]  + i] = db_min;
    x_l[_Index_["x"     ]  + i] = db_min;
    x_l[_Index_["Vx"    ]  + i] = db_min;
    x_l[_Index_["netT"  ]  + i] = db_min;
    x_l[_Index_["Mx"    ]  + i] = db_min;
    x_l[_Index_["My"    ]  + i] = db_min;
    x_l[_Index_["Mz"    ]  + i] = db_min;
  }

  // Upper Bounds
  for (int i = 0 ; i <= N ; i++)
  {
    x_u[_Index_["p"     ]  + i] = db_max;
    x_u[_Index_["q"     ]  + i] = db_max;
    x_u[_Index_["r"     ]  + i] = db_max;
    x_u[_Index_["phi"   ]  + i] = +2.0 * PI;
    x_u[_Index_["theta" ]  + i] = +2.0 * PI;
    x_u[_Index_["psi"   ]  + i] = +2.0 * PI;
    x_u[_Index_["z"     ]  + i] = db_max;
    x_u[_Index_["Vz"    ]  + i] = db_max;
    x_u[_Index_["y"     ]  + i] = db_max;
    x_u[_Index_["Vy"    ]  + i] = db_max;
    x_u[_Index_["x"     ]  + i] = db_max;
    x_u[_Index_["Vx"    ]  + i] = db_max;
    x_u[_Index_["netT"  ]  + i] = db_max;
    x_u[_Index_["Mx"    ]  + i] = db_max;
    x_u[_Index_["My"    ]  + i] = db_max;
    x_u[_Index_["Mz"    ]  + i] = db_max;
  }

  // set bounds on constraints for ineuality constraints
  return true;
}

// returns the initial point for the problem
bool QUAD_ROTOR_NLP::get_starting_point
(
  Index n,           //
  bool init_x,       //
  Number * x,        //
  bool init_z,       //
  Number * z_L,      //
  Number * z_U,      //
  Index m,           //
  bool init_lambda,  //
  Number * lambda
)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  std::vector<Number > time(N + 1);

  // set the physical time bounds
  for (Index i = 0; i <= N; i++)
  {
    time[i] = ((tf - t0) * (T[i]) + (tf + t0)) / 2.0;
    //std::cout << "time[" << i << "]" << " : " << time[i] << std::endl;
  }
  //myfile << "Initialization\n";
  for (int i = 0 ; i <= N ; i++)
  {
    x[_Index_["p"     ]  + i] = 0.0;
    x[_Index_["q"     ]  + i] = 0.0;
    x[_Index_["r"     ]  + i] = 0.0;
    x[_Index_["phi"   ]  + i] = 0.0;
    x[_Index_["theta" ]  + i] = 0.0;
    x[_Index_["psi"   ]  + i] = 0.0;
    x[_Index_["z"     ]  + i] = 0.0;
    x[_Index_["Vz"    ]  + i] = 0.0;
    x[_Index_["y"     ]  + i] = 0.0;
    x[_Index_["Vy"    ]  + i] = 0.0;
    x[_Index_["x"     ]  + i] = 0.0;
    x[_Index_["Vx"    ]  + i] = 0.0;
    x[_Index_["netT"  ]  + i] = 0.0;
    x[_Index_["Mx"    ]  + i] = 0.0;
    x[_Index_["My"    ]  + i] = 0.0;
    x[_Index_["Mz"    ]  + i] = 0.0;
  }
  x[n - 1] = t_f;  // initial guess for finial time (take some high values)
  return true;
}

// returns the value of the objective function
//J = (1 / 2) * integral(t_0, t_f, (X'.Q.X+U'.R.U))
bool QUAD_ROTOR_NLP::eval_f
(
  Index n,           //
  const Number * x,  //
  bool new_x,        //
  Number & obj_value //
)
{
  assert(n == ((no_of_stt_var + no_of_ctrl_var) * (N + 1)));

  Index sz_X = (no_of_stt_var * N) + no_of_stt_var;
  Index sz_U = (no_of_ctrl_var * N) + no_of_ctrl_var;

  //define weights
  std::vector<Number> w; // weights
  w = compute_integral_weights(N);

  obj_value = 0.0;
  for (Index k = 0 ; k <= N ; k++)
  {
    std::vector<Number > X(sz_X ), U(sz_U), QX(sz_X), RU(sz_U);
    std::vector<std::vector<Number > > Q(sz_X, std::vector<decimal > (sz_X)), R(sz_U, std::vector<decimal > (sz_U));

    Number XTQX = 0.0, UTRU = 0.0;

    // define X
    Index i = 0;
    X[  i] = x[_Index_["p"    ]  + k];
    X[++i] = x[_Index_["q"    ]  + k];
    X[++i] = x[_Index_["r"    ]  + k];
    X[++i] = x[_Index_["phi"  ]  + k];
    X[++i] = x[_Index_["theta"]  + k];
    X[++i] = x[_Index_["psi"  ]  + k];
    X[++i] = x[_Index_["z"    ]  + k];
    X[++i] = x[_Index_["Vz"   ]  + k];
    X[++i] = x[_Index_["y"    ]  + k];
    X[++i] = x[_Index_["Vy"   ]  + k];
    X[++i] = x[_Index_["x"    ]  + k];
    X[++i] = x[_Index_["Vx"   ]  + k];
    X[++i] = x[_Index_["netT" ]  + k];
    X[++i] = x[_Index_["Mx"   ]  + k];
    X[++i] = x[_Index_["My"   ]  + k];
    X[++i] = x[_Index_["Mz"   ]  + k];

    //define U
    for (Index i = 0 ; i < sz_U ; i++)
    {
      U[i] = x[i + sz_X];
    }
    //define Q
    for (Index i = 0 ; i < sz_X ; i++)
    {
      for (Index j = 0 ; j < sz_X ; j++)
      {
        if (i == j)
        {
          Q[i][j] = 1.0 / (x[i] * x[i]);
        } else
        {
          Q[i][j] = 0.0;
        }
      }
    }
    //define R
    for (Index i = 0 ; i < sz_U ; i++)
    {
      for (Index j = 0 ; j < sz_U ; j++)
      {
        if (i == j)
        {
          R[i][j] = 1.0 / (u[i] * u[i]);
        } else
        {
          R[i][j] = 0.0;
        }
      }
    }
    QX = multiply_M_V(Q, X, sz_X);
    RU = multiply_M_V(R, U, sz_U);
    XTQX = multiply_V_V(X, QX, sz_X);
    UTRU = multiply_V_V(U, RU, sz_U);


    //obj_value += (XTQX + UTRU) / 2.0;
    obj_value += ((XTQX + UTRU) * w(k));
    X.clear();
    U.clear();

    QX.clear();
    RU.clear();

    Q.clear();
    R.clear();

  }
  obj_value /= ((t_f - t_0) / 2.0);

  w.clear();
  return true;
}

/*
Number grad_at_x
(
  Number Obj_func(Number * X, Index n), //
  Number * X, //
  Index pos,//
  Index n, //
  Number h
)
{
  X[pos] = X[pos] + h;
  Number f_x_p_h = Obj_func(X, n);
  X[pos] = X[pos] - h;

  X[pos] = X[pos] - h;
  Number f_x_m_h = Obj_func(X, n);
  X[pos] = X[pos] + h;

  Number grad = (f_x_p_h - f_x_m_h) / (2.0 * h);
  return grad;
}
*/


// return the gradient of the objective function grad_{x} f(x)
bool QUAD_ROTOR_NLP::eval_grad_f
(
  Index n, //
  const Number * x, //
  bool new_x,       //
  Number * grad_f
)
{
  assert(n == ((no_of_stt_var + no_of_ctrl_var) * (N + 1)));

  for (Index i = 0; i <= n - 2; i++)
  {
    grad_f[i] = 0;
  }

  grad_f[n - 1] = 1;
  return true;
}

// return the value of the constraints: g(x)
bool QUAD_ROTOR_NLP::eval_g
(
  Index n,          //
  const Number * x, //
  bool new_x,       //
  Index m,          //
  Number * g
)
{
  assert(n == ((no_of_stt_var + no_of_ctrl_var) * (N + 1)));
  assert(m == ((N + 1) * no_of_stt_var));

  Index sz_X = (no_of_stt_var * N) + no_of_stt_var;
  Index sz_U = (no_of_ctrl_var * N) + no_of_ctrl_var;
  Number A[no_of_stt_var][no_of_stt_var], B[no_of_stt_var][no_of_ctrl_var];

  //Number tf = 1.25;

  std::vector<Number > C(N + 1);

  std::vector<std::vector<Number > > D(N + 1 , std::vector <Number > (N + 1)), DX(N + 1 , std::vector <Number > (no_of_stt_var));
  std::vector<std::vector<Number > > X, AX, BU;

  C = compute_c<Number, Index>(N + 1);
  D = formulate_differentiation_matrix<Number, Index>(C, T, N + 1);

  // form X
  for (Index k = 0 ; k <= N ; k++)
  {
    std::vector<Number > k_X(no_of_stt_var);
    Index i = 0;
    k_X[  i] = x[_Index_["p"    ]  + k];
    k_X[++i] = x[_Index_["q"    ]  + k];
    k_X[++i] = x[_Index_["r"    ]  + k];
    k_X[++i] = x[_Index_["phi"  ]  + k];
    k_X[++i] = x[_Index_["theta"]  + k];
    k_X[++i] = x[_Index_["psi"  ]  + k];
    k_X[++i] = x[_Index_["z"    ]  + k];
    k_X[++i] = x[_Index_["Vz"   ]  + k];
    k_X[++i] = x[_Index_["y"    ]  + k];
    k_X[++i] = x[_Index_["Vy"   ]  + k];
    k_X[++i] = x[_Index_["x"    ]  + k];
    k_X[++i] = x[_Index_["Vx"   ]  + k];
    // push kthX to X
    X.push_back(k_X);
    k_X.clear();
  }

  // matrix multiply D and X
  DX = multiply_M_M(D, X, (N + 1), no_of_stt_var);

  // form AX
  // read A
  std::ofstream myfile;
  myfile.open("./Inputs/A.txt");
  for (integer i = 0 ; i < no_of_stt_var ; i++)
  {
    for (integer j = 0 ; j < no_of_stt_var ; j++)
    {
      myfile >> A[i][j];
    }
  }
  myfile.close();

  for (Index k = 0 ; k <= N ; k++)
  {
    std::vector<Number > k_X(no_of_stt_var);
    Index i = 0;
    k_X[  i] = x[_Index_["p"    ]  + k];
    k_X[++i] = x[_Index_["q"    ]  + k];
    k_X[++i] = x[_Index_["r"    ]  + k];
    k_X[++i] = x[_Index_["phi"  ]  + k];
    k_X[++i] = x[_Index_["theta"]  + k];
    k_X[++i] = x[_Index_["psi"  ]  + k];
    k_X[++i] = x[_Index_["z"    ]  + k];
    k_X[++i] = x[_Index_["Vz"   ]  + k];
    k_X[++i] = x[_Index_["y"    ]  + k];
    k_X[++i] = x[_Index_["Vy"   ]  + k];
    k_X[++i] = x[_Index_["x"    ]  + k];
    k_X[++i] = x[_Index_["Vx"   ]  + k];

    AX.push_back(multiply_M_V(A, k_X, no_of_stt_var));
    k_X.clear();
  }

  // form BX
  // read B
  std::ofstream myfile;
  myfile.open("./Inputs/B.txt");
  for (integer i = 0 ; i < no_of_stt_var ; i++)
  {
    for (integer j = 0 ; j < no_of_ctrl_var ; j++)
    {
      myfile >> B[i][j];
    }
  }
  myfile.close();
  for (Index k = 0 ; k <= N ; k++)
  {
    std::vector<Number > k_U(no_of_ctrl_var);
    Index i = 0;
    k_U[++i] = x[_Index_["netT" ]  + k];
    k_U[++i] = x[_Index_["Mx"   ]  + k];
    k_U[++i] = x[_Index_["My"   ]  + k];
    k_U[++i] = x[_Index_["Mz"   ]  + k];

    BU.push_back(multiply_M_V(B, k_U, N + 1, no_of_ctrl_var));
    k_U.clear();
  }

  // define constraints
  Index nth = 0 ;
  Number c1 = 2.0 / (x[n - 1] - t_0);
  for (Index i = 0; i <= N ; i++)
  {
    for (Index j = 0; j < no_of_stt_var ; j++)
    {
      g[nth] = (c1 * DX[i][j]) - AX[i][j] - BU[i][j];
      nth += 1;
    }
  }
  //additional constraints
  //k = N;
  // g[nth++] = x[_Index_["p"    ]  + k] - /*final value of "p"    */;
  // g[nth++] = x[_Index_["q"    ]  + k] - /*final value of "q"    */;
  // g[nth++] = x[_Index_["r"    ]  + k] - /*final value of "r"    */;
  // g[nth++] = x[_Index_["phi"  ]  + k] - /*final value of "phi"  */;
  // g[nth++] = x[_Index_["theta"]  + k] - /*final value of "theta"*/;
  // g[nth++] = x[_Index_["psi"  ]  + k] - /*final value of "psi"  */;
  // g[nth++] = x[_Index_["z"    ]  + k] - /*final value of "z"    */;
  // g[nth++] = x[_Index_["Vz"   ]  + k] - /*final value of "Vz"   */;
  // g[nth++] = x[_Index_["y"    ]  + k] - /*final value of "y"    */;
  // g[nth++] = x[_Index_["Vy"   ]  + k] - /*final value of "Vy"   */;
  // g[nth++] = x[_Index_["x"    ]  + k] - /*final value of "x"    */;
  // g[nth++] = x[_Index_["Vx"   ]  + k] - /*final value of "Vx"   */;
  // g[nth++] = x[_Index_["netT" ]  + k] - /*final value of "netT" */;
  // g[nth++] = x[_Index_["Mx"   ]  + k] - /*final value of "Mx"   */;
  // g[nth++] = x[_Index_["My"   ]  + k] - /*final value of "My"   */;
  // g[nth++] = x[_Index_["Mz"   ]  + k] - /*final value of "Mz"   */;

  return true;
}


// return the structure or values of the Jacobian
bool QUAD_ROTOR_NLP::eval_jac_g
(
  Index n,          //
  const Number * x, //
  bool new_x,       //
  Index m,          //
  Index nele_jac,   //
  Index * iRow,     //
  Index * jCol,     //
  Number * values   //
)
{
  assert(n == ((no_of_stt_var + no_of_ctrl_var) * (N + 1)));
  assert(m == ((N + 1) * no_of_stt_var));

  if (values == NULL) {
    // return the structure of the Jacobian

    // this particular Jacobian is dense
    Index nnz = 0;
    for (Index i = 0; i < m; i++) {
      for (Index j = 0; j < n; j++) {
        iRow[nnz] = i;
        jCol[nnz] = j;
        nnz += 1;
      }
    }
  }
  // else
  // {
  //    // return the values of the Jacobian of the constraintss
  // }

  return true;
}

/*
//return the structure or values of the Hessian
bool QUAD_ROTOR_NLP::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
   )
{
   if( values == NULL )
   {
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.

      // the hessian for this problem is actually dense
      Index idx = 0;
      for( Index row = 0; row < 4; row++ )
      {
         for( Index col = 0; col <= row; col++ )
         {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      assert(idx == nele_hess);
   }
   else
   {
      // return the values. This is a symmetric matrix, fill the lower left
      // triangle only

      // fill the objective portion

   }

   return true;
}
*/

void QUAD_ROTOR_NLP::finalize_solution
(
  SolverReturn status,       //
  Index n,                   //
  const Number * x,          //
  const Number * z_L,        //
  const Number * z_U,        //
  Index m,                   //
  const Number * g,          //
  const Number * lambda,     //
  Number obj_value,          //
  const IpoptData * ip_data, //
  IpoptCalculatedQuantities * ip_cq
)
{
  // here is where we would store the solution to variables, or write to a
  // file,
  // etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;

  // std::cout << std::endl << "Time t " << std::endl;
  // for (Index i = 0; i < n / 3; i++) {
  //   std::cout << "t[" << i << "] = " << i << std::endl;
  // }
  // myfile << "Writing this to a file.\n";
  string state_var[no_of_stt_var + no_of_ctrl_var] = {"p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx", "netT", "Mx", "My", "Mz"};
  std::vector<Number > time(N + 1);
  //Number tf = 1.25;
  for (Index i = 0; i <= N; i++)
  {
    time[i] = (x[n - 1] / 2.0) * (T[i] + 1.0);
    std::cout << "t[" << i << "]" << " : " << time[i] << std::endl;
  }

  for (Index i = 0 ; i < no_of_stt_var ; i++)
  {
    std::ofstream myfile;

    myfile.open(state_var[i] + ".txt");
    //std::cout << state_var[i] << "\n";
    for (Index k = 0; k <= N_; k++) {
      //std::cout << x[i] << std::endl;
      myfile << time[k] << "," << x[_Index_[state_var[i]]  + k] << "\n";
    }
    myfile.close();
  }

  std::cout << x[n - 1] << std::endl;

  // std::cout << std::endl << std::endl << "Solution of the bound
  // multipliers,
  // z_L and z_U" << std::endl;
  // for ( Index i = 0; i < n; i++ )
  // {
  //    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  // }
  // for ( Index i = 0; i < n; i++ )
  // {
  //    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  // }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::cout << std::endl
            << "Final value of the constraints:" << std::endl;
  for (Index i = 0; i < m; i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }

  // Analytical Solution
  /*std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;
  */
}
