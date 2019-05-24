// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16
// Authors: Rishabh Vashistha rishabhsharma1906@gmail.com
//
/**
 *
 *  Direct Trajectory Optimization by a Chebyshev
 *             Pseudospectral Method
 *            Quadrotor Modelling and Dynamics
 *
 *  X = [p q r phi theta psi z  Vz  y  Vy   x  Vx]'
 *       1 2 3 4   5     6   7  8   9  10  11  12
 *  U = [netT Mx My Mz]'
 *  J = integral(t_0,t_f,(X'.Q.X+U'.R.U))
 *  Q and R are taken as identity matrices
 */
#include "quad_rotor_nlp.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <map>
#include <utility>
#include <string>
#include <cstring>
#include <set>
#include "Differentiation_Matrix.cpp"
#include "Matrix_Multiplication.cpp"
#include "Integral_Weights.cpp"
using namespace Ipopt;
#define PI acos(-1)
#define nxt_indx (((idx+1) * N) + (idx+1))
#define rndm ((mx - mn) * ( (double)rand() / (double)RAND_MAX ) + mn)
#define sz_n ((n_s + n_c) * (N + 1))
#define sz_m (((N + 1) * n_s) + n_s)
// define size of stepsize for finite difference scheme to find gradient
const Number step_size = 1e-8;
const Number t_0 = 0.00;
const Number t_f = 12.00;
const Number mx = 999.999;
const Number mn = -999.999;
const Index n_s = 12 ;  // define the length of X
const Index n_c = 4 ;   // define the length of U
const Index add_frc_const = 4;

// constants of the dynamic system
//****************************************************************
const Number mass = 1.104; // mass
const Number grav = 9.81;  // gravity
// const Number density = 1.225;
const Number Ix  = 0.008562874765838073;
const Number Iy  = 0.008788914621963906;
const Number Iz  = 0.015570395039175332;
// const Number momentArm   = 0.225; %half of quadcopter diagonal
const Number tow_x = 0.0;
const Number tow_y = 0.0;
const Number tow_z = 0.0;
const Number tow_wx = 0.0;
const Number tow_wy = 0.0;
const Number tow_wz = 0.0;
const Number f_wx = 0.0;
const Number f_wy = 0.0;
const Number f_wz = 0.0;
const Number f_t = mass * grav;
//*****************************************************************

// const Number db_min = std::numeric_limits<double>::min();
// const Number db_max = std::numeric_limits<double>::max();

const Number db_min = -99999.00;
const Number db_max = 99999.00;

// *  X = [p q r phi theta psi z Vz y Vy x Vx]' // these are the pertubrations
// *  U = [netT Mx My Mz]'  // these are the pertubrations
// *  J = integral(t_0,t_f,(X'.Q.X+U'.R.U))

std::map <string, int> get_indices(int N)
{
  int idx = 0;
  std::map <string, int> mp;
  mp.insert(make_pair("p"    , idx));
  mp.insert(make_pair("q"    , nxt_indx));
  idx += 1;
  mp.insert(make_pair("r"    , nxt_indx));
  idx += 1;
  mp.insert(make_pair("phi"  , nxt_indx));
  idx += 1;
  mp.insert(make_pair("theta", nxt_indx));
  idx += 1;
  mp.insert(make_pair("psi"  , nxt_indx));
  idx += 1;
  mp.insert(make_pair("z"    , nxt_indx));
  idx += 1;
  mp.insert(make_pair("Vz"   , nxt_indx));
  idx += 1;
  mp.insert(make_pair("y"    , nxt_indx));
  idx += 1;
  mp.insert(make_pair("Vy"   , nxt_indx));
  idx += 1;
  mp.insert(make_pair("x"    , nxt_indx));
  idx += 1;
  mp.insert(make_pair("Vx"   , nxt_indx));
  idx += 1;
  mp.insert(make_pair("netT" , nxt_indx));
  idx += 1;
  mp.insert(make_pair("Mx"   , nxt_indx));
  idx += 1;
  mp.insert(make_pair("My"   , nxt_indx));
  idx += 1;
  mp.insert(make_pair("Mz"   , nxt_indx));
  return mp;
};

//Default Constructor
QUAD_ROTOR_NLP::QUAD_ROTOR_NLP
(
  Index N_
)
{
  N = N_;
  // define time stemps
  std::ofstream file;
  file.open("output.txt");

  T = define_time_stamps<Number, Index>(N + 1);
  // test if time stamps are generated correctly

  for (Index i = 0; i <= N; i++)
  {
    //std::cout << "T[" << i << "]" << " : " << T[i] << "\n";
    file << T[i] << "\n";
  }
  file.close();

  // get initial indices of state variables
  INDX = get_indices(N);
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
  n = sz_n;

  // size of constraints
  m = sz_m;

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
  assert(n == sz_n);
  assert(m == sz_m);

  // Lower Bounds
  for (int i = 0 ; i <= N ; i++)
  {
    x_l[INDX["p"     ]  + i] = db_min;
    x_l[INDX["q"     ]  + i] = db_min;
    x_l[INDX["r"     ]  + i] = db_min;
    x_l[INDX["phi"   ]  + i] = db_min;
    x_l[INDX["theta" ]  + i] = db_min;
    x_l[INDX["psi"   ]  + i] = db_min;
    x_l[INDX["z"     ]  + i] = db_min;
    x_l[INDX["Vz"    ]  + i] = db_min;
    x_l[INDX["y"     ]  + i] = db_min;
    x_l[INDX["Vy"    ]  + i] = db_min;
    x_l[INDX["x"     ]  + i] = db_min;
    x_l[INDX["Vx"    ]  + i] = db_min;
    x_l[INDX["netT"  ]  + i] = db_min;
    x_l[INDX["Mx"    ]  + i] = db_min;
    x_l[INDX["My"    ]  + i] = db_min;
    x_l[INDX["Mz"    ]  + i] = db_min;
  }
  /*
  for (int i = 0 ; i < n ; i++)
  {
    x_l[i] = db_min;
  }
  */
  // Upper Bounds

  for (int i = 0 ; i <= N ; i++)
  {
    x_u[INDX["p"     ]  + i] = db_max;
    x_u[INDX["q"     ]  + i] = db_max;
    x_u[INDX["r"     ]  + i] = db_max;
    x_u[INDX["phi"   ]  + i] = db_max;
    x_u[INDX["theta" ]  + i] = db_max;
    x_u[INDX["psi"   ]  + i] = db_max;
    x_u[INDX["z"     ]  + i] = db_max;
    x_u[INDX["Vz"    ]  + i] = db_max;
    x_u[INDX["y"     ]  + i] = db_max;
    x_u[INDX["Vy"    ]  + i] = db_max;
    x_u[INDX["x"     ]  + i] = db_max;
    x_u[INDX["Vx"    ]  + i] = db_max;
    x_u[INDX["netT"  ]  + i] = db_max;
    x_u[INDX["Mx"    ]  + i] = db_max;
    x_u[INDX["My"    ]  + i] = db_max;
    x_u[INDX["Mz"    ]  + i] = db_max;
  }

  /*
  for (int i = 0 ; i < n ; i++)
  {
  x_u[i] = db_max;
  }
  */
  // set bounds on constraints for inequality constraints
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
    time[i] = ((t_f - t_0) * (T[i]) + (t_f + t_0)) / 2.0;
    //std::cout << "time[" << i << "]" << " : " << time[i] << std::endl;
  }

  std::ifstream f1;
  f1.open("./Inputs/x.txt");
  for (integer i = 0 ; i < n ; i++)
  {
    f1 >> x[i];
    //x[i] += 0.5;
    //myfile << x[i] << " ";
  }
  f1.close();
  /* take input of guess_stt

  string guess_stt = "z";

  std::ifstream rd_z;
  string file_name = "./Inputs/guess_" + guess_stt + ".txt";
  const char* fn = file_name.c_str();
  rd_z.open(fn);
  for (Index k = N; k >= 0; k--) {
    //std::cout << x[i] << std::endl;
    rd_z >> x[INDX["z"]  + k];
  }
  rd_z.close();
  */
  //Index cnt_dwn=1;
  string state_var[n_s + n_c] = {"p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx", "netT", "Mx", "My", "Mz"};
  for (Index i = 0 ; i < n_s + n_c ; i++)
  {
    string file_name = "./Outputs/" + state_var[i] + "_inp.txt";
    std::ofstream myfile;
    const char* fn = file_name.c_str();
    myfile.open(fn);
    //std::cout << state_var[i] << "\n";
    for (Index k = N; k >= 0; k--) {
      //std::cout << x[i] << std::endl;
      myfile << time[k] << "," << x[INDX[state_var[i]]  + k] << "\n";
    }
    myfile.close();
  }
  //myfile.close();
  /*
  Number prtb = 1.0;
  Number del_prtb = 0.1;
  for (int i = 0 ; i <= N ; i++)
  {
    x[INDX["p"     ]  + i] = prtb - i * del_prtb;
    myfile << "p[" << i << "] : " << x[INDX["p"     ]  + i];
    x[INDX["q"     ]  + i] = prtb - i * del_prtb;
    x[INDX["r"     ]  + i] = prtb - i * del_prtb;
    x[INDX["phi"   ]  + i] = prtb - i * del_prtb;
    x[INDX["theta" ]  + i] = prtb - i * del_prtb;
    x[INDX["psi"   ]  + i] = prtb - i * del_prtb;
    x[INDX["z"     ]  + i] = prtb - i * del_prtb;
    x[INDX["Vz"    ]  + i] = prtb - i * del_prtb;
    x[INDX["y"     ]  + i] = prtb - i * del_prtb;
    x[INDX["Vy"    ]  + i] = prtb - i * del_prtb;
    x[INDX["x"     ]  + i] = prtb - i * del_prtb;
    x[INDX["Vx"    ]  + i] = prtb - i * del_prtb;
    x[INDX["netT"  ]  + i] = prtb - i * del_prtb;
    x[INDX["Mx"    ]  + i] = prtb - i * del_prtb;
    x[INDX["My"    ]  + i] = prtb - i * del_prtb;
    x[INDX["Mz"    ]  + i] = prtb - i * del_prtb;
  }*/


  //x[n - 1] = t_f;  // initial guess for finial time (take some high values)
  return true;
}

// change me : define computation of objective function to be used in finite difference scheme
Number QUAD_ROTOR_NLP::Obj_func
(
  Number* x,
  Index N
)
{

  //define weights
  std::vector<Number> w(N + 1); // weights
  w = compute_integral_weights<Number, Index>(N + 1);

  Number obj_value = 0.0;
  for (Index k = N ; k >= 0 ; k--)
  {
    std::vector<Number > X(n_s ), U(n_c), QX(n_s), RU(n_c);
    std::vector<std::vector<Number > > Q(n_s, std::vector<Number > (n_s)), R(n_c, std::vector<Number > (n_c));

    Number XTQX = 0.0, UTRU = 0.0;

    // define X
    Index i = 0;
    X[  i] = x[INDX["p"    ]  + k];
    X[++i] = x[INDX["q"    ]  + k];
    X[++i] = x[INDX["r"    ]  + k];
    X[++i] = x[INDX["phi"  ]  + k];
    X[++i] = x[INDX["theta"]  + k];
    X[++i] = x[INDX["psi"  ]  + k];
    X[++i] = x[INDX["z"    ]  + k];
    X[++i] = x[INDX["Vz"   ]  + k];
    X[++i] = x[INDX["y"    ]  + k];
    X[++i] = x[INDX["Vy"   ]  + k];
    X[++i] = x[INDX["x"    ]  + k];
    X[++i] = x[INDX["Vx"   ]  + k];

    //define U
    i = 0;
    U[  i] = x[INDX["netT" ]  + k];
    U[++i] = x[INDX["Mx"   ]  + k];
    U[++i] = x[INDX["My"   ]  + k];
    U[++i] = x[INDX["Mz"   ]  + k];

    //define Q
    for (Index i = 0 ; i < n_s ; i++)
    {
      for (Index j = 0 ; j < n_s ; j++)
      {

        if (i == j)
        {
          Q[i][j] = 1.0;
          //Q[i][j] = 1.0 / (X[i] * X[i]);
        } else
        {
          Q[i][j] = 0.0;
        }
      }
    }
    //define R
    for (Index i = 0 ; i < n_c ; i++)
    {
      for (Index j = 0 ; j < n_c ; j++)
      {
        if (i == j)
        {
          R[i][j] = 1.0;
          //R[i][j] = 1.0 / (U[i] * U[i]);
        } else
        {
          R[i][j] = 0.0;
        }
      }
    }
    QX = multiply_M_V<Number, Index>(Q, X, n_s);
    RU = multiply_M_V<Number, Index>(R, U, n_c);
    XTQX = multiply_V_V<Number, Index>(X, QX, n_s);
    UTRU = multiply_V_V<Number, Index>(U, RU, n_c);

    obj_value += ((XTQX + UTRU) * w[k]);
    X.clear();
    U.clear();

    QX.clear();
    RU.clear();

    Q.clear();
    R.clear();

  }
  obj_value *= ((t_f - t_0) / 2.0);
  w.clear();
  return obj_value;
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
  assert(n == sz_n);

  // declare X a copy array of same size of x
  Number X[n];
  // make a copy of x in X
  for (Index k = 0; k < n; k++)
  {
    X[k] = x[k];
  }
  obj_value = Obj_func(X, N);
  //cout<<obj_value<<"\n";
  return true;
}


Number QUAD_ROTOR_NLP::grad_at_x
(
  Number * x, //
  Index pos,//
  Index N, //
  Number h
)
{
  string state_var[n_s + n_c] = {"p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx", "netT", "Mx", "My", "Mz"};
  Index len =  n_s + n_c;

  for (Index nth = 0 ; nth < len ; nth++)
  {
    x[INDX[state_var[nth]]  + pos] += h;
  }
  Number f_x_p_h = Obj_func(x, N);

  for (Index nth = 0 ; nth < len ; nth++)
  {
    x[INDX[state_var[nth]]  + pos] -= h;
  }

  for (Index nth = 0 ; nth < len ; nth++)
  {
    x[INDX[state_var[nth]]  + pos] -= h;
  }

  Number f_x_m_h = Obj_func(x, N);

  Number grad = (f_x_p_h - f_x_m_h) / (2.0 * h);
  for (Index nth = 0 ; nth < len ; nth++)
  {
    x[INDX[state_var[nth]]  + pos] += h;
  }
  return grad;
}

// return the gradient of the objective function grad_{x} f(x)
bool QUAD_ROTOR_NLP::eval_grad_f
(
  Index n, //
  const Number * x, //
  bool new_x,       //
  Number * grad_f
)
{
  assert(n == sz_n);

  // Approx Gradient using FDS
  // declare X a copy array of same size of x
  Number X[n];
  // make a copy of x in X
  for (Index k = 0; k < n; k++)
  {
    X[k] = x[k];
  }
  Number h = step_size;

  for (Index at = 0; at < n; at++)
  {
    Number val = grad_at_x(X, at, N, h);
    grad_f[at] = val;
  }

  return true;
}

// return the value of the constraints
bool QUAD_ROTOR_NLP::eval_g
(
  Index n,          //
  const Number * x, //
  bool new_x,       //
  Index m,          //
  Number * g
)
{
  assert(n == sz_n);
  assert(m == sz_m);

  std::vector<Number > C(N + 1), X(n + 1), DX(N + 1);

  std::vector<std::vector<Number > > D(N + 1 , std::vector <Number > (N + 1));

  C = compute_c<Number, Index>(N + 1);
  D = formulate_differentiation_matrix<Number, Index>(C, T, N + 1);

  // string state_var[n_s + n_c] = {"p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx", "netT", "Mx", "My", "Mz"};
  // Index len =  n_s + n_c;

  Index nth = 0;
  Number c1 = 2.0 / (t_f - t_0);

/*
  constraint for p
  p_dot = (((Iy-Iz)*r*q)/Ix) + ((tow_x + tow_wx)/Ix);
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["p"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (((Iy - Iz) * x[INDX["r"] + k] * x[INDX["q"] + k]) / Ix)
             - ((x[INDX["Mx"] + k]) / Ix);
    // - ((tow_x + tow_wx) / Ix);
    nth += 1;
  }

  /*
  constraint for q
  q_dot = (((Iz-Ix)*p*r)/Iy) + ((tow_y + tow_wy)/Iy);
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["q"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (((Iz - Ix) * x[INDX["p"] + k] * x[INDX["r"] + k]) / Iy)
             - ((x[INDX["My"] + k]) / Iy);
    // - ((tow_y + tow_wy) / Iy);
    nth += 1;
  }


  /*
  constraint for r
  r_dot = (((Ix-Iy)*p*q)/Iz) + ((tow_z + tow_wz)/Iz);
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["r"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (((Ix - Iy) * x[INDX["p"] + k] * x[INDX["q"] + k]) / Iz)
             - ((x[INDX["Mz"] + k]) / Iz);
    // - ((tow_z + tow_wz) / Iz);
    nth += 1;
  }

  /*
  constraint for phi
  phi_dot = p + (r*cos(phi)*tan(theta)) + q*(sin(phi))*(tan(theta));
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["phi"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - x[INDX["p"] + k]
             - (x[INDX["r"] + k] * cos(x[INDX["phi"] + k]) * tan(x[INDX["theta"] + k]))
             - (x[INDX["q"] + k] * (sin(x[INDX["phi"] + k])) * (tan(x[INDX["theta"] + k])));
    nth += 1;
  }

  /*
  constraint for theta
  theta_dot = (q*cos(phi)) - (r*sin(phi))
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["theta"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["q"] + k] * cos(x[INDX["phi"] + k]))
             + (x[INDX["r"] + k] * sin(x[INDX["phi"] + k]));
    nth += 1;
  }

  /*
  constraint for psi
  psi_dot = ((r*cos(phi))/cos(theta)) + ((q*sin(phi))/cos(theta))
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["psi"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - ((x[INDX["r"] + k] * cos(x[INDX["phi"] + k])) / cos(x[INDX["theta"] + k]))
             - ((x[INDX["q"] + k] * sin(x[INDX["phi"] + k])) / cos(x[INDX["theta"] + k]));
    nth += 1;
  }

  /*
  constraint for x
  x_dot = (w*((sin(phi)*sin(psi))+(cos(phi)*cos(psi)*sin(theta)))) - (v*((cos(phi)*sin(psi))-(cos(psi)*sin(phi)*sin(theta)))) + (u*cos(psi)*cos(theta))
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["x"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["Vz"] + k] * ((sin(x[INDX["phi"] + k]) * sin(x[INDX["psi"] + k])) + (cos(x[INDX["phi"] + k]) * cos(x[INDX["psi"] + k]) * sin(x[INDX["theta"] + k]))))
             + (x[INDX["Vy"] + k] * ((cos(x[INDX["psi"] + k]) * sin(x[INDX["psi"] + k])) - (cos(x[INDX["psi"] + k]) * sin(x[INDX["phi"] + k]) * sin(x[INDX["theta"] + k]))))
             - (x[INDX["Vx"] + k] * cos(x[INDX["psi"] + k]) * cos(x[INDX["theta"] + k]));
    nth += 1;
  }
  /*
  constraint for y
  y_dot = (v*((cos(phi)*cos(psi))+(sin(phi)*sin(psi)*sin(theta)))) - (w*((cos(psi)*sin(phi))-(cos(phi)*sin(psi)*sin(theta)))) + (u*cos(theta)*sin(psi))
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["y"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["Vy"] + k] * ((cos(x[INDX["phi"] + k]) * cos(x[INDX["psi"] + k])) + (sin(x[INDX["phi"] + k]) * sin(x[INDX["psi"] + k]) * sin(x[INDX["theta"] + k]))))
             + (x[INDX["Vz"] + k] * ((cos(x[INDX["psi"] + k]) * sin(x[INDX["phi"] + k])) - (cos(x[INDX["phi"] + k]) * sin(x[INDX["psi"] + k]) * sin(x[INDX["theta"] + k]))))
             - (x[INDX["Vx"] + k] * cos(x[INDX["theta"] + k]) * sin(x[INDX["psi"] + k]));
    nth += 1;
  }
  /*
  constraint for z
  z_dot = (w*cos(phi)*cos(theta)) - (u*sin(theta)) + (v*cos(theta)*sin(phi))
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["z"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["Vz"] + k] * cos(x[INDX["phi"] + k]) * cos(x[INDX["theta"] + k]))
             + (x[INDX["Vx"] + k] * sin(x[INDX["theta"] + k]))
             - (x[INDX["Vy"] + k] * cos(x[INDX["theta"] + k]) * sin(x[INDX["phi"] + k]));
    nth += 1;
  }

  /*
  constraint for Vx(u)
  Vx(u)_dot = (r*v) - (q*w) - (g*sin(theta)) + (f_wx/m)
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["Vx"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["r"] + k] * x[INDX["Vy"] + k])
             + (x[INDX["q"] + k] * x[INDX["Vz" ] + k])
             + (grav * sin(x[INDX["theta"] + k]))
             - (f_wx / mass);
    nth += 1;
  }

  /*
  constraint for Vy(v)
  Vy(v)_dot = (p*w) - (r*u) + (g*sin(phi)*cos(theta)) + (f_wy/m)
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["Vy"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["p"] + k] * x[INDX["Vz"] + k])
             + (x[INDX["r"] + k] * x[INDX["Vx"] + k])
             - (grav * sin(x[INDX["phi"] + k]) * cos(x[INDX["theta"] + k]))
             - (f_wy / mass);
    nth += 1;
  }

  /*
  constraint for Vz(w)
  Vz(w)_dot = (q*u) - (p*v) + (g*cos(theta)*cos(phi)) + ((f_wz-f_t)/m)
  */
  for (Index k = 0 ; k <= N ; k++)
  {
    X[k] = x[INDX["Vz"] + k];
  }
  DX = multiply_M_V<Number, Index>(D, X, (N + 1));
  for (Index k = 0 ; k <= N ; k++)
  {
    g[nth] = (c1*DX[k])
             - (x[INDX["q"] + k] * x[INDX["Vx"] + k])
             + (x[INDX["p"] + k] * x[INDX["Vy"] + k])
             - (grav * cos(x[INDX["theta"] + k]) * cos(x[INDX["phi"] + k]))
             - ((f_wz - (x[INDX["netT"] + k] + f_t)) / mass);
    // - ((f_wz - f_t) / mass);
    nth += 1;
  }

  C.clear();
  D.clear();
  DX.clear();
  X.clear();

  //additional constraints initial values

  g[nth++] = x[INDX["p"    ] + N] - 0.00; /*intial value of "p"    */
  g[nth++] = x[INDX["q"    ] + N] - 0.00; /*intial value of "q"    */
  g[nth++] = x[INDX["r"    ] + N] - 0.00; /*intial value of "r"    */
  g[nth++] = x[INDX["phi"  ] + N] - 1.00; /*intial value of "phi"  */
  g[nth++] = x[INDX["theta"] + N] - 0.00; /*intial value of "theta"*/
  g[nth++] = x[INDX["psi"  ] + N] - 0.00; /*intial value of "psi"  */
  g[nth++] = x[INDX["z"    ] + N] - 1.00; /*intial value of "z"    */
  g[nth++] = x[INDX["Vz"   ] + N] - 0.00; /*intial value of "Vz"   */
  g[nth++] = x[INDX["y"    ] + N] - 0.00; /*intial value of "y"    */
  g[nth++] = x[INDX["Vy"   ] + N] - 0.00; /*intial value of "Vy"   */
  g[nth++] = x[INDX["x"    ] + N] - 0.00; /*intial value of "x"    */
  g[nth++] = x[INDX["Vx"   ] + N] - 0.00; /*intial value of "Vx"   */
  // g[nth++] = x[INDX["netT" ]+N] - 1.0; /*intial value of "netT" */
  // g[nth++] = x[INDX["Mx"   ]+N] - 1.0; /*intial value of "Mx"   */
  // g[nth++] = x[INDX["My"   ]+N] - 1.0; /*intial value of "My"   */
  // g[nth++] = x[INDX["Mz"   ]+N] - 1.0; /*intial value of "Mz"   */

  // additional constraints final values

  // g[nth++] = x[INDX["p"    ]] - 0.00; /*final value of "p"    */
  // g[nth++] = x[INDX["q"    ]] - 0.00; /*final value of "q"    */
  // g[nth++] = x[INDX["r"    ]] - 0.00; /*final value of "r"    */
  // g[nth++] = x[INDX["phi"  ]] - 0.00; /*final value of "phi"  */
  // g[nth++] = x[INDX["theta"]] - 0.00; /*final value of "theta"*/
  // g[nth++] = x[INDX["psi"  ]] - 0.00; /*final value of "psi"  */
  // g[nth++] = x[INDX["z"    ]] - 0.00; /*final value of "z"    */
  // g[nth++] = x[INDX["Vz"   ]] - 0.00; /*final value of "Vz"   */
  // g[nth++] = x[INDX["y"    ]] - 0.00; /*final value of "y"    */
  // g[nth++] = x[INDX["Vy"   ]] - 0.00; /*final value of "Vy"   */
  // g[nth++] = x[INDX["x"    ]] - 0.00; /*final value of "x"    */
  // g[nth++] = x[INDX["Vx"   ]] - 0.00; /*final value of "Vx"   */

  // g[nth++] = x[INDX["netT" ]]; /*final value of "netT" */
  // g[nth++] = x[INDX["Mx"   ]]; /*final value of "Mx"   */
  // g[nth++] = x[INDX["My"   ]]; /*final value of "My"   */
  // g[nth++] = x[INDX["Mz"   ]]; /*final value of "Mz"   */

 /*
   add some external force say from time stamps frm to to
   */
  /*Index frm = N;
  Index to = N - add_frc_const ;
  Number ext_frc = 2.5;
  for (Index k = frm ; k <=to ; k++)
  {
    g[nth] = (c1 * DX[k])
             - (x[INDX["r"] + k] * x[INDX["Vy"] + k])
             + (x[INDX["q"] + k] * x[INDX["Vz" ] + k])
             + (grav * sin(x[INDX["theta"] + k]))
             - (f_wx / mass)
             + ext_frc;
    nth += 1;
  }
  */
  //assert(nth == m);
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
  assert(n == sz_n);
  assert(m == sz_m);

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
  string state_var[n_s + n_c] = {"p", "q" , "r", "phi", "theta", "psi", "z", "Vz", "y", "Vy", "x", "Vx", "netT", "Mx", "My", "Mz"};
  std::vector<Number > time(N + 1);
  //Number tf = 1.25;
  for (Index i = 0; i <= N; i++)
  {
    time[i] = (t_f / 2.0) * (T[i] + 1.0);
    //std::cout << "t[" << i << "]" << " : " << time[i] << std::endl;
  }

  for (Index i = 0 ; i < n_s + n_c ; i++)
  {
    string file_name = "./Outputs/" + state_var[i] + ".txt";
    std::ofstream myfile;
    const char* fn = file_name.c_str();
    myfile.open(fn);
    //std::cout << state_var[i] << "\n";
    for (Index k = N; k >= 0; k--) {
      //std::cout << x[i] << std::endl;
      myfile << time[k] << "," << x[INDX[state_var[i]]  + k] << "\n";
    }
    //myfile.close();
  }

  //std::cout << x[n - 1] << std::endl;

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


  // std::ofstream mf1;
  // mf1.open("./Outputs/final_X.txt");
  // for (Index i = 0 ; i < n_s ; i++)
  // {
  //   //std::cout << state_var[i] << "\n";
  //   for (Index k = N; k >= 0; k--) {
  //     //std::cout << x[i] << std::endl;
  //     mf1 << x[INDX[state_var[i]]  + k] << " ";
  //   } mf1 << endl;
  // }
  // mf1.close();

  std::ofstream mf2;
  mf2.open("./Outputs/final_U.txt");
  mf2 << "[ ";
  for (Index i = n_s ; i < n_s + n_c ; i++)
  {
    //std::cout << state_var[i] << "\n";
    for (Index k = N; k >= 0; k--)
    {
      //std::cout << x[i] << std::endl;
      mf2 << x[INDX[state_var[i]]  + k] << " ";
    }
    if (i < n_s + n_c - 1)
    {
      mf2 << " ;" << endl;
    }
  }
  mf2 << "]";
  mf2.close();

  // std::cout << std::cout << std::endl
  //           << "Final value of the constraints:" << std::endl;
  // for (Index i = 0; i < m; i++) {
  //   std::cout << "g(" << i << ") = " << g[i] << std::endl;
  // }
  time.clear();
  // Analytical Solution
  /*std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;
  */
}
