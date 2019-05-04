// to test Integral Weight over Chebyshev Polynomial
// g++ -I/home/dell/Quad_Rotor test_integral.cpp -o test_integral
// ./test_integral
// python plot_test_graph.py
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Differentiation_Matrix.cpp"
#include "Integral_Weights.cpp"
#define PI acos(-1)
#define f(val) cos(val)
using namespace std;

int main()
{
	ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
	freopen("test_input.txt", "r", stdin);
	freopen("test_output.txt", "w", stdout);
	int n;

	cin >> n;
	std::vector<double > X(n + 1), T(n + 1), t(n + 1), w(n + 1);

	T = define_time_stamps<double, int>(n + 1);
	decimal t_0 = 0, t_f = 3;

	for (integer i = 0; i <= n; i++)
	{
		t[i] = ((t_f - t_0) * (T[i]) + (t_f + t_0)) / 2.0;
		//myfile << "t[" << i << "]" << " : " << t[i] << std::endl;
	}

	std::ofstream myfile;
	myfile.open("test_output_1.txt");
	for (int i = 0 ; i <= n ; i++) {
		X[i] = f(t[i]);
		//X[i] = t[i] * t[i] - 1;
		//X[i] = exp(t[i]);
		//X[i] = 1;
		//X[i] = t[i];
		//X[i] = t[i] * t[i];
		myfile << t[i] << "," << X[i] << "\n";
	}
	myfile.close();

	w = compute_integral_weights<double, int>(n + 1);

	myfile.open("test_output_2.txt");

	for (int i = 0 ; i <= n ; i++) {
		//cout << ((tf / 2.0) * (T[i] + 1)) << "," << ((pf / 2.0) * (P[i] + 1)) << "\n";
		//cout << ((tf / 2.0) * (T[i] + 1)) << "," << P[i] << "\n";
		myfile << t[i] << "," << ((t_f - t_0) / 2.0)*X[i]*w[i] << "\n";
	}
	myfile.close();

	double ans = 0.0;
	for (int i = 0 ; i <= n ; i++) {
		ans += X[i] * w[i];
	}
	cout << ((t_f - t_0) / 2.0)*ans;
	return 0;
}