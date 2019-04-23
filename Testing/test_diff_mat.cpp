// to test Differentiation Matrix over Chebyshev polynomial
//  g++ -I/home/dell/sliding_stone test_diff_mat.cpp -o test_diff_mat

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Differentiation_Matrix.cpp"
#define PI acos(-1)
#define f(val) sin(val)
using namespace std;
int main() {
	ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
	freopen("test_input.txt", "r", stdin);
	//freopen("test_output.txt", "w", stdout);
	int n;
	cin >> n;
	std::vector<double > X(n + 1), P(n + 1), C(n + 1), T(n + 1), t(n + 1);
	std::vector < std::vector<double > > D(n + 1 , std::vector <double > (n + 1));

	T = define_time_stamps<double, int>(n + 1);

	//cout << endl;
	// for (int i = 0 ; i <=n ; i++) {
	//cout <<double(i)<< ","<<T[i] << "\n";
	//} //cout << endl;

	//decimal t_0 = (-3.0*PI) , t_f = 3.0*PI ;
	decimal t_0 = -2.0, t_f = 2.0;

	for (integer i = 0; i <= n; i++)
	{
		t[i] = ((t_f - t_0) * (T[i]) + (t_f + t_0)) / 2.0;
		//myfile << "t[" << i << "]" << " : " << t[i] << std::endl;
	}

	std::ofstream myfile;
	myfile.open("test_output_1.txt");
	for (int i = 0 ; i <= n ; i++) {
		X[i] = f(t[i]);
		//X[i] = t[i];
		//X[i] =t[i]*t[i];
		myfile<< t[i] << "," << X[i] << "\n";
	}
	myfile.close();
	C = compute_c<double, int>(n + 1);
	// cout << endl;
	// for (int i = 0 ; i < n ; i++) {
	// 	cout << i << "," << C[i] << "\n";
	// } cout << endl;

	D = formulate_differentiation_matrix<double, int>(C, T, n + 1);
	std::cout.precision(4);
	/*
	for (int k = 0 ; k <=n; k++)
	{
		double sum=0;
		for (int j = 0 ; j <=n ; j++)
		{
			sum+=D[j][k];
			//cout << D[j][k] << " ";
		} cout <<sum<< endl<<endl;
	}
	*/
	P = multiply_D_X<double, int>(D, X, n + 1);
	// scale P
	/*
	decimal x_0 = f(t_0), x_f = f(t_f);
	for (integer i = 0; i <= n; i++)
	{
		P[i] = ((x_f - x_0) * (P[i]) + (x_f + x_0)) / 2.0;
		//std::cout << "t[" << i << "]" << " : " << t[i] << std::endl;
	}*/
	//std::ofstream myfile;
	myfile.open("test_output_2.txt");
	for (int i = 0 ; i <= n ; i++) {
		//cout << ((tf / 2.0) * (T[i] + 1)) << "," << ((pf / 2.0) * (P[i] + 1)) << "\n";
		//cout << ((tf / 2.0) * (T[i] + 1)) << "," << P[i] << "\n";
		myfile << t[i] << "," << (2.0/(t_f-t_0))*P[i] << "\n";
	}
	myfile.close();
	return 0;
}


