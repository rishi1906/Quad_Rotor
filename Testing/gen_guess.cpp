#include <iostream>
#include <vector>
#include "../Differentiation_Matrix.cpp"
using namespace std;
const double N = 41;
int main() {

	ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
	// freopen("test_input.txt", "r", stdin);
	freopen("../Inputs/guess_z.txt", "w", stdout);
	//double t_0 = 0.0, t_f = 1.0;
	double z[41];
	z[40] = 0.0;
	z[0] = 1.0;
	cout << z[0] << "\n";
	std::vector<double > t;
	t = define_time_stamps<double, int>(N);
	for (int i = 1 ; i < N ; i++)
	{
		//z(t)=z(0)+(z(tf)-z(0))/(tf-t0)*(t-t0)
		z[i] = z[0] + ((z[40] - z[0]) * (t[i] - t[0])) / (t[N - 1] - t[0]);
		cout << z[i] << "\n";
	}

	return 0;
}