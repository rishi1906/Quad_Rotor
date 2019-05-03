#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <set>
using namespace std;
#define nxt_indx (((idx+1) * N) + (idx+1))

struct comp
{
	template<typename T>
	bool operator()(const T& l, const T& r) const
	{
		if (l.second != r.second)
			return l.second < r.second;

		return l.first < r.first;
	}
};

map <string, int> INDX;

void get_indices(int N)
{
	int idx = 0;
	INDX.insert(make_pair("p"    , idx));
	INDX.insert(make_pair("q"    , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("r"    , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("phi"  , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("theta", nxt_indx));
	idx += 1;
	INDX.insert(make_pair("psi"  , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("z"    , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("Vz"   , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("y"    , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("Vy"   , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("x"    , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("Vx"   , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("netT" , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("Mx"   , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("My"   , nxt_indx));
	idx += 1;
	INDX.insert(make_pair("Mz"   , nxt_indx));
};

int main() {

	// ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
	// freopen("test_input.txt", "r", stdin);
	// freopen("test_output.txt", "w", stdout);
	int N;
	cin >> N;
	get_indices(N);
	// create a empty vector of pairs
	//std::set<std::pair<std::string, int>, comp> set;
	//set.insert(INDX.begin(), INDX.end());

	/*for (auto const &pair : set) {
		std::cout << '{' << pair.first << "," << pair.second << '}' << '\n';
	}
	*/
	//cout << "1";
	//
	for (auto it = INDX.cbegin(); it != INDX.cend(); ++it)
	{
		std::cout << it->first << " " << it->second << "\n";
	}
	//cout << INDX["_p_"];
	// double max = 999.99;
	// double min = -999.99;
	// cout << (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
	// cout << (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
	// cout << (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;

	double A[12][12];
	std::ifstream f1;
	f1.open("../Inputs/A.txt");
	for (int i = 0 ; i < 12 ; i++)
	{
		for (int j = 0 ; j < 12 ; j++)
		{
			f1 >> A[i][j];

		} //cout << endl;
	}
	f1.close();
	for (int i = 0 ; i < 12 ; i++)
	{
		for (int j = 0 ; j < 12 ; j++)
		{
			cout << A[i][j] << " ";
		} cout << endl;
	}

	return (0);
}