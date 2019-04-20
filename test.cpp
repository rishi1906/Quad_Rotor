#include <iostream>
#include<fstream>
#include <unordered_map>
#include <string>
#include <utility>
#include <set>
using namespace std;
#define nxt_indx (((++idx) * N) + idx)

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

unordered_map <string, int> INDX;

void get_indices(int N)
{
	int idx = 0;
	INDX.insert(make_pair("p"    , idx));
	INDX.insert(make_pair("q"    , nxt_indx));
	INDX.insert(make_pair("r"    , nxt_indx));
	INDX.insert(make_pair("phi"  , nxt_indx));
	INDX.insert(make_pair("theta", nxt_indx));
	INDX.insert(make_pair("psi"  , nxt_indx));
	INDX.insert(make_pair("z"    , nxt_indx));
	INDX.insert(make_pair("Vz"   , nxt_indx));
	INDX.insert(make_pair("y"    , nxt_indx));
	INDX.insert(make_pair("Vy"   , nxt_indx));
	INDX.insert(make_pair("x"   , nxt_indx));
	INDX.insert(make_pair("Vx"   , nxt_indx));
	INDX.insert(make_pair("netT" , nxt_indx));
	INDX.insert(make_pair("Mx"   , nxt_indx));
	INDX.insert(make_pair("My"   , nxt_indx));
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
	std::set<std::pair<std::string,int>, comp> set(INDX.begin(), INDX.end());

	for (auto const &pair: set) {
	   std::cout << '{' << pair.first << "," << pair.second << '}' << '\n';
	}//cout << "1";
	//cout << INDX["_p_"];

	return (0);
}