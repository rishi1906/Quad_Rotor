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

std::unordered_map <string, int> get_indices(int N)
{
	int idx = 0;
	std::unordered_map <string, int> mp;
	mp.insert(make_pair("p"    , idx));
	mp.insert(make_pair("q"    , nxt_indx_N));
	mp.insert(make_pair("r"    , nxt_indx_N));
	mp.insert(make_pair("phi"  , nxt_indx_N));
	mp.insert(make_pair("theta", nxt_indx_N));
	mp.insert(make_pair("psi"  , nxt_indx_N));
	mp.insert(make_pair("z"    , nxt_indx_N));
	mp.insert(make_pair("Vz"   , nxt_indx_N));
	mp.insert(make_pair("y"    , nxt_indx_N));
	mp.insert(make_pair("Vy"   , nxt_indx_N));
	mp.insert(make_pair("x"    , nxt_indx_N));
	mp.insert(make_pair("Vx"   , nxt_indx_N));
	idx = ((12 * N) + 12);
	mp.insert(make_pair("netT" , nxt_indx_4));
	mp.insert(make_pair("Mx"   , nxt_indx_4));
	mp.insert(make_pair("My"   , nxt_indx_4));
	mp.insert(make_pair("Mz"   , nxt_indx_4));
	return mp;
};