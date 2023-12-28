// 510C.cpp
#define PROBLEMSET
// Created by Ashish Negi

/********   All Required Header Files ********/
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <queue>
#include <deque>
#include <bitset>
#include <iterator>
#include <list>
#include <stack>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <utility>
#include <limits>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "bits/stdc++.h"
#include <numeric>		// contains inbuilt gcd(a, b) function
// #include "ext/pb_ds/assoc_container.hpp"
// #include "ext/pb_ds/tree_policy.hpp"

using namespace std;

/*******  All Required define Pre-Processors and typedef Constants *******/
const string nl = "\n";
#define scd(t) scanf("%d",&t)
#define scld(t) scanf("%ld",&t)
#define sclld(t) scanf("%lld",&t)
#define scc(t) scanf("%c",&t)
#define scs(t) scanf("%s",t)
#define scf(t) scanf("%f",&t)
#define sclf(t) scanf("%lf",&t)
#define inp(t) cin>>t
#define inpp(t,u) cin>>t>>u
#define inppp(t,u,v) cin>>t>>u>>v
#define out(t) cout<<t<<nl
#define outt(t,u) cout<<t<<" "<<u<<nl
#define outtt(t,u,v) cout<<t<<" "<<u<<" "<<v<<nl
#define outf(t, p) cout << fixed;	cout << setprecision(p);	out(t);
#define mems(a, b) memset(a, (b), sizeof(a))
#define lpj(i, j, k) for (int i=j ; i<k ; i+=1)
#define rlpj(i, j, k) for (int i=j ; i>=k ; i-=1)
#define lp(i, j) lpj(i, 0, j)
#define rlp(i, j) rlpj(i, j, 0)
#define inpv(a,n) lp(i, n) inp(a[i])
#define inpvv(a,n) lp(i, n)lp(j, n)	inp(a[i][j]);
#define outv(a,n) lp(i, n) {	if(i != 0 ){	cout << " ";}	cout << a[i];	}	cout<<nl;
#define outvv(a,n) lp(i, n){	lp(j, n){	if(j != 0 ){	cout << " ";}	cout << a[i][j];	}	cout<<nl;	}
#define outy() out("YES")
#define outn() out("NO")
#define all(cont) cont.begin(), cont.end()
#define rall(cont) cont.end(), cont.begin()
#define each(it, l) for (auto it = l.begin(); it != l.end(); it++)
#define IN(A, B, C) assert( B <= A && A <= C)
#define CNT(a, x)	count(all(a), x)
#define mpr make_pair
#define pbk push_back
#define INF (int)1e9
#define EPS 1e-9
#define PI 3.1415926535897932384626433832795
// #define MOD 1000000007
#define read(type) readInt<type>()
#define clk_start()	time_req = clock();
#define clk_end()	cout << "time taken to solve (in seconds) = "; outf((float)(clock() - time_req)/(float)CLOCKS_PER_SEC, 6);
#define rotl(a, i)	rotate(a.begin(), a.begin() + i, a.end())
#define rotr(a, i)	rotate(a.begin(), a.begin() + a.size() - i, a.end())
#define tern(cond, y, n)	((cond)?(y):(n))
#define ternyn(cond)	TERN(cond, outy(), outn())
const double pi = acos(-1.0);
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<unsigned long long int> vi64;
typedef vector<string> vs;
typedef vector<pii> vpii;
typedef vector<vi> vvi;
typedef map<int, int> mpii;
typedef set<int> seti;
typedef multiset<int> mseti;
typedef long int int32;
typedef unsigned long int uint32;
typedef long long int int64;
typedef unsigned long long int  uint64;

/****** Template of some basic operations *****/
template<typename T, typename U> inline void amin(T &x, U y) { if (y < x) x = y; }
template<typename T, typename U> inline void amax(T &x, U y) { if (x < y) x = y; }
/**********************************************/



/******* Debugging Class Template *******/
#ifdef DEBUG

#define debug(args...)     (Debugger()) , args
#define dbg(var1) cerr<<#var1<<" = "<<(var1)<<nl;
#define dbg2(var1,var2) cerr<<#var1<<" = "<<(var1)<<", "<<#var2<<" = "<<(var2)<<nl;
#define dbg3(var1,var2,var3) cerr<<#var1<<" = "<<(var1)<<", "<<#var2<<" = "<<(var2)<<", "<<#var3<<" = "<<(var3)<<nl;
#define dbg4(var1,var2,var3,var4) cerr<<#var1<<" = "<<(var1)<<", "<<#var2<<" = "<<(var2)<<", "<<#var3<<" = "<<(var3)<<", "<<#var4<<" = "<<(var4)<<nl;

class Debugger
{
public:
	Debugger(const std::string& _separator = " - ") :
		first(true), separator(_separator) {}

	template<typename ObjectType> Debugger& operator , (const ObjectType& v)
	{
		if (!first)
std: cerr << separator;
		std::cerr << v;
		first = false;
		return *this;
	}
~Debugger() {  std: cerr << nl;}

private:
	bool first;
	std::string separator;
};

#else
#define debug(args...)                  // Just strip off all debug tokens
#define dbg(args...)                  	// Just strip off all debug tokens
#define dbg2(args...)                  	// Just strip off all debug tokens
#define dbg3(args...)                  	// Just strip off all debug tokens
#define dbg4(args...)                  	// Just strip off all debug tokens
#endif

/************** Macros ****************/
#ifndef ONLINE_JUDGE
#define ONLINE_JUDGE
#endif	/*	ONLINE_JUDGE	*/
// #define SUBLIME_TEXT
// #define DEBUG
// #define CLOCK
// #define MULT_TC

/** Conditional variables/ constants **/
#ifdef CLOCK
clock_t time_req;
#endif /* CLOCK */

/***** Global variables/constants *****/
const int NMAX = 26;
int n, m;

/******* User-defined Functions *******/
vector<vector<int>> adj; // adjacency list of graph
vector<bool> visited;
vector<bool> inProcess; // To track vertices in the process of being explored
vector<int> ans;

bool hasCycle = false;

void dfs(int v) {
	visited[v] = true;
	inProcess[v] = true;

	for (int u : adj[v]) {
		if (!visited[u]) {
			dfs(u);
		} else if (inProcess[u]) {
			// Cycle detected
			hasCycle = true;
		}
	}

	inProcess[v] = false;
	ans.push_back(v);
}

void topological_sort() {
	visited.assign(NMAX, false);
	inProcess.assign(NMAX, false);
	ans.clear();
	hasCycle = false;

	for (int i = 0; i < NMAX; ++i) {
		if (!visited[i] && !hasCycle) {
			dfs(i);
		}
	}

	if (hasCycle) {
		cout << "Impossible" << endl;
		ans.clear(); // If a cycle is detected, clear the result vector
	} else {
		reverse(ans.begin(), ans.end());
	}
}

/**************************************/
void solve()
{
	adj.resize(26, vector<int>());
	// lp(i, 26) {
	// 	cout << i << nl;
	// 	lp(j, adj[i].size()) {
	// 		cout << adj[i][j] << " ";
	// 	}
	// 	cout << nl;
	// }
	inp(n);
	vector<string> s(n);
	lp(i, n) {
		inp(s[i]);
	}

	string last_string = s[0];

	lpj(i, 1, n) {
		dbg2(last_string, s[i]);
		if (s[i] > last_string) {
			int j = 0;
			while (j < last_string.length() && j < s[i].length() && last_string[j] == s[i][j]) {
				j++;
			}

			if (j == s[i].length()) {
			}
			else {
				dbg2(last_string[j], s[i][j]);
				if (find(adj[s[i][j] - 'a'].begin(), adj[s[i][j] - 'a'].end(), last_string[j] - 'a') != adj[s[i][j] - 'a'].end()) {
					out("Impossible");
					return;
				}


				adj[last_string[j] - 'a'].pbk(s[i][j] - 'a');
			}
		}
		else {
			int j = 0;
			while (j < last_string.length() && j < s[i].length() && last_string[j] == s[i][j]) {
				j++;
			}

			if (j < last_string.length() && j < s[i].length()) {
				if (find(adj[s[i][j] - 'a'].begin(), adj[s[i][j] - 'a'].end(), last_string[j] - 'a') != adj[s[i][j] - 'a'].end()) {
					out("Impossible");
					return;
				}

				dbg2(last_string[j], s[i][j]);
				adj[last_string[j] - 'a'].pbk(s[i][j] - 'a');
			}
			else {
				out("Impossible");
				return;
			}


		}
		last_string = s[i];
	}

	topological_sort();

	lp(i, ans.size()) {
		cout << (char)(ans[i] + (int)('a'));
	}
	cout << nl;

	// lp(i, 26) {
	// 	cout << (char)(i + (int)('a')) << ":" << nl;
	// 	lp(j, adj[i].size()) {
	// 		cout << (char)(adj[i][j] + (int)('a')) << " ";
	// 	}
	// 	cout << nl << nl;
	// }
	return;
}

/********** Main()  function **********/
int main()
{
#if !defined(ONLINE_JUDGE) || defined(SUBLIME_TEXT)
#if !defined(PROBLEMSET)
	freopen("../IO/input.txt", "r", stdin);
	freopen("../IO/output.txt", "w", stdout);
	freopen("../IO/log.txt", "w", stderr);
#else
	freopen("../../IO/input.txt", "r", stdin);
	freopen("../../IO/output.txt", "w", stdout);
	freopen("../../IO/log.txt", "w", stderr);
#endif
#endif

	std::ios::sync_with_stdio(false);
	cin.tie(0);

#ifdef MULT_TC
	int tc;
	inp(tc);

	while (tc--) {
#endif	/* MULT_TC */
#ifdef CLOCK
		clk_start();
#endif /* CLOCK */
		solve();
#ifdef CLOCK
		clk_end();
#endif /* CLOCK */
#ifdef MULT_TC
	}
#endif	/* MULT_TC */
	return 0;
}
/********  Main() Ends Here *************/
