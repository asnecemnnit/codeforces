// Created by Ashish Negi

/************************* All Required Header Files *************************/
#include "bits/stdc++.h"
#include "ext/pb_ds/assoc_container.hpp"
#include "ext/pb_ds/tree_policy.hpp"
#include <algorithm>
#include <assert.h>
#include <bitset>
#include <deque>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <numeric>
#include <numeric> // contains inbuilt gcd(a, b) function
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <utility>
#include <vector>

using namespace std;
using namespace __gnu_pbds;

/********** All Required Define Pre-Processors and Typedef Constants *********/
const string nl = "\n";
#define scd(t) scanf("%d", &t)
#define scld(t) scanf("%ld", &t)
#define sclld(t) scanf("%lld", &t)
#define scc(t) scanf("%c", &t)
#define scs(t) scanf("%s", t)
#define scf(t) scanf("%f", &t)
#define sclf(t) scanf("%lf", &t)
#define inp(t) cin >> t
#define inpp(t, u) cin >> t >> u
#define inppp(t, u, v) cin >> t >> u >> v
#define out(t) cout << t << nl
#define outt(t, u) cout << t << " " << u << nl
#define outtt(t, u, v) cout << t << " " << u << " " << v << nl
#define outf(t, p)                                                             \
  cout << fixed;                                                               \
  cout << setprecision(p);                                                     \
  out(t);
#define mems(a, b) memset(a, (b), sizeof(a))
#define lpj(i, j, k) for (int i = j; i < k; i += 1)
#define rlpj(i, j, k) for (int i = j; i >= k; i -= 1)
#define lp(i, j) lpj(i, 0, j)
#define rlp(i, j) rlpj(i, j, 0)
#define inpv(a, n) lp(i, n) inp(a[i])
#define inpvv(a, n) lp(i, n) lp(j, n) inp(a[i][j]);
#define outv(a, n)                                                             \
  lp(i, n) {                                                                   \
    if (i != 0) {                                                              \
      cout << " ";                                                             \
    }                                                                          \
    cout << a[i];                                                              \
  }                                                                            \
  cout << nl;
#define outvv(a, n)                                                            \
  lp(i, n) {                                                                   \
    lp(j, n) {                                                                 \
      if (j != 0) {                                                            \
        cout << " ";                                                           \
      }                                                                        \
      cout << a[i][j];                                                         \
    }                                                                          \
    cout << nl;                                                                \
  }
#define outy() out("YES")
#define outn() out("NO")
#define all(cont) cont.begin(), cont.end()
#define rall(cont) cont.end(), cont.begin()
#define each(it, l) for (auto it = l.begin(); it != l.end(); it++)
#define IN(A, B, C) assert(B <= A && A <= C)
#define CNT(a, x) count(all(a), x)
#define mpr make_pair
#define pbk push_back
#define INF (int)1e9
#define EPS 1e-9
#define PI 3.1415926535897932384626433832795
// #define MOD 1000000007
#define read(type) readInt<type>()
#define clk_start() time_req = clock();
#define clk_end()                                                              \
  cout << "time taken to solve (in seconds) = ";                               \
  outf((float)(clock() - time_req) / (float)CLOCKS_PER_SEC, 6);
#define rotl(a, i) rotate(a.begin(), a.begin() + i, a.end())
#define rotr(a, i) rotate(a.begin(), a.begin() + a.size() - i, a.end())
#define tern(cond, y, n) ((cond) ? (y) : (n))
#define ternyn(cond) TERN(cond, outy(), outn())
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
typedef unsigned long long int uint64;

/****** Template of some basic operations *****/
template <typename T, typename U> inline void amin(T &x, U y) {
  if (y < x)
    x = y;
}
template <typename T, typename U> inline void amax(T &x, U y) {
  if (x < y)
    x = y;
}

/****************************** Miscellaneous ********************************/
// Custom hash to be used in conjunction with unordered_map, unordered_set
// Read blog post : Blowing up unordered_map, and how to stop getting hacked on
// it (https://codeforces.com/blog/entry/62393) Usage: unordered_map<int, int,
// custom_hash>, unordered_set<int, custom_hash>
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM =
        chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};

// Ordered set (to be used for ordered sets (multiset))
// Usage:
// ordered_set ost;
// ost.insert(x)
// ost.order_of_key(x)
typedef tree<int, null_type, less<int>, rb_tree_tag,
             tree_order_statistics_node_update>
    ordered_set;

/******************** Template of Fast I/O Methods ***************************/
template <typename T> inline void write(T x) {
  int i = 20;
  char buf[21];
  // buf[10] = 0;
  buf[20] = '\n';

  do {
    buf[--i] = x % 10 + '0';
    x /= 10;
  } while (x);
  do {
    putchar(buf[i]);
  } while (buf[i++] != '\n');
}
template <typename T> inline T readInt() {
  T n = 0, s = 1;
  char p = getchar();
  if (p == '-')
    s = -1;
  while ((p < '0' || p > '9') && p != EOF && p != '-')
    p = getchar();
  if (p == '-')
    s = -1, p = getchar();
  while (p >= '0' && p <= '9') {
    n = (n << 3) + (n << 1) + (p - '0');
    p = getchar();
  }

  return n * s;
}

/************************************ RNG ************************************/
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
inline int64_t random_long(long long l = LLONG_MIN, long long r = LLONG_MAX) {
  uniform_int_distribution<int64_t> generator(l, r);
  return generator(rng);
}

/************************* Debugging Class Template **************************/
#ifdef DEBUG

#define debug(args...) (Debugger()), args
#define dbg(var1) cerr << #var1 << " = " << (var1) << nl;
#define dbg2(var1, var2)                                                       \
  cerr << #var1 << " = " << (var1) << ", " << #var2 << " = " << (var2) << nl;
#define dbg3(var1, var2, var3)                                                 \
  cerr << #var1 << " = " << (var1) << ", " << #var2 << " = " << (var2) << ", " \
       << #var3 << " = " << (var3) << nl;
#define dbg4(var1, var2, var3, var4)                                           \
  cerr << #var1 << " = " << (var1) << ", " << #var2 << " = " << (var2) << ", " \
       << #var3 << " = " << (var3) << ", " << #var4 << " = " << (var4) << nl;

class Debugger {
public:
  Debugger(const std::string &_separator = " - ")
      : first(true), separator(_separator) {}

  template <typename ObjectType> Debugger &operator,(const ObjectType &v) {
    if (!first)
    std:
      cerr << separator;
    std::cerr << v;
    first = false;
    return *this;
  }
  ~Debugger() {
  std:
    cerr << nl;
  }

private:
  bool first;
  std::string separator;
};

#else
#define debug(args...) // Just strip off all debug tokens
#define dbg(args...)   // Just strip off all debug tokens
#define dbg2(args...)  // Just strip off all debug tokens
#define dbg3(args...)  // Just strip off all debug tokens
#define dbg4(args...)  // Just strip off all debug tokens
#endif

/********************************** Macros ***********************************/
#ifndef ONLINE_JUDGE
#define ONLINE_JUDGE
#endif /*	ONLINE_JUDGE	*/
// #define SUBLIME_TEXT
// #define DEBUG
// #define CLOCK
#define MULT_TC

/********************** Conditional variables/ constants *********************/
#ifdef CLOCK
clock_t time_req;
#endif /* CLOCK */

/************************* Global Variables/Constants ************************/
const int NMAX = 3e5;
int n, m;

/************************** User-Defined Functions ***************************/

/******************************** Solve **************************************/
void solve() { return; }

/****************************  Main() Starts Here ****************************/
int main() {
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
#endif /* MULT_TC */
#ifdef CLOCK
    clk_start();
#endif /* CLOCK */
    solve();
#ifdef CLOCK
    clk_end();
#endif /* CLOCK */
#ifdef MULT_TC
  }
#endif /* MULT_TC */
  return 0;
}
/***************************** Main() Ends Here *****************************/
