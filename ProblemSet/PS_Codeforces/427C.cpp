// 427C.cpp
#define PROBLEMSET
// Created by Ashish Negi

/************************* All Required Header Files *************************/
#include "bits/stdc++.h"
#include "ext/pb_ds/assoc_container.hpp"
#include "ext/pb_ds/tree_policy.hpp"

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
typedef long double uld64;

/********************* Template of some basic operations *********************/
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
// Template for efficiently writing integers to the output stream.
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

// Template for efficiently reading integers from the input stream.
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
// Initialize a 64-bit Mersenne Twister pseudo-random number generator using the
// current time as seed
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

// Template for generating a random 64-bit integer in the specified range [l, r]
inline int64_t random_long(long long l = LLONG_MIN, long long r = LLONG_MAX) {
  // Create a uniform distribution for generating random integers in the
  // specified range
  uniform_int_distribution<int64_t> generator(l, r);

  // Generate and return a random 64-bit integer
  return generator(rng);
}

/**************************** Modular Arithmetic *****************************/
const int mod = 1e9 + 7;

template <const int MOD> struct modular_int {
  int x;
  static vector<int> inverse_list;
  const static int inverse_limit;
  const static bool is_prime;
  modular_int() { x = 0; }
  template <typename T> modular_int(const T z) {
    x = (z % MOD);
    if (x < 0)
      x += MOD;
  }
  modular_int(const modular_int<MOD> *z) { x = z->x; }
  modular_int(const modular_int<MOD> &z) { x = z.x; }
  modular_int operator-(const modular_int<MOD> &m) const {
    modular_int<MOD> U;
    U.x = x - m.x;
    if (U.x < 0)
      U.x += MOD;
    return U;
  }
  modular_int operator+(const modular_int<MOD> &m) const {
    modular_int<MOD> U;
    U.x = x + m.x;
    if (U.x >= MOD)
      U.x -= MOD;
    return U;
  }
  modular_int &operator-=(const modular_int<MOD> &m) {
    x -= m.x;
    if (x < 0)
      x += MOD;
    return *this;
  }
  modular_int &operator+=(const modular_int<MOD> &m) {
    x += m.x;
    if (x >= MOD)
      x -= MOD;
    return *this;
  }
  modular_int operator*(const modular_int<MOD> &m) const {
    modular_int<MOD> U;
    U.x = (x * 1ull * m.x) % MOD;
    return U;
  }
  modular_int &operator*=(const modular_int<MOD> &m) {
    x = (x * 1ull * m.x) % MOD;
    return *this;
  }
  template <typename T>
  friend modular_int operator+(const T &l, const modular_int<MOD> &p) {
    return (modular_int<MOD>(l) + p);
  }
  template <typename T>
  friend modular_int operator-(const T &l, const modular_int<MOD> &p) {
    return (modular_int<MOD>(l) - p);
  }
  template <typename T>
  friend modular_int operator*(const T &l, const modular_int<MOD> &p) {
    return (modular_int<MOD>(l) * p);
  }
  template <typename T>
  friend modular_int operator/(const T &l, const modular_int<MOD> &p) {
    return (modular_int<MOD>(l) / p);
  }
  int value() const { return x; }
  modular_int operator^(const modular_int<MOD> &cpower) const {
    modular_int<MOD> ans;
    ans.x = 1;
    modular_int<MOD> curr(this);
    int power = cpower.x;

    if (curr.x <= 1) {
      if (power == 0)
        curr.x = 1;
      return curr;
    }
    while (power > 0) {
      if (power & 1) {
        ans *= curr;
      }
      power >>= 1;
      if (power)
        curr *= curr;
    }
    return ans;
  }
  modular_int operator^(long long power) const {
    modular_int<MOD> ans;
    ans.x = 1;
    modular_int<MOD> curr(this);
    if (curr.x <= 1) {
      if (power == 0)
        curr.x = 1;
      return curr;
    }
    // Prime Mods
    if (power >= MOD && is_prime) {
      power %= (MOD - 1);
    }

    while (power > 0) {
      if (power & 1) {
        ans *= curr;
      }
      power >>= 1;
      if (power)
        curr *= curr;
    }
    return ans;
  }
  modular_int operator^(int power) const {
    modular_int<MOD> ans;
    ans.x = 1;
    modular_int<MOD> curr(this);

    if (curr.x <= 1) {
      if (power == 0)
        curr.x = 1;
      return curr;
    }
    while (power > 0) {
      if (power & 1) {
        ans *= curr;
      }
      power >>= 1;
      if (power)
        curr *= curr;
    }
    return ans;
  }
  template <typename T> modular_int &operator^=(T power) {
    modular_int<MOD> res = (*this) ^ power;
    x = res.x;
    return *this;
  }
  template <typename T> modular_int pow(T x) { return (*this) ^ x; }
  pair<long long, long long> gcd(const int a, const int b) const {
    if (b == 0)
      return {1, 0};
    pair<long long, long long> c = gcd(b, a % b);
    return {c.second, c.first - (a / b) * c.second};
  }
  inline void init_inverse_list() const {
    vector<int> &dp = modular_int<MOD>::inverse_list;
    dp.resize(modular_int<MOD>::inverse_limit + 1);
    int n = modular_int<MOD>::inverse_limit;
    dp[0] = 1;
    if (n > 1)
      dp[1] = 1;
    for (int i = 2; i <= n; ++i) {
      dp[i] = (dp[MOD % i] * 1ull * (MOD - MOD / i)) % MOD;
    }
  }
  modular_int<MOD> get_inv() const {
    if (modular_int<MOD>::inverse_list.size() <
        modular_int<MOD>::inverse_limit + 1)
      init_inverse_list();

    if (this->x <= modular_int<MOD>::inverse_limit) {
      return modular_int<MOD>::inverse_list[this->x];
    }
    pair<long long, long long> G = gcd(this->x, MOD);
    return modular_int<MOD>(G.first);
  }
  modular_int<MOD> operator-() const {
    modular_int<MOD> v(0);
    v -= (*this);
    return v;
  }
  modular_int operator/(const modular_int<MOD> &m) const {
    modular_int<MOD> U(this);
    U *= m.get_inv();
    return U;
  }
  modular_int &operator/=(const modular_int<MOD> &m) {
    (*this) *= m.get_inv();
    return *this;
  }
  bool operator==(const modular_int<MOD> &m) const { return x == m.x; }
  bool operator<(const modular_int<MOD> &m) const { return x < m.x; }
  template <typename T> bool operator==(const T &m) const {
    return (*this) == (modular_int<MOD>(m));
  }
  template <typename T> bool operator<(const T &m) const {
    return x < (modular_int<MOD>(m)).x;
  }
  template <typename T> bool operator>(const T &m) const {
    return x > (modular_int<MOD>(m)).x;
  }
  template <typename T>
  friend bool operator==(const T &x, const modular_int<MOD> &m) {
    return (modular_int<MOD>(x)).x == m.x;
  }
  template <typename T>
  friend bool operator<(const T &x, const modular_int<MOD> &m) {
    return (modular_int<MOD>(x)).x < m.x;
  }
  template <typename T>
  friend bool operator>(const T &x, const modular_int<MOD> &m) {
    return (modular_int<MOD>(x)).x > m.x;
  }
  friend istream &operator>>(istream &is, modular_int<MOD> &p) {
    int64_t val;
    is >> val;
    p.x = (val % MOD);
    if (p.x < 0)
      p.x += MOD;
    return is;
  }
  friend ostream &operator<<(ostream &os, const modular_int<MOD> &p) {
    return os << p.x;
  }
};

using mint = modular_int<mod>;
template <const int MOD> vector<int> modular_int<MOD>::inverse_list;
template <const int MOD> const int modular_int<MOD>::inverse_limit = -1;
template <const int MOD> const bool modular_int<MOD>::is_prime = true;
template <> //-> useful if computing inverse fact = i_f[i-1] / i;
const int modular_int<mod>::inverse_limit = 1000;

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
// #define VSCODE
// #define DEBUG
// #define CLOCK
// #define MULT_TC

/********************** Conditional variables/ constants *********************/
#ifdef CLOCK
clock_t time_req;
#endif /* CLOCK */

/************************* Global Variables/Constants ************************/
const int NMAX = 3e5 + 1;
int n, m;
vector<bool> visited;
vi pos, st;
vvi res;
vector<vector<int>> adj;

/************************** User-Defined Functions ***************************/
int dfs(int u) {
  int low = st.size();
  pos[u] = st.size();
  st.push_back(u);
  for (int v : adj[u]) {
    if (not visited[v]) {
      low = min(low, ~pos[v] ? pos[v] : dfs(v));
    }
  }
  if (low == pos[u]) {
    res.emplace_back(st.begin() + low, st.end());
    for (int v : res.back()) {
      visited[v] = true;
    }
    st.resize(low);
  }
  return low;
};

vvi scc_fun(vector<vector<int>> adj) {

  pos.resize(n);
  lp(i, n) { pos[i] = -1; }
  for (int i = 0; i < n; i += 1) {
    if (not visited[i]) {
      dfs(i);
    }
  }
  reverse(all(res));
  return res;
}
/******************************** Solve **************************************/
void solve() {
  mint num_ways = 1;
  uint64 min_money = 0;

  inp(n);

  vector<uint64> cost(n);
  inpv(cost, n);

  adj.resize(n);
  visited.resize(n);

  inp(m);
  lp(i, m) {
    int u, v;
    inpp(u, v);
    u--;
    v--;
    adj[u].pbk(v);
  }

  vvi res = scc_fun(adj);

  lp(i, res.size()) {
    uint64 min_money_j = 1e9 + 1;
    lp(j, res[i].size()) { min_money_j = min(min_money_j, cost[res[i][j]]); }
    mint num_ways_j = 0;
    lp(j, res[i].size()) {
      num_ways_j += ((min_money_j == cost[res[i][j]]) ? 1 : 0);
    }
    min_money += min_money_j;
    num_ways *= num_ways_j;
  }

  outt(min_money, num_ways);

  return;
}

/****************************  Main() Starts Here ****************************/
int main() {
#if !defined(ONLINE_JUDGE) || defined(SUBLIME_TEXT) || defined(VSCODE)
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
