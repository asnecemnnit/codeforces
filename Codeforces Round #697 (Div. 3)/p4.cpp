// p4.cpp
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

using namespace std;

/*******  All Required define Pre-Processors and typedef Constants *******/
#define SCD(t) scanf("%d",&t)
#define SCLD(t) scanf("%ld",&t)
#define SCLLD(t) scanf("%lld",&t)
#define SCC(t) scanf("%c",&t)
#define SCS(t) scanf("%s",t)
#define SCF(t) scanf("%f",&t)
#define SCLF(t) scanf("%lf",&t)
#define INP(t) cin>>t
#define INP2(t,u) cin>>t>>u
#define INP3(t,u,v) cin>>t>>u>>v
#define OUT(t) cout<<t<<endl
#define OUT2(t,u) cout<<t<<" "<<u<<endl
#define OUT3(t,u,v) cout<<t<<" "<<u<<" "<<v<<endl
#define MEM(a, b) memset(a, (b), sizeof(a))
#define FOR(i, j, k, in) for (int i=j ; i<k ; i+=in)
#define RFOR(i, j, k, in) for (int i=j ; i>=k ; i-=in)
#define REP(i, j) FOR(i, 0, j, 1)
#define RREP(i, j) RFOR(i, j, 0, 1)
#define all(cont) cont.begin(), cont.end()
#define rall(cont) cont.end(), cont.begin()
#define FOREACH(it, l) for (auto it = l.begin(); it != l.end(); it++)
#define IN(A, B, C) assert( B <= A && A <= C)
#define MP make_pair
#define PB push_back
#define INF (int)1e9
#define EPS 1e-9
#define PI 3.1415926535897932384626433832795
#define MOD 1000000007
#define read(type) readInt<type>()
const double pi=acos(-1.0);
typedef pair<int, int> PII;
typedef vector<int> VI;
typedef vector<unsigned long long int> VI64;
typedef vector<string> VS;
typedef vector<PII> VII;
typedef vector<VI> VVI;
typedef map<int,int> MPII;
typedef set<int> SETI;
typedef multiset<int> MSETI;
typedef long int int32;
typedef unsigned long int uint32;
typedef long long int int64;
typedef unsigned long long int  uint64;

/****** Template of some basic operations *****/
template<typename T, typename U> inline void amin(T &x, U y) { if(y < x) x = y; }
template<typename T, typename U> inline void amax(T &x, U y) { if(x < y) x = y; }
/**********************************************/

/****** Template of Fast I/O Methods *********/
template <typename T> inline void write(T x)
{
	int i = 20;
	char buf[21];
	// buf[10] = 0;
	buf[20] = '\n';

	do
	{
		buf[--i] = x % 10 + '0';
		x/= 10;
	}while(x);
	do
	{
		putchar(buf[i]);
	} while (buf[i++] != '\n');
}

template <typename T> inline T readInt()
{
	T n=0,s=1;
	char p=getchar();
	if(p=='-')
		s=-1;
	while((p<'0'||p>'9')&&p!=EOF&&p!='-')
		p=getchar();
	if(p=='-')
		s=-1,p=getchar();
	while(p>='0'&&p<='9') {
		n = (n<< 3) + (n<< 1) + (p - '0');
		p=getchar();
	}

	return n*s;
}
/************************************/


/******* Debugging Class Template *******/
#define DEBUG

#ifdef DEBUG

    #define debug(args...)     (Debugger()) , args

    class Debugger
    {
        public:
        Debugger(const std::string& _separator = " - ") :
        first(true), separator(_separator){}

        template<typename ObjectType> Debugger& operator , (const ObjectType& v)
        {
            if(!first)
                std:cerr << separator;
            std::cerr << v;
            first = false;
            return *this;
        }
        ~Debugger() {  std:cerr << endl;}

        private:
        bool first;
        std::string separator;
    };

#else
    #define debug(args...)                  // Just strip off all debug tokens
#endif

/**************************************/

/******** User-defined Function *******/
// // Returns the value of minimum points
// int knapSackRec(int W, VI wt,
//                 VI val, int i,
//                 int** dp)
// {
//     // base condition
//     if (i < 0)
//         return INT_MAX;
//     if (dp[i][W] != INT_MAX)
//         return dp[i][W];
//
//     if (wt[i] > W) {
//
//         // Store the value of function call
//         // stack in table before return
//         dp[i][W] = knapSackRec(W, wt,
//                                val, i - 1,
//                                dp);
//         return dp[i][W];
//     }
//     else {
//         // Store value in a table before return
//         dp[i][W] = min(val[i]
//                       + knapSackRec(W - wt[i],
//                                    wt, val,
//                                    i - 1, dp),
//                        knapSackRec(W, wt, val,
//                                    i - 1, dp));
//
//         // Return value of table after storing
//         return dp[i][W];
//     }
// }
//
// // Returns the minimum number of convenience points that Polycarp will lose
// int knapSack(int W, VI wt, VI val, int n)
// {
// 	// double pointer to declare the
//     // table dynamically
//     int** dp;
//     dp = new int*[n];
//
//     // loop to create the table dynamically
//     for (int i = 0; i < n; i++)
//         dp[i] = new int[W + 1];
//
//     // loop to initially filled the
//     // table with -1
//     for (int i = 0; i < n; i++)
//         for (int j = 0; j < W + 1; j++)
//             dp[i][j] = INT_MAX;
//     return knapSackRec(W, wt, val, n - 1, dp);
// }

// Returns the minimum value that
// can be put in a knapsack of capacity W
int knapSack(int W, VI wt, VI val, int n)
{
    int i, w;
    VVI K(n + 1, VI(W + 1, INT_MAX));

    // Build table K[][] in bottom up manner
    for(i = 1; i <= n; i++)
    {
        for(w = 0; w <= W; w++)
        {
			// K[i][w] = K[i - 1][w];

			if(w <= wt[i - 1])
				K[i][w] = K[i - 1][w];

			else
				K[i][w] = min(val[i - 1] +
							K[i - 1][w - wt[i - 1]],
							K[i - 1][w]);


        }
    }
    return K[n][W];
}

/**************************************/


/********** Main()  function **********/
int main()
{
    #define ONLINE_JUDGE
	#ifndef ONLINE_JUDGE
	freopen("input.txt","r",stdin);
	//freopen("output.txt","w",stdout);
	#endif


	int tc;
	INP(tc);

	while(tc--){
		int n, m;
		INP2(n,m);

		VI a(n);
		VI b(n);

		int memory_occupied = 0;

		REP(i,n)
		{
			INP(a[i]);
			memory_occupied += a[i];
		}

		REP(i,n)
		{
			INP(b[i]);
		}

		if(memory_occupied < m)
			OUT(-1);
		else
		{
			OUT(knapSack(m, a, b, n));
		}
	}
	return 0;
}
/********  Main() Ends Here *************/
