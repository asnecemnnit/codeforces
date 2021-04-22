// p1.cpp
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
#define INPV(a,n) REP(i, n) INP(a[i])
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
bool found = false;
/******** User-defined Function *******/

/*		vector<int> is a palindrome		*/
bool isVecPalindrome(VI a)
{
	int n = a.size();
	if(n==1)
		return true;

	REP(i,n/2)
	{
		if(a[i] != a[n-1-i])
		{
			return false;
		}
	}

	return true;
}

/*		string is a palindrome		*/
bool isStrPalindrome(string a)
{
	int n = a.length();
	if(n==1)
		return true;

	REP(i,n/2)
	{
		if(a[i] != a[n-1-i])
		{
			return false;
		}
	}

	return true;
}

/*		number is a power of 2		*/
bool isPowerOfTwo(int n)
{
	if(floor(log2(n)) == ceil(log2(n)))
		return true;
	else
		return false;
}

// int countPeaks(VI a)
// {
// 	int n = a.size();
// 	if(n==1 or n==2)
// 		return 0;
// 	int count = 0;
// 	for(int i=1; i<n-1; i++)
// 	{
// 		if(a[i]>a[i-1] and a[i]>a[i+1])
// 		{
// 			count++;
// 			i++;
// 		}
// 	}
//
// 	return count;
// }
//
// // Generating permutation using Heap Algorithm
// VI heapPermutation(VI a, int size, int n, int k)
// {
// 	VI ans(n);
//     // if size becomes 1 then prints the obtained
//     // permutation
//     if (size == 1) {
//         if(countPeaks(a)==k)
// 		{
// 			// REP(i,n)
// 			// {
// 			// 	cout<<a[i]<<" ";
// 			// }
// 			// cout<<endl;
// 			found = true;
// 			return a;
// 		}
//     }
//
//     for (int i = 0; i < size; i++) {
//         ans = heapPermutation(a, size - 1, n, k);
// 		if(found)
// 			return ans;
//
//         // if size is odd, swap 0th i.e (first) and
//         // (size-1)th i.e (last) element
//         if (size % 2 == 1)
//             swap(a[0], a[size - 1]);
//
//         // If size is even, swap ith and
//         // (size-1)th i.e (last) element
//         else
//             swap(a[i], a[size - 1]);
//     }
//
// 	return ans;
// }
// /**************************************/


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
		int n,k;
		INP2(n,k);
		found = false;

		if(n%2==0)
		{
			if(k>=n/2)
			{
				OUT(-1);
				continue;
			}
		}
		else
		{
			if(k>n/2)
			{
				OUT(-1);
				continue;
			}
		}

		VI p(n);

		REP(i,n)
		{
			p[i] = i+1;
		}
		// VI ans = heapPermutation(p, n, n, k);
		for(int i=1; k>0; i+=2)
		{
			int temp = p[i];
			p[i] = p[i+1];
			p[i+1] = temp;
			k--;
		}
		REP(i,n)
		{
			cout<<p[i]<<" ";
		}
		cout<<endl;


	}
	return 0;
}
/********  Main() Ends Here *************/
