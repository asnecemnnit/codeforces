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

/******** User-defined Function *******/

/*		check vector<int> is a palindrome		*/
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

/*		check string is a palindrome		*/
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

/*		check number is a power of 2		*/
bool isPowerOfTwo(int n)
{
	if(floor(log2(n)) == ceil(log2(n)))
		return true;
	else
		return false;
}

/*		count number of digits in an integer		*/
int countDigits(int n)
{
	int ans=0;

	while(n>0)
	{
		n /= 10;
		ans++;
	}

	return ans;
}

// Function to convert a map<key,value> to a multimap<value,key>
multimap<int, int> invert(map<int, int> & mymap)
{
	multimap<int, int> multiMap;

	map<int, int> :: iterator it;
	for (it=mymap.begin(); it!=mymap.end(); it++)
	{
		multiMap.insert(make_pair(it->second, it->first));
	}

	return multiMap;
}
/**************************************/
void solve()
{
	int n;
	INP(n);
	VI a(n);
	INPV(a,n);

	// sort(all(a));
	// int count = 0;
	// for(int i=0, j=a.size()-1; i<j; i++, j--)
	// {
	// 	if(a[i]==a[j])
	// 		break;
	// 	count += 2;
	// }
	//
	// OUT(a.size()-count);

	map<int, int> count;

	REP(i,n)
	{
		count[a[i]]++;
	}


	// map<int, int>::iterator itr;
	// for (itr = count.begin(); itr != count.end(); ++itr) {
    //     cout << '\t' << itr->first
    //          << '\t' << itr->second << '\n';
    // }
    // cout << endl;

	multimap<int, int> newmap = invert(count);
	// multimap<int, int>::iterator itr;
	// for (itr2 = newmap.begin(); itr != newmap.end(); ++itr) {
    //     cout << itr->second<<" occurs " << itr->first <<" times "<<endl;
    // }
    // cout << endl;

	int ans=n;
	while(newmap.size()>=2)
	{
		// cout<<"remaining size of multimap = "<<newmap.size()<<endl;
		// multimap<int, int>::iterator lmn;
		// for (lmn = newmap.begin(); lmn != newmap.end(); ++lmn) {
	    //     cout << lmn->second<<" occurs " << lmn->first <<" times "<<endl;
	    // }
	    // cout << endl;


		auto itr1 = newmap.rbegin();
		auto cnt1 = itr1->first;
		auto x1 = itr1->second;

		typedef multimap<int, int>::iterator iterator;
		pair<iterator, iterator> iterpair = newmap.equal_range(cnt1);

		// Erase specific pair
		iterator it = iterpair.first;
		for (; it != iterpair.second; ++it)
		{
		    if (it->second == x1) {
		        newmap.erase(it);
		        break;
		    }
		}


		auto itr2 = newmap.rbegin();
		auto cnt2 = itr2->first;
		auto x2 = itr2->second;
		iterpair = newmap.equal_range(cnt2);

		// Erase specific pair
		it = iterpair.first;
		for (; it != iterpair.second; ++it)
		{
		    if (it->second == x2) {
		        newmap.erase(it);
		        break;
		    }
		}




		cnt1--;
        cnt2--;
        ans -= 2;
		// cout<<x1<<" erased, count left : "<<cnt1<<endl;
		// cout<<x2<<" erased, count left : "<<cnt2<<endl;
		// cout<<"remaining size of multimap = "<<newmap.size()<<endl;

		if (cnt1) {
            newmap.insert((pair <int, int> (cnt1, x1)));
			// cout<<x1<<" inserted, count left : "<<cnt1<<endl;
        }

        if (cnt2) {
            newmap.insert((pair <int, int> (cnt2, x2)));
			// cout<<x2<<" inserted, count left : "<<cnt2<<endl;
        }

		// cout<<"remaining size of multimap = "<<newmap.size()<<endl;

	}

	OUT(ans);

}

/********** Main()  function **********/
int main()
{
    #define ONLINE_JUDGE
	#ifndef ONLINE_JUDGE
	freopen("input.txt","r",stdin);
	//freopen("output.txt","w",stdout);
	#endif

	std::ios::sync_with_stdio(false);

	int tc;
	INP(tc);

	while(tc--){
		solve();
	}
	return 0;
}
/********  Main() Ends Here *************/
