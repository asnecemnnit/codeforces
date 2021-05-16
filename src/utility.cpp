// Created by Ashish Negi

/******** User-defined Function *******/

////////////////////////////////////////////////////////////////////////////////

/*	precomputes factorial till	M	*/
#define M 100000
mint fact[M + 1], i_f[M + 1];

void precomputeFactorial() {
	fact[0] = 1;
	i_f[0] = 1;

	for (int i = 1; i < M; ++i) {
		fact[i] = fact[i - 1] * i;
		i_f[i] = i_f[i - 1] / i;
	}
}

// Combinations
/*	Computes ncr using fact and i_f	*/
mint ncr(int n, int r) {
	if (n < r || r < 0 || n < 0) {
		return 0;
	}
	return fact[n] * i_f[r] * i_f[n - r];
}

/*	Computes ncr from scratch */
mint Ncr(int n, int r) {
	mint ans = 1;
	for (int i = 1; i <= r; ++i) {
		ans *= n;
		ans /= i;
		--n;
	}
	return ans;
}

// Permutations
/*	Computes npr using fact and i_f	*/
mint npr(int n, int r) {
	if (n < r || r < 0 || n < 0) {
		return 0;
	}
	return fact[n] * i_f[n - r];
}

/*	Computes npr from scratch */
mint Npr(int n, int r) {
	mint ans = 1;
	for (int i = 1; i <= r; ++i) {
		ans *= n;
		// ans /= i;
		--n;
	}
	return ans;
}
////////////////////////////////////////////////////////////////////////////////
/*  returns gcd of a & b    */
int gcd(int a, int b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);

}
////////////////////////////////////////////////////////////////////////////////

/*	check vector<int> is a palindrome	*/
bool isVecPalindrome(VI a)
{
	int n = a.size();
	if (n == 1)
		return true;

	REP(i, n / 2)
	{
		if (a[i] != a[n - 1 - i])
		{
			return false;
		}
	}

	return true;
}

/*	check string is a palindrome	*/
bool isStrPalindrome(string a)
{
	int n = a.length();
	if (n == 1)
		return true;

	REP(i, n / 2)
	{
		if (a[i] != a[n - 1 - i])
		{
			return false;
		}
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////

/*	check number is a power of 2	*/
bool isPowerOfTwo(int n)
{
	if (floor(log2(n)) == ceil(log2(n)))
		return true;
	else
		return false;
}

/*	count number of digits in an integer	*/
int countDigits(int n)
{
	int ans = 0;

	while (n > 0)
	{
		n /= 10;
		ans++;
	}

	return ans;
}

////////////////////////////////////////////////////////////////////////////////
/* Count Inversions in an array O(n^2). Can be optimized using merge sort.
 * Inversion : i<j and a[i]>a[j] */
int getInvCount(vector<int> arr, int n)
{
	int inv_count = 0;
	for (int i = 0; i < n - 1; i++)
		for (int j = i + 1; j < n; j++)
			if (arr[i] > arr[j])
				inv_count++;

	return inv_count;
}

/* Minimum Adjacent Swaps needed to get another target permutation string s2 of string s1 */
int getMinAdjSwaps(string s1, string s2)
{
	vector<vector<int>> Q(26);
	vector<int> P(s1.size());

	for ( int i = 0; i < s1.size(); ++i )
		Q[s2[i] - 'a'].push_back(i); // basically, Q is a vector [0 .. 25] of lists

	vector<int> temp(26, 0);
	for ( int i = 0; i < s1.size(); ++i )
		P[i] = 1 + Q[s1[i] - 'a'][ temp[s1[i] - 'a']++ ];



	return getInvCount(P, s1.size());
}


////////////////////////////////////////////////////////////////////////////////

/* Returns prime bool vector using sieve of Eratosthenes algorithm */
vector<bool> SieveOfEratosthenes(int n = 1000000)
{
	// Create a boolean array
	// "prime[0..n]" and initialize
	// all entries it as true.
	// A value in prime[i] will
	// finally be false if i is
	// Not a prime, else true.
	vector<bool> prime(n + 1);
	memset(prime, true, sizeof(prime));

	for (int p = 2; p * p <= n; p++)
	{
		// If prime[p] is not changed,
		// then it is a prime
		if (prime[p] == true)
		{
			// Update all multiples
			// of p greater than or
			// equal to the square of it
			// numbers which are multiple
			// of p and are less than p^2
			// are already been marked.
			for (int i = p * p; i <= n; i += p)
				prime[i] = false;
		}
	}

	// Print all prime numbers
	// for (int p = 2; p <= n; p++)
	//     if (prime[p])
	//         cout << p << " ";

	return prime;
}

// isPrime[] : isPrime[i] is true if number is prime
// prime[] : stores all prime number less than N
// SPF[] that store smallest prime factor of number
// [for Exp : smallest prime factor of '8' and '16'
// is '2' so we put SPF[8] = 2 , SPF[16] = 2 ]
const long long MAX_SIZE = 1000001;
vector<long long >isprime(MAX_SIZE , true);
vector<long long >prime;
vector<long long >SPF(MAX_SIZE);
/* function generate all prime number less then N in O(n) */
void manipulated_seive(int N)
{
	// 0 and 1 are not prime
	isprime[0] = isprime[1] = false ;

	// Fill rest of the entries
	for (long long int i = 2; i < N ; i++)
	{
		// If isPrime[i] == True then i is
		// prime number
		if (isprime[i])
		{
			// put i into prime[] vector
			prime.push_back(i);

			// A prime number is its own smallest
			// prime factor
			SPF[i] = i;
		}

		// Remove all multiples of  i*prime[j] which are
		// not prime by making isPrime[i*prime[j]] = false
		// and put smallest prime factor of i*Prime[j] as prime[j]
		// [ for exp :let  i = 5 , j = 0 , prime[j] = 2 [ i*prime[j] = 10 ]
		// so smallest prime factor of '10' is '2' that is prime[j] ]
		// this loop run only one time for number which are not prime
		for (long long int j = 0;
		        j < (int)prime.size() &&
		        i * prime[j] < N && prime[j] <= SPF[i];
		        j++)
		{
			isprime[i * prime[j]] = false;

			// put smallest prime factor of i*prime[j]
			SPF[i * prime[j]] = prime[j] ;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/*	Returns the length of the longest
    increasing subsequence (LIS) in nums. O(nlogn)	*/
int lengthOfLIS(vector<int>& nums)
{
	int n = nums.size();
	vector<int> pilesTop;
	pilesTop.reserve(n);
	for (int i = 0; i < n; i++)
	{
		// O(logn)
		auto left = lower_bound(pilesTop.begin(), pilesTop.end(), nums[i]);
		if (left == pilesTop.end())
			pilesTop.push_back(nums[i]);
		else
			*left = nums[i];

		// cout<<"i = "<<i<<", *left "<<*left<<endl;
		// for(auto x:pilesTop)
		// {
		//     cout<<x<<" ";
		// }
		// cout<<endl;
	}
	return pilesTop.size();
}

/////////////////////////////////////////////////////////////////////////////

/*	Returns the length of the longest
    increasing subsequence (LIS) in nums. O(nlogn)	*/

// Binary search (note boundaries in the caller)
int CeilIndex(vector<int>& v, int l, int r, int key)
{
	while (r - l > 1) {
		int m = l + (r - l) / 2;
		if (v[m] >= key)
			r = m;
		else
			l = m;
	}

	return r;
}

int LongestIncreasingSubsequenceLength(vector<int>& v)
{
	if (v.size() == 0)
		return 0;

	vector<int> tail(v.size(), 0);
	int length = 1; // always points empty slot in tail

	tail[0] = v[0];
	for (int i = 1; i < v.size(); i++) {

		// new smallest value
		if (v[i] < tail[0])
			tail[0] = v[i];

		// v[i] extends largest subsequence
		else if (v[i] > tail[length - 1])
			tail[length++] = v[i];

		// v[i] will become end candidate of an existing
		// subsequence or Throw away larger elements in all
		// LIS, to make room for upcoming grater elements
		// than v[i] (and also, v[i] would have already
		// appeared in one of LIS, identify the location
		// and replace it)
		else
			tail[CeilIndex(tail, -1, length - 1, v[i])] = v[i];
	}

	return length;
}

////////////////////////////////////////////////////////////////////////////////

/*	Returns the length of the longest
    increasing subsequence (LIS) in nums. O(n^2)	*/
int lengthOfLIS(vector<int>& nums)
{
	int n = nums.size();
	vector<int> lis(n);

	lis[0] = 1;

	/* Compute optimized LIS values in
	   bottom up manner */
	for (int i = 1; i < n; i++ )
	{
		lis[i] = 1;
		for (int j = 0; j < i; j++ )
			if (nums[i] > nums[j] && lis[i] < lis[j] + 1)
				lis[i] = lis[j] + 1;
	}

	// Return maximum value in lis[]
	return *max_element(lis.begin(), lis.end());
}

////////////////////////////////////////////////////////////////////////////////

/* Returns length of longest common subsequence for X[0..m-1], Y[0..n-1] */
int lcs(string X, string Y, int m, int n )
{
	vector<int> L(m + 1, vector<int>(n + 1));
	int i, j;

	/* Following steps build L[m+1][n+1] in
	   bottom up fashion. Note that L[i][j]
	   contains length of LCS of X[0..i-1]
	   and Y[0..j-1] */
	for (i = 0; i <= m; i++)
	{
		for (j = 0; j <= n; j++)
		{
			if (i == 0 || j == 0)
				L[i][j] = 0;

			else if (X[i - 1] == Y[j - 1])
				L[i][j] = L[i - 1][j - 1] + 1;

			else
				L[i][j] = max(L[i - 1][j], L[i][j - 1]);
		}
	}

	/* L[m][n] contains length of LCS
	for X[0..n-1] and Y[0..m-1] */
	return L[m][n];
}

////////////////////////////////////////////////////////////////////////////////

/* Returns length of longest
   common substring of X[0..m-1]
   and Y[0..n-1] */
int LCSubStr(string X, string Y, int m, int n)
{
	// Create a table to store
	// lengths of longest
	// common suffixes of substrings.
	// Note that LCSuff[i][j] contains
	// length of longest common suffix
	// of X[0..i-1] and Y[0..j-1].

	vector<int> LCSuff(m + 1, vector<int>(n + 1));
	int result = 0; // To store length of the
	// longest common substring

	/* Following steps build LCSuff[m+1][n+1] in
	    bottom up fashion. */
	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			// The first row and first column
			// entries have no logical meaning,
			// they are used only for simplicity
			// of program
			if (i == 0 || j == 0)
				LCSuff[i][j] = 0;

			else if (X[i - 1] == Y[j - 1]) {
				LCSuff[i][j] = LCSuff[i - 1][j - 1] + 1;
				result = max(result, LCSuff[i][j]);
			}
			else
				LCSuff[i][j] = 0;
		}
	}
	return result;
}

////////////////////////////////////////////////////////////////////////////////

/*	Find number of times a string occurs as
	a subsequence in given string	*/

// A Dynamic Programming based C++ program to find the
// number of times the second string occurs in the first
// string, whether continuous or discontinuous

// Iterative DP function to find the number of times
// the second string occurs in the first string,
// whether continuous or discontinuous
int count(string a, string b)
{
	int m = a.length();
	int n = b.length();

	// Create a table to store results of sub-problems
	vector<vector<int>> lookup(m + 1, vector<int>(n + 1, 0));

	// If first string is empty
	for (int i = 0; i <= n; ++i)
		lookup[0][i] = 0;

	// If second string is empty
	for (int i = 0; i <= m; ++i)
		lookup[i][0] = 1;

	// Fill lookup[][] in bottom up manner
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			// If last characters are same, we have two
			// options -
			// 1. consider last characters of both strings
			//    in solution
			// 2. ignore last character of first string
			if (a[i - 1] == b[j - 1])
				lookup[i][j] = lookup[i - 1][j - 1] +
				               lookup[i - 1][j];

			else
				// If last character are different, ignore
				// last character of first string
				lookup[i][j] = lookup[i - 1][j];
		}
	}

	return lookup[m][n];
}

///////////////////////////////////////////////////////////////////////
// Returns the maximum value that
// can be put in a knapsack of capacity W
// Time Complexity: O(N*W)
// Auxiliary Space: O(N*W)
// DP top-down memoization

// Returns the value of maximum profit
int knapSackRec(int W, vector<int> wt,
                vector<int> val, int i,
                int** dp)
{
	// base condition
	if (i < 0)
		return 0;
	if (dp[i][W] != -1)
		return dp[i][W];

	if (wt[i] > W) {

		// Store the value of function call
		// stack in table before return
		dp[i][W] = knapSackRec(W, wt,
		                       val, i - 1,
		                       dp);
		return dp[i][W];
	}
	else {
		// Store value in a table before return
		dp[i][W] = max(val[i]
		               + knapSackRec(W - wt[i],
		                             wt, val,
		                             i - 1, dp),
		               knapSackRec(W, wt, val,
		                           i - 1, dp));

		// Return value of table after storing
		return dp[i][W];
	}
}

int knapSack(int W, vector<int> wt, vector<int> val, int n)
{
	// double pointer to declare the
	// table dynamically
	int** dp;
	dp = new int*[n];

	// loop to create the table dynamically
	for (int i = 0; i < n; i++)
		dp[i] = new int[W + 1];

	// loop to initially filled the
	// table with -1
	for (int i = 0; i < n; i++)
		for (int j = 0; j < W + 1; j++)
			dp[i][j] = -1;
	return knapSackRec(W, wt, val, n - 1, dp);
}
//////////////////////////////////////////////////////////////////////////
// Returns the maximum value that
// can be put in a knapsack of capacity W
// Time Complexity: O(N*W)
// Auxiliary Space: O(N*W)
// DP bottom-up
int knapSack(int W, vector<int> wt, vector<int> val, int n)
{
	int i, w;
	vector<vector<int>> K(n + 1, vector<int>(W + 1));

	// Build table K[][] in bottom up manner
	for (i = 0; i <= n; i++)
	{
		for (w = 0; w <= W; w++)
		{
			if (i == 0 || w == 0)
				K[i][w] = 0;
			else if (wt[i - 1] <= w)
				K[i][w] = max(val[i - 1] +
				              K[i - 1][w - wt[i - 1]],
				              K[i - 1][w]);
			else
				K[i][w] = K[i - 1][w];
		}
	}
	return K[n][W];
}
//////////////////////////////////////////////////////////////////////////

/*  Unbounded Knapsack (Repetition of items allowed)    */
// Returns the maximum value with knapsack of
// W capacity
int unboundedKnapsack(int W, int n,
                      vector<int> val, vector<int> wt)
{
	// dp[i] is going to store maximum value
	// with knapsack capacity i.
	vector<int> dp(W + 1, 0);

	// Fill dp[] using above recursive formula
	for (int i = 0; i <= W; i++)
		for (int j = 0; j < n; j++)
			if (wt[j] <= i)
				dp[i] = max(dp[i], dp[i - wt[j]] + val[j]);

	return dp[W];
}

/////////////////////////////////////////////////////////////////////////////
/*  In Fractional Knapsack, we can break items for maximizing the total value of knapsack   */

// Structure for an item which stores weight and
// corresponding value of Item
struct Item {
	int value, weight;

	// Constructor
	Item(int value, int weight)
	{
		this->value = value;
		this->weight = weight;
	}
};

// Comparison function to sort Item according to val/weight
// ratio
bool cmp(struct Item a, struct Item b)
{
	double r1 = (double)a.value / (double)a.weight;
	double r2 = (double)b.value / (double)b.weight;
	return r1 > r2;
}

// Main greedy function to solve problem
double fractionalKnapsack(int W, struct Item arr[], int n)
{
	//    sorting Item on basis of ratio
	sort(arr, arr + n, cmp);

	//    Uncomment to see new order of Items with their
	//    ratio
	/*
	for (int i = 0; i < n; i++)
	{
	    cout << arr[i].value << "  " << arr[i].weight << " :
	"
	         << ((double)arr[i].value / arr[i].weight) <<
	endl;
	}
	*/

	int curWeight = 0; // Current weight in knapsack
	double finalvalue = 0.0; // Result (value in Knapsack)

	// Looping through all Items
	for (int i = 0; i < n; i++) {
		// If adding Item won't overflow, add it completely
		if (curWeight + arr[i].weight <= W) {
			curWeight += arr[i].weight;
			finalvalue += arr[i].value;
		}

		// If we can't add current Item, add fractional part
		// of it
		else {
			int remain = W - curWeight;
			finalvalue += arr[i].value
			              * ((double)remain
			                 / (double)arr[i].weight);
			break;
		}
	}

	// Returning final value
	return finalvalue;
}

////////////////////////////////////////////////////////////////////////////////

/*	At each iteration the vertex v is selected
	which has the smallest distance d[v] among
	all the unmarked vertices. If the distance
	to selected vertex v is equal to infinity,
	the algorithm stops. Otherwise the vertex
	is marked, and all the edges going out from
	this vertex are checked. If relaxation along
	the edge is possible (i.e. distance d[to]can
	be improved), the distance d[to] and predecessor
	p[to] are updated	*/
const int INF = 1000000000;
vector<vector<pair<int, int>>> adj;

/*	Dijkstra Algorithm for finding shortest paths from
	source to all vertices in the given graph	*/
void dijkstra(int s, vector<int> & d, vector<int> & p) {
	int n = adj.size();
	d.assign(n, INF);
	p.assign(n, -1);
	vector<bool> u(n, false);

	d[s] = 0;
	for (int i = 0; i < n; i++) {
		int v = -1;
		for (int j = 0; j < n; j++) {
			if (!u[j] && (v == -1 || d[j] < d[v]))
				v = j;
		}

		if (d[v] == INF)
			break;

		u[v] = true;
		for (auto edge : adj[v]) {
			int to = edge.first;
			int len = edge.second;

			if (d[v] + len < d[to]) {
				d[to] = d[v] + len;
				p[to] = v;
			}
		}
	}
}

/*	p[] stores the predecessors of all vertices
	(except starting vertex s). The path to any vertex t
	can be restored in the following way	*/
vector<int> restore_path(int s, int t, vector<int> const& p) {
	vector<int> path;

	for (int v = t; v != s; v = p[v])
		path.push_back(v);
	path.push_back(s);

	reverse(path.begin(), path.end());
	return path;
}


// Number of vertices in the graph
#define V 9

// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(int dist[], bool sptSet[])
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < V; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

// A utility function to print the constructed distance array
void printSolution(int dist[])
{
	printf("Vertex \t\t Distance from Source\n");
	for (int i = 0; i < V; i++)
		printf("%d \t\t %d\n", i, dist[i]);
}

// Function that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
void dijkstra(int graph[V][V], int src)
{
	int dist[V]; // The output array.  dist[i] will hold the shortest
	// distance from src to i

	bool sptSet[V]; // sptSet[i] will be true if vertex i is included in shortest
	// path tree or shortest distance from src to i is finalized

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < V; i++)
		dist[i] = INT_MAX, sptSet[i] = false;

	// Distance of source vertex from itself is always 0
	dist[src] = 0;

	// Find shortest path for all vertices
	for (int count = 0; count < V - 1; count++) {
		// Pick the minimum distance vertex from the set of vertices not
		// yet processed. u is always equal to src in the first iteration.
		int u = minDistance(dist, sptSet);

		// Mark the picked vertex as processed
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the picked vertex.
		for (int v = 0; v < V; v++)

			// Update dist[v] only if is not in sptSet, there is an edge from
			// u to v, and total weight of path from src to  v through u is
			// smaller than current value of dist[v]
			if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
			        && dist[u] + graph[u][v] < dist[v])
				dist[v] = dist[u] + graph[u][v];
	}

	// print the constructed distance array
	printSolution(dist);
}

////////////////////////////////////////////////////////////////////////////////

// a structure to represent a weighted edge in graph
struct Edge {
	int src, dest, weight;
};

// a structure to represent a connected, directed and
// weighted graph
struct Graph {
	// V-> Number of vertices, E-> Number of edges
	int V, E;

	// graph is represented as an array of edges.
	struct Edge* edge;
};

// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
	struct Graph* graph = new Graph;
	graph->V = V;
	graph->E = E;
	graph->edge = new Edge[E];
	return graph;
}

// A utility function used to print the solution
void printArr(int dist[], int n)
{
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < n; ++i)
		printf("%d \t\t %d\n", i, dist[i]);
}

// The main function that finds shortest distances from src to
// all other vertices using Bellman-Ford algorithm.  The function
// also detects negative weight cycle
void BellmanFord(struct Graph* graph, int src)
{
	int V = graph->V;
	int E = graph->E;
	int dist[V];

	// Step 1: Initialize distances from src to all other vertices
	// as INFINITE
	for (int i = 0; i < V; i++)
		dist[i] = INT_MAX;
	dist[src] = 0;

	// Step 2: Relax all edges |V| - 1 times. A simple shortest
	// path from src to any other vertex can have at-most |V| - 1
	// edges
	for (int i = 1; i <= V - 1; i++) {
		for (int j = 0; j < E; j++) {
			int u = graph->edge[j].src;
			int v = graph->edge[j].dest;
			int weight = graph->edge[j].weight;
			if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
				dist[v] = dist[u] + weight;
		}
	}

	// Step 3: check for negative-weight cycles.  The above step
	// guarantees shortest distances if graph doesn't contain
	// negative weight cycle.  If we get a shorter path, then there
	// is a cycle.
	for (int i = 0; i < E; i++) {
		int u = graph->edge[i].src;
		int v = graph->edge[i].dest;
		int weight = graph->edge[i].weight;
		if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
			printf("Graph contains negative weight cycle");
			return; // If negative cycle is detected, simply return
		}
	}

	printArr(dist, V);

	return;
}

////////////////////////////////////////////////////////////////////////////////

/* 	find a spanning tree of this graph which connects
 	all vertices and has the least weight (i.e. the sum
	of weights of edges is minimal). A spanning tree is
	a set of edges such that any vertex can reach any
	other by exactly one simple path. The spanning tree
	with the least weight is called a minimum spanning tree.

	The minimum spanning tree is built gradually by adding
	edges one at a time. At first the spanning tree consists
	only of a single vertex (chosen arbitrarily). Then the
	minimum weight edge outgoing from this vertex is selected
	and added to the spanning tree. After that the spanning
	tree already consists of two vertices. Now select and
	add the edge with the minimum weight that has one end
	in an already selected vertex (i.e. a vertex that is
	already part of the spanning tree), and the other end
	in an unselected vertex. And so on, i.e. every time we
	select and add the edge with minimal weight that connects
	one selected vertex with one unselected vertex. The
	process is repeated until the spanning tree contains
	all vertices (or equivalently until we have n−1 edges)	*/

// Dense graphs : O(n2) time and O(n) memory
vector<vector<pair<int, int>>> adj;
vector<pair<int, pair<int, int>>> v;
vector<pair<int, int>> vals;

/* 	find a spanning tree of this graph which connects
 	all vertices and has the least weight (i.e. the sum
	of weights of edges is minimal). A spanning tree is
	a set of edges such that any vertex can reach any
	other by exactly one simple path. The spanning tree
	with the least weight is called a minimum spanning tree.

	The minimum spanning tree is built gradually by adding
	edges one at a time. At first the spanning tree consists
	only of a single vertex (chosen arbitrarily). Then the
	minimum weight edge outgoing from this vertex is selected
	and added to the spanning tree. After that the spanning
	tree already consists of two vertices. Now select and
	add the edge with the minimum weight that has one end
	in an already selected vertex (i.e. a vertex that is
	already part of the spanning tree), and the other end
	in an unselected vertex. And so on, i.e. every time we
	select and add the edge with minimal weight that connects
	one selected vertex with one unselected vertex. The
	process is repeated until the spanning tree contains
	all vertices (or equivalently until we have n−1 edges)	*/

// Dense graphs : O(n2) time and O(n) memory

void prim(int n) {
	int total_weight = 0;
	vector<bool> selected(n, false);

	// first = to, second = w
	vector<pair<int, int>> min_e(n, make_pair(-1, INF));
	min_e[0].second = 0;

	for (int i = 0; i < n; ++i) {
		int v = -1;
		for (int j = 0; j < n; ++j) {
			if (!selected[j] && (v == -1 || min_e[j].second < min_e[v].second))
				v = j;
		}

		// INF means there is no edge
		if (min_e[v].second == INF) {
			cout << "No MST!" << endl;
			exit(0);
		}

		selected[v] = true;
		total_weight += min_e[v].second;
		// if (min_e[v].first != -1)
		//     cout << v << " " << min_e[v].first << endl;

		for (int j = 0; j < adj[v].size(); ++j) {
			if (adj[v][j].second < min_e[adj[v][j].first].second)
				min_e[adj[v][j].first] = make_pair(v, adj[v][j].second);
		}

		// for (int i=0; i<n; ++i) {
		// 	debug("i",i,"to",min_e[i].first,"weight",min_e[i].second);
		// }
	}

	cout << total_weight << endl;

	return;
}

// Sparse graphs : O(mlogn) time
const int INF = 1000000000;
struct Edge {
	int w = INF, to = -1;
	bool operator<(Edge const& other) const {
		return make_pair(w, to) < make_pair(other.w, other.to);
	}
};
int n;
vector<vector<Edge>> adj;

void prim() {
	int total_weight = 0;
	vector<Edge> min_e(n);
	min_e[0].w = 0;
	set<Edge> q;
	q.insert({0, 0});
	vector<bool> selected(n, false);
	for (int i = 0; i < n; ++i) {
		if (q.empty()) {
			cout << "No MST!" << endl;
			exit(0);
		}

		int v = q.begin()->to;
		selected[v] = true;
		total_weight += q.begin()->w;
		q.erase(q.begin());

		if (min_e[v].to != -1)
			cout << v << " " << min_e[v].to << endl;

		for (Edge e : adj[v]) {
			if (!selected[e.to] && e.w < min_e[e.to].w) {
				q.erase({min_e[e.to].w, e.to});
				min_e[e.to] = {e.w, v};
				q.insert({e.w, e.to});
			}
		}
	}

	cout << total_weight << endl;
}

////////////////////////////////////////////////////////////////////////////////
/* Finding MST using Kruskal's algorithm

	1. Sort all the edges in non-decreasing order of their weight.
	2. Pick the smallest edge. Check if it forms a cycle with the
	spanning tree formed so far. If cycle is not formed, include
	this edge. Else, discard it.
	3. Repeat step#2 until there are (V-1) edges in the spanning tree. */

// a structure to represent a
// weighted edge in graph
class Edge {
public:
	int src, dest, weight;
};

// a structure to represent a connected,
// undirected and weighted graph
class Graph {
public:

	// V-> Number of vertices, E-> Number of edges
	int V, E;

	// graph is represented as an array of edges.
	// Since the graph is undirected, the edge
	// from src to dest is also edge from dest
	// to src. Both are counted as 1 edge here.
	Edge* edge;
};

// Creates a graph with V vertices and E edges
Graph* createGraph(int V, int E)
{
	Graph* graph = new Graph;
	graph->V = V;
	graph->E = E;

	graph->edge = new Edge[E];

	return graph;
}

// A structure to represent a subset for union-find
class subset {
public:
	int parent;
	int rank;
};

// A utility function to find set of an element i
// (uses path compression technique)
int find(subset subsets[], int i)
{
	// find root and make root as parent of i
	// (path compression)
	if (subsets[i].parent != i)
		subsets[i].parent
		    = find(subsets, subsets[i].parent);

	return subsets[i].parent;
}

// A function that does union of two sets of x and y
// (uses union by rank)
void Union(subset subsets[], int x, int y)
{
	int xroot = find(subsets, x);
	int yroot = find(subsets, y);

	// Attach smaller rank tree under root of high
	// rank tree (Union by Rank)
	if (subsets[xroot].rank < subsets[yroot].rank)
		subsets[xroot].parent = yroot;
	else if (subsets[xroot].rank > subsets[yroot].rank)
		subsets[yroot].parent = xroot;

	// If ranks are same, then make one as root and
	// increment its rank by one
	else {
		subsets[yroot].parent = xroot;
		subsets[xroot].rank++;
	}
}

// Compare two edges according to their weights.
// Used in qsort() for sorting an array of edges
int myComp(const void* a, const void* b)
{
	Edge* a1 = (Edge*)a;
	Edge* b1 = (Edge*)b;
	return a1->weight > b1->weight;
}

// The main function to construct MST using Kruskal's
// algorithm
// Greedy approach
// Step #2 uses the Union-Find algorithm to detect cycles
void KruskalMST(Graph* graph)
{
	int V = graph->V;
	Edge result[V]; // Tnis will store the resultant MST
	int e = 0; // An index variable, used for result[]
	int i = 0; // An index variable, used for sorted edges

	// Step 1: Sort all the edges in non-decreasing
	// order of their weight. If we are not allowed to
	// change the given graph, we can create a copy of
	// array of edges
	qsort(graph->edge, graph->E, sizeof(graph->edge[0]),
	      myComp);

	// Allocate memory for creating V ssubsets
	subset* subsets = new subset[(V * sizeof(subset))];

	// Create V subsets with single elements
	for (int v = 0; v < V; ++v)
	{
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}

	// Number of edges to be taken is equal to V-1
	while (e < V - 1 && i < graph->E)
	{
		// Step 2: Pick the smallest edge. And increment
		// the index for next iteration
		Edge next_edge = graph->edge[i++];

		int x = find(subsets, next_edge.src);
		int y = find(subsets, next_edge.dest);

		// If including this edge does't cause cycle,
		// include it in result and increment the index
		// of result for next edge
		if (x != y) {
			result[e++] = next_edge;
			Union(subsets, x, y);
		}
		// Else discard the next_edge
	}

	// print the contents of result[] to display the
	// built MST
	cout << "Following are the edges in the constructed "
	     "MST\n";
	int minimumCost = 0;
	for (i = 0; i < e; ++i)
	{
		cout << result[i].src << " -- " << result[i].dest
		     << " == " << result[i].weight << endl;
		minimumCost = minimumCost + result[i].weight;
	}
	// return;
	cout << "Minimum Cost Spanning Tree: " << minimumCost
	     << endl;
}

// Kruskal Algorithm (STL approach)
// Step #2 uses the Union-Find algorithm to detect cycles
// Creating shortcut for an integer pair
typedef  pair<int, int> iPair;

// Structure to represent a graph
struct Graph
{
	int V, E;
	vector< pair<int, iPair> > edges;

	// Constructor
	Graph(int V, int E)
	{
		this->V = V;
		this->E = E;
	}

	// Utility function to add an edge
	void addEdge(int u, int v, int w)
	{
		edges.push_back({w, {u, v}});
	}

	// Function to find MST using Kruskal's
	// MST algorithm
	int kruskalMST();
};

// To represent Disjoint Sets
struct DisjointSets
{
	int *parent, *rnk;
	int n;

	// Constructor.
	DisjointSets(int n)
	{
		// Allocate memory
		this->n = n;
		parent = new int[n + 1];
		rnk = new int[n + 1];

		// Initially, all vertices are in
		// different sets and have rank 0.
		for (int i = 0; i <= n; i++)
		{
			rnk[i] = 0;

			//every element is parent of itself
			parent[i] = i;
		}
	}

	// Find the parent of a node 'u'
	// Path Compression
	int find(int u)
	{
		/* Make the parent of the nodes in the path
		   from u--> parent[u] point to parent[u] */
		if (u != parent[u])
			parent[u] = find(parent[u]);
		return parent[u];
	}

	// Union by rank
	void merge(int x, int y)
	{
		x = find(x), y = find(y);

		/* Make tree with smaller height
		   a subtree of the other tree  */
		if (rnk[x] > rnk[y])
			parent[y] = x;
		else // If rnk[x] <= rnk[y]
			parent[x] = y;

		if (rnk[x] == rnk[y])
			rnk[y]++;
	}
};

/* Functions returns weight of the MST*/

int Graph::kruskalMST()
{
	int mst_wt = 0; // Initialize result

	// Sort edges in increasing order on basis of cost
	sort(edges.begin(), edges.end());

	// Create disjoint sets
	DisjointSets ds(V);

	// Iterate through all sorted edges
	vector< pair<int, iPair> >::iterator it;
	for (it = edges.begin(); it != edges.end(); it++)
	{
		int u = it->second.first;
		int v = it->second.second;

		int set_u = ds.find(u);
		int set_v = ds.find(v);

		// Check if the selected edge is creating
		// a cycle or not (Cycle is created if u
		// and v belong to same set)
		if (set_u != set_v)
		{
			// Current edge will be in the MST
			// so print it
			cout << u << " - " << v << endl;

			// Update MST weight
			mst_wt += it->first;

			// Merge two sets
			ds.merge(set_u, set_v);
		}
	}

	return mst_wt;
}

////////////////////////////////////////////////////////////////////////////////

/*	KMP Algorithm for pattern matching	*/

// Prints occurrences of txt[] in pat[]
void KMPSearch(string pat, string txt)
{
	int M = pat.length();
	int N = txt.length();

	// create lps[] that will hold the longest prefix suffix
	// values for pattern
	int lps[M];
	vector<int> lps(M);

	// Preprocess the pattern (calculate lps[] array)
	computeLPSArray(pat, M, lps);

	int i = 0; // index for txt[]
	int j = 0; // index for pat[]
	while (i < N) {
		if (pat[j] == txt[i]) {
			j++;
			i++;
		}

		if (j == M) {
			printf("Found pattern at index %d ", i - j);
			j = lps[j - 1];
		}

		// mismatch after j matches
		else if (i < N && pat[j] != txt[i]) {
			// Do not match lps[0..lps[j-1]] characters,
			// they will match anyway
			if (j != 0)
				j = lps[j - 1];
			else
				i = i + 1;
		}
	}
}

// Fills lps[] for given patttern pat[0..M-1]
void computeLPSArray(string pat, int M, vector<int> &lps)
{
	// length of the previous longest prefix suffix
	int len = 0;

	lps[0] = 0; // lps[0] is always 0

	// the loop calculates lps[i] for i = 1 to M-1
	int i = 1;
	while (i < M) {
		if (pat[i] == pat[len]) {
			len++;
			lps[i] = len;
			i++;
		}
		else // (pat[i] != pat[len])
		{
			// This is tricky. Consider the example.
			// AAACAAAA and i = 7. The idea is similar
			// to search step.
			if (len != 0) {
				len = lps[len - 1];

				// Also, note that we do not increment
				// i here
			}
			else // if (len == 0)
			{
				lps[i] = 0;
				i++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

vector<vector<pair<int, int>>> adj;
vector<pair<int, pair<int, int>>> v;
vector<pair<int, int>> vals;

/*	Add an edge (weighted, if isWeighted = true) from u to v in a graph	*/
void addEdge(int u,
             int v, bool isWeighted , int wt)
{
	if (isWeighted)
	{
		adj[u].push_back(make_pair(v, wt));
		adj[v].push_back(make_pair(u, wt));
	}
	else
	{
		// ignore wt
		adj[u].push_back(make_pair(v, 0));
		adj[v].push_back(make_pair(v, 0));
	}

	return;
}

/*	Create an undirected graph of n vertices and
	returns vector of adj vectors	*/
void createGraph(int n, bool isWeighted)
{

	for (int i = 0; i < v.size(); i++)
	{
		addEdge(v[i].second.first, v[i].second.second,
		        isWeighted, v[i].first);
	}

	return;

}


////////////////////////////////////////////////////////////////////////////////
/*	MO’s Algorithm (Query Square Root Decomposition)
	The preprocessing part takes O(m Log m) time.
	Processing all queries takes
	O(n * √n) + O(m * √n) = O((m+n) * √n) time.	*/

// Variable to represent block size. This is made global
// so compare() of sort can use it.
int block;

// Structure to represent a query range
struct Query
{
	int L, R;
};

// Function used to sort all queries so that all queries
// of the same block are arranged together and within a block,
// queries are sorted in increasing order of R values.
bool compare(Query x, Query y)
{
	// Different blocks, sort by block.
	if (x.L / block != y.L / block)
		return x.L / block < y.L / block;

	// Same block, sort by R value
	return x.R < y.R;
}

// Prints sum of all query ranges. m is number of queries
// n is size of array a[].
void queryResults(int a[], int n, Query q[], int m)
{
	// Find block size √n
	block = (int)sqrt(n);

	// Sort all queries so that queries of same blocks
	// are arranged together.
	sort(q, q + m, compare);

	// Initialize current L, current R and current sum
	int currL = 0, currR = 0;
	int currSum = 0;

	// Traverse through all queries
	for (int i = 0; i < m; i++)
	{
		// L and R values of current range
		int L = q[i].L, R = q[i].R;

		// Remove extra elements of previous range. For
		// example if previous range is [0, 3] and current
		// range is [2, 5], then a[0] and a[1] are subtracted
		while (currL < L)
		{
			currSum -= a[currL];
			currL++;
		}

		// Add Elements of current Range
		while (currL > L)
		{
			currSum += a[currL - 1];
			currL--;
		}
		while (currR <= R)
		{
			currSum += a[currR];
			currR++;
		}

		// Remove elements of previous range.  For example
		// when previous range is [0, 10] and current range
		// is [3, 8], then a[9] and a[10] are subtracted
		while (currR > R + 1)
		{
			currSum -= a[currR - 1];
			currR--;
		}

		// Print sum of current range
		cout << "Sum of [" << L << ", " << R
		     << "] is "  << currSum << endl;
	}
}

////////////////////////////////////////////////////////////////////////////////

/*	Binary Index Tree (BITree)
	query and update operations in O(log n) time.	*/

vector<int> BITree;

// Constructs and returns a Binary Indexed Tree for given
// array of size n.
void constructBITree(vector<int> arr, int n)
{
	// Create and initialize BITree[] as 0
	BITree.resize(n + 1, 0);

// Store the actual values in BITree[] using update()
	for (int i = 0; i < n; i++)
		updateBIT(n, i, arr[i]);

	// Uncomment below lines to see contents of BITree[]
	//for (int i=1; i<=n; i++)
	//     cout << BITree[i] << " ";

}

// Returns sum of arr[0..index]. This function assumes
// that the array is preprocessed and partial sums of
// array elements are stored in BITree[].
int getSum(int index)
{
	int sum = 0; // Iniialize result

	// index in BITree[] is 1 more than the index in arr[]
	index = index + 1;

	// Traverse ancestors of BITree[index]
	while (index > 0)
	{
		// Add current element of BITree to sum
		sum += BITree[index];

		// Move index to parent node in getSum View
		index -= index & (-index);
	}
	return sum;
}

// Updates a node in Binary Index Tree (BITree) at given index
// in BITree. The given value 'val' is added to BITree[i] and
// all of its ancestors in tree.
void updateBIT(int n, int index, int val)
{
	// index in BITree[] is 1 more than the index in arr[]
	index = index + 1;

	// Traverse all ancestors and add 'val'
	while (index <= n)
	{
		// Add 'val' to current node of BI Tree
		BITree[index] += val;

		// Update index to that of parent in update View
		index += index & (-index);
	}
}

////////////////////////////////////////////////////////////////////////////////