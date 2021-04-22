// Created by Ashish Negi

/******** User-defined Function *******/

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

/* Returns prime bool vector using sieve of Eratosthenes algorithm */
vector<bool> SieveOfEratosthenes(int n=1000000)
{
    // Create a boolean array
    // "prime[0..n]" and initialize
    // all entries it as true.
    // A value in prime[i] will
    // finally be false if i is
    // Not a prime, else true.
	vector<bool> prime(n+1);
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
    for (long long int i=2; i<N ; i++)
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
        for (long long int j=0;
             j < (int)prime.size() &&
             i*prime[j] < N && prime[j] <= SPF[i];
             j++)
        {
            isprime[i*prime[j]]=false;

            // put smallest prime factor of i*prime[j]
            SPF[i*prime[j]] = prime[j] ;
        }
    }
}

/* Returns length of longest common subsequence for X[0..m-1], Y[0..n-1] */
int lcs(string X, string Y, int m, int n )
{
    vector<int> L(m+1, vector<int>(n+1));
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

	vector<int> LCSuff(m+1, vector<int>(n+1));
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
vector<vector<pair<int,int>>> adj;
vector<pair<int,pair<int,int>>> v;
vector<pair<int,int>> vals;

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
    vector<pair<int,int>> min_e(n,make_pair(-1,INF));
	min_e[0].second = 0;

    for (int i=0; i<n; ++i) {
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

vector<vector<pair<int,int>>> adj;
vector<pair<int,pair<int,int>>> v;
vector<pair<int,int>> vals;

/*	Add an edge (weighted, if isWeighted = true) from u to v in a graph	*/
void addEdge(int u,
	int v, bool isWeighted , int wt)
{
	if(isWeighted)
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

	for(int i=0; i<v.size(); i++)
	{
		addEdge(v[i].second.first, v[i].second.second,
			isWeighted, v[i].first);
	}

	return;

}
