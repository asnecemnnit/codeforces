// Created by Ashish Negi

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
