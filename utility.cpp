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

void SieveOfEratosthenes(int n)
{
    // Create a boolean array
    // "prime[0..n]" and initialize
    // all entries it as true.
    // A value in prime[i] will
    // finally be false if i is
    // Not a prime, else true.
    bool prime[n + 1];
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
    for (int p = 2; p <= n; p++)
        if (prime[p])
            cout << p << " ";
}
