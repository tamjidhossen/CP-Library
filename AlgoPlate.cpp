#include<bits/stdc++.h>
using namespace std;
#define endl '\n'
#define pb push_back
#define int long long
#define lld long double
#define ull unsigned long long
#define size(x) (int) x.size()
#define all(x) (x).begin(), (x).end() 

#ifndef ONLINE_JUDGE
#define debug(x) cerr << #x <<" "; _print(x); cerr << endl;
#define debugin(x) cerr << #x <<" "; _print(x); cerr << "; ";
#else
#define debug(x)
#define debugin(x)
#endif

void _print(int t) {cerr << t;}void _print(string t) {cerr << t;}void _print(char t) {cerr << t;}
void _print(lld t) {cerr << t;}void _print(double t) {cerr << t;}void _print(ull t) {cerr << t;}
template <class T, class V> void _print(pair <T, V> p);
template <class T> void _print(vector <T> v);template <class T> void _print(set <T> v);
template <class T, class V> void _print(map <T, V> v);template <class T> void _print(multiset <T> v);
template <class T, class V> void _print(pair <T, V> p) {cerr << "{"; _print(p.first); cerr << ","; _print(p.second); cerr << "}";}
template <class T> void _print(vector <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(deque <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(set <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(multiset <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T, class V> void _print(map <T, V> v) {cerr << "[ "; for (auto i : v) {_print(i); cerr << " ";} cerr << "]";}





///////////////////////////////////////////////// int128 /////////////////////////////////////////////////
/*

// int128 bit for numbers larger than 1e18. Will support numbers till 1e36
// Typedef to ell -> extra long long
typedef __int128 ell;

// For printing
std::ostream&
operator<<( std::ostream& dest, __int128_t value ) {
	std::ostream::sentry s(dest);
	if (s) {
		__uint128_t tmp = value < 0 ? -value : value; char buffer[ 128 ];
		char* d = std::end( buffer );
		do {	-- d; *d = "0123456789"[ tmp % 10 ]; tmp /= 10;} while ( tmp != 0 );
		if ( value < 0 ) {-- d; *d = '-';}
		int len = std::end( buffer ) - d;
		if ( dest.rdbuf()->sputn( d, len ) != len ) {dest.setstate( std::ios_base::badbit );}
	}
	return dest;
}

// For reading _int128 to_read = read()
__int128 read() {
	__int128 x = 0, f = 1;
	char ch = getchar();
	while (ch < '0' || ch > '9') {if (ch == '-') f = -1; ch = getchar();}
	while (ch >= '0' && ch <= '9') {x = x * 10 + ch - '0'; ch = getchar();}
	return x * f;
}

// For debugging
void _print(ell t) {cerr << t;}

*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////







#define MOD 1000000007
///////////////////////////////////////////////// Number Theory /////////////////////////////////////////////////
/*          Modular Arithmetic
    (a+b) mod m = ((a mod m)+(b mod m)) mod m

    (a−b) mod m = ((a mod m)−(b mod m)) mod m

    (a⋅b) mod m = ((a mod m)⋅(b mod m)) mod m

    (a/b) mod m = ((a mod m)⋅(b^-1 mod m)) mod m
*/
// __builtin_popcountll(n); -> counts number of set bits
// set bit -> [ n | (1 << ith) ];
// unset bit -> [ n & (~(1 << ith)) ];
// to uppercase -> [ 'a' ^ 32 ];
bool isPowerof2(int n){ return n & (n-1); }
int toggleBit(int n, int ith) { return n ^ (1 << ith);}
int numOfDigits(int n) { return floor(log10(n) + 1); } 
//--------------------Binary Exponentiation----------------------//
/*
    # Inverse : (n^-1)%m = (n^m-2)%m
        Ex: Inverse of x would be,
        BinaryExponentiation(x, MOD-2);
    # For normal constrain just divide by x(number)
*/ 
int BinaryExponentiation(int x, int y){ 
    int res = 1;
    while(y > 0){
        if(y & 1) res *= x; // MOD
        y >>= 1; // -> y /= 2;
        x *= x; // MOD
    } // MOD for larger numbers
    return res;
}
//----------------------------------------------------------------//
//----------------------------Sieve------------------------------//
const int N = 1e7+7;
bitset<N> marked;
vector<int>primes;
void sieve(int n){
    marked.reset(); //sets every bit to 0
    // vector<bool>marked(N, false);
    marked[1] = true; // 1 not prime
    for(int i = 2; i*i <= n; i++){
        if(marked[i]) continue;
        for(int j = i*i; j <= n; j += i){
            marked[j] = true; // marked[i] == 0 means prime
        }
    }
    for(int i = 2; i < N; i++) if(!marked[i]) primes.pb(i);
}
//----------------------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////










///////////////////////////////////////////////// Combinatorics /////////////////////////////////////////////////
/*
    # All non empty substring of length n: (n*(n+1))/2;
    # All non empty subsequence of length n (doesn't has to be consequcative): 2^n - 1 // either choose one or not

*/
//--------------------Pre Calculate Factorial----------------------//
int Factorial[N+1];
void factorialPreCal() {
    int ithFact = 1;
    for(int i = 1; i <= N; i++) {
        ithFact *= i;
        // ithFact %= MOD;
        Factorial[i] = ithFact;
    }
}
//----------------------------------------------------------------//

//------------------------------nCr-------------------------------//
// nCr = [n! / {r! * (n-r)!}]
int nCr(int n, int r) { 
    if(r >= n) return 0LL;
    int numerator = Factorial[n];
    int denominator = Factorial[r] * Factorial[n-r];
    int ans = numerator / denominator;
    return ans;

    // int numerator = Factorial[n] % MOD;
    // int denominator = (Factorial[r] % MOD * Factorial[n-r] % MOD) % MOD;
    // denominator = BinaryExponentiation(denominator, MOD-2); // for inverse
    // int ans = (numerator % MOD * denominator % MOD) % MOD;
    // return ans;
}
//-----------------------------------------------------------------//

//------------------------------nPr-------------------------------//
// nPr = [n! / (n-r)!]
int nPr(int n, int r) { 
    if(r >= n) return 0LL;
    int numerator = Factorial[n];
    int denominator = Factorial[n-r];
    int ans = numerator / denominator;
    return ans;

    // int numerator = Factorial[n] % MOD;
    // int denominator = Factorial[n-r] % MOD;
    // denominator = BinaryExponentiation(denominator, MOD-2); // for inverse
    // int ans = (numerator % MOD * denominator % MOD) % MOD;
    // return ans;
}
//-----------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////











void solve()
{
    int n = 10;
    vector<int>v(n);

    string name = "HelloWorld";
    set<char> my_name(name.begin(), name.end()); // stores all elements of string

    deque<int> dq; // it's just a duble ended vector
    
    ////////Function inside main()/////////
    vector<int> v = {1, 2, 3};
    auto func = [&](vector<int>v) -> int { //[&] means all memory of main() is accessable
        int sum = 0;                       // () parameter can be empty
        for(auto i: v) sum += i;
        return sum;
    };
    cout<<func(v)<<endl; // returns 6;
    //////////////////////////////////////

    /////////////BITSET///////////////////
    bitset<10>bits(11);
    cout<<bits<<endl; // 0000000101 <- accessable as a string
    cout<<bits[3]<<endl; // returns 4th bit from right
    cout<<bits.to_ullong()<<endl; // converts to unsigned ll
    cout<< bits.to_string()<<endl
        << bits.to_string('*')<<endl
        << bits.to_string('O', 'X')<<endl;
    /////////////////////////////////////

    /////////sort vector of pair ////////
    vector<pair<int,int>> a;
    a.push_back({1,2});
    a.push_back({2,3});
    a.push_back({1,4});
    a.push_back({2,1});

    sort(a.begin(),a.end(),[&](pair<int,int>x, pair<int,int>y){
        if(x.first==y.first){
            return (x.second<y.second);
        }
        return (x.first<y.first);
    });
    /////////////////////////////////////

    ///////////isPerfectSquare///////////
    auto isPrefectSquare = [&](int value) -> bool { // lambda function
        return (int)sqrtl(value) * (int)sqrtl(value) == value;
    };
    ////////////////////////////////////

    ///////////Binary Search////////////
    int lo = 0, hi = n-1, ans, target;
    while(lo <= hi){
        int mid = lo + (hi - lo)/2;
        if(v[mid] == target){
            ans = mid;
            break;
        }
        else if(v[mid] < target) lo = mid + 1;
        else if(v[mid] > target) hi = mid - 1;
    }
    ///////////////////////////////////
}

int32_t main()
{
    ios_base::sync_with_stdio(false); cin.tie(NULL);

    int t = 1; 
    cin>>t; 
    while(t--) solve();

    return 0;
}