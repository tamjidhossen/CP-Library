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
#define debugGrid(x) cerr << #x <<" "; __print(x); cerr << endl;
#else
#define debug(x)
#define debugin(x)
#define debugGrid(x)
#endif

void _print(int t) {cerr << t;}void _print(string t) {cerr << t;}void _print(char t) {cerr << t;}
void _print(lld t) {cerr << t;}void _print(double t) {cerr << t;}void _print(ull t) {cerr << t;}
template <class T, class V> void _print(pair <T, V> p);
template <class T> void _print(vector <T> v);template <class T> void _print(set <T> v);
template <class T, class V> void _print(map <T, V> v);template <class T> void _print(multiset <T> v);
template <class T, class V> void _print(pair <T, V> p) {cerr << "{"; _print(p.first); cerr << ","; _print(p.second); cerr << "}";}
template <class T> void _print(vector <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void __print(vector<vector<T>> v) {cerr << endl; for(auto i:v) { for(auto j:i) {_print(j); cerr << " ";} cerr << endl;}}
template <class T> void _print(deque <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(set <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(multiset <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T, class V> void _print(map <T, V> v) {cerr << "[ "; for (auto i : v) {_print(i); cerr << " ";} cerr << "]";}


// -------------------------------------- Shorter Debug Code -------------------------------------- //
#ifndef ONLINE_JUDGE
#define debug(args...) cerr << "(" << #args << "):", _print(args);
#define debugGrid(x) cerr << #x <<" "; __print(x); cerr << endl;
#else
#define debug(args...)
#define debugGrid(x)
#endif

template<typename A, typename B> ostream& operator<<(ostream& os, const pair<A, B>&p) {return os<<'(' << p.first << ", " << p.second << ')';}
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream& os, const T_container& v) {os << '{'; string sep;for(const T& x: v) os << sep << x, sep = ", "; return os << '}';}
void _print() {cerr<<endl;}
template <typename Head, typename... Tail> void _print(Head H, Tail... T) {cerr << " " << H << ","; _print(T...);}
template <typename T> void __print(vector<vector<T>> v) {cerr << endl; for(auto i:v) { for(auto j:i) {cerr << j << " ";} cerr << endl;}}
// ------------------------------------------------------------------------------------------------- //



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
/*          Some Observations
    if  a mod k = x, one of the following holds:
        a mod 2k = x
        a mod 2k = x+k
    if a % k == b % k then (a-b)%k==0
*/


/*          Modular Arithmetic
    (a+b) mod m = ((a mod m)+(b mod m)) mod m

    (a−b) mod m = ((a mod m)−(b mod m)) mod m

    (a⋅b) mod m = ((a mod m)⋅(b mod m)) mod m

    (a/b) mod m = ((a mod m)⋅(b^-1 mod m)) mod m
*/


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


//------------------------- Number of Divisors ------------------------//
// N can be expressed as a product of its prime factors: N = p₁^e₁ × p₂^e₂ × ... × pₖ^eₖ
// The total number of divisors of N is given by:(e₁ + 1) × (e₂ + 1) × ... × (eₖ + 1)
int numberOfDivisors(int num) { // O(sqrt(n))
    if (num == 1) return 1;
    int total = 1;
    for (int i = 2; i * i <= num; i++) {
        if (num % i == 0) {
            int e = 0;
            do {
                e++;
                num /= i;
            } while (num % i == 0);
            total *= e + 1;
        }
    }
    if (num > 1) total *= 2;
    return total;
}
//--------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////









///////////////////////////////////////////////// Bit Manipulation /////////////////////////////////////////////////
// __builtin_popcountll(n); -> counts number of set bits
// check bit (n & (1LL << ith))
// n set bit number -> (1LL << n) - 1 
// set bit -> [ n | (1LL << ith) ];
// unset bit -> [ n & (~(1LL << ith)) ];
// to uppercase -> [ 'a' ^ 32 ];
// flip the kth bit -> X = X ^ (1LL << k);

// Use 1LL when shifting bits : (1LL << x)

/*
Some properties of bitwise operations:

    a|b = a⊕b + a&b
    a⊕(a&b) = (a|b)⊕b
    b⊕(a&b) = (a|b)⊕a
    (a&b)⊕(a|b) = a⊕b

Addition:

    a+b = a|b + a&b
    a+b = a⊕b + 2(a&b)

Subtraction:

    a-b = (a⊕(a&b))-((a|b)⊕a)
    a-b = ((a|b)⊕b)-((a|b)⊕a)
    a-b = (a⊕(a&b))-(b⊕(a&b))
    a-b = ((a|b)⊕b)-(b⊕(a&b))

*/
bool isPowerof2(int n){ return (!(n & (n - 1)) && n); }
int toggleBit(int n, int ith) { return n ^ (1 << ith);}
int numOfDigits(int n) { return floor(log10(n) + 1); } 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////











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

    // via inverse factorial
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











//////////////////////////////////////////////////////Graph///////////////////////////////////////////////////
/*                       |notes|
# for DFS in,
    (i)  acyclic graph(tree)     -> pass parameter(node, parent node)        
    (ii) cyclic graph            -> use visited array

# Single Source Shortest path [Dijkstra]
    (i)  weight =  1 [Can be done with BFS]
    (ii) weight >= 1 [Dijkstra]

# All Pair Shortest Path [Floyd Warshell]

# Tree Diameter (Max Distance Between two nodes)
    (i)   DFS from root and find the farthest node from root (say, x)
    (ii)  Now, from x DFS and find the farthest node (say, y)
    (iii) Distance between x and y is the diameter

# Path Between any two Nodes (let's say, x and y)
    (i)   Keep parent of ith node in vector while doing BFS from x to y
    (ii)  Backtrack from y to x by accessing last one's parent

# Bipartite Graph (a graph whose vertices can be divided into two disjoint and independent sets)
    (i)   two adjacent nodes can't have the same color
    (ii)  can have cycle
    (iii) bipartite graphs don't have odd length cycles, can have even lengh cycles; (means number of edges making a cycle)
*/

const int n = 10;
//------------------------------DFS--------------------------------//
vector<vector<int>>adj; //resize in solve() function
void DFS(int curr, int par) { // current and parent node
    for(auto i: adj[curr]) {
        if(i != par) {
            DFS(i, curr);
        }
    }
}
//-----------------------------------------------------------------//

//-------------------------------BFS-------------------------------//
vector<vector<int>>adj; // resize in solve() function
vector<bool>vis; // assign all false inside solve() function
void BFS(int initial) {
    queue<int>q;
    q.push(initial);
    while(!q.empty()) {
        int curr = q.front();
        cout<<curr<<" ";
        vis[curr] = true;
        for(auto i: adj[curr]) if(!vis[i]) q.push(i);
        q.pop();
    }
}
//-----------------------------------------------------------------//

//------------------------Connected Componnent----------------------------------//
vector<vector<int>>adj;
vector<bool>vis;
void DFS(int curr) { // current node
    vis[curr] = true;
    for(auto i: adj[curr]) {
        if(!vis[i]) {
            DFS(i);
        }
    }
}
void connectedComponnent() { 
    // adj.resize(n+1, vector<int>()); -> in solve()
    // vis.assign(n+1, false); -> in solve()
    int cc = 0; // number of connected component
    for(int i = 1; i <= n; i++) {
        if(!vis[i]) {
            cc++; DFS(i);
        }
    }
    cout << cc << endl;
}
//--------------------------------------------------------------------------------//



//------------------------Single Source Shortest Path with BFS(w = 1)----------------------------------//
vector<vector<int>>adj;
vector<int>dist; // from 1st to other nodes
void BFS(int initial) {
    queue<int>q;
    q.push(initial);
    dist[initial] = 0; 
    while(!q.empty()) {
        int curr = q.front();
        for(auto i: adj[curr]) {
            if(dist[i] > dist[curr]+1) { // here, curr is the parent of i
                dist[i] = dist[curr] + 1;      
                q.push(i);
            }      
        }
        q.pop();
    }
}
void SSSP_Unweighted_BFS(){ 
    // dist.assign(n+1, INT_MAX); -> in solve()
    // adj.resize(n+1, vector<int>()); -> in solve()
    int begNode = 1;
    BFS(begNode);
    int tar; cin>>tar;
    cout<<"Dist from "<<begNode<<" to "<<tar<<": "<<dist[tar]<<endl;
}
//-----------------------------------------------------------------------------------------------------//



//----------------------------------------------Tree Diameter------------------------------------------//
vector<vector<int>>adj;
int x, y;
int mx1 = 0, mx2 = 0;
void DFS1(int curr, int par, int lvl) {
    if(lvl > mx1) {
        mx1 = lvl; x = curr;
    }
    for(auto i: adj[curr]) {
        if(i != par) {
            DFS1(i, curr, lvl+1);
        }
    }
}
void DFS2(int curr, int par, int lvl) {
    if(lvl > mx2) {
        mx2 = lvl; y = curr;
    }
    for(auto i: adj[curr]) {
        if(i != par) {
            DFS2(i, curr, lvl+1);
        }
    }
}
void TreeDiameter() {
    // adj.resize(n+1, vector<int>()); -> in solve()
    DFS1(1, -1, 0);
    DFS2(x, -1, 0);
    cout<<"Distance between "<<x<<" and "<<y<<" is "<<mx2<<endl;
}
//-----------------------------------------------------------------------------------------------------//



//-----------------------------------Path Between Two Nodes--------------------------------------------//
vector<vector<int>>adj;
vector<int>parent;
vector<bool>vis;
int x, y;
void BFS(int initial) {
    queue<int>q;
    q.push(initial);
    while(!q.empty()) {
        int curr = q.front();
        vis[curr] = true;
        if(curr == y) return;
        for(auto i: adj[curr]) {
            if(!vis[i]) {
                q.push(i);
                parent[i] = curr;
            }
        }
        q.pop();
    }
}
void backtrack(vector<int>&v, int last) {
    while(last != -1) {
        v.pb(last);
        last = parent[last];
    }
    reverse(v.begin(), v.end());
}
void Path_Between_Two_Nodes(int x, int y) {
    // adj.resize(n+1, vector<int>()); -> in solve()
    BFS(x);
    vector<int>ans;
    backtrack(ans, y);
    for(int i = 0; i < size(ans); i++) cout<<ans[i]<<" ";
    cout<<endl;
} 
//-----------------------------------------------------------------------------------------------------//




//-------------------------------------Bipartite Graph-------------------------------------------------//
vector<vector<int>>adj;
vector<int>color;
bool dfs(int curr, int col) {
    color[curr] = col;
    for(auto i: adj[curr]) {
        if(color[i] == -1) {
            if(!dfs(i, col^1)) { // if false for next vertex
                return false;
            }
        } else {
            if(color[curr] == color[i]) return false;
        }
    }
    return true;
}
void Bipartite_Graph_Coloring() {
    // adj.resize(n+1, vector<int>()); -> in solve
    // color.assign(n+1, -1);-> in solve
    if(dfs(1, 0)) cout<<"Bipartite\n";
    else cout<<"Not Bipartite\n";
}




//-------------------------------------DSU-----------------------------------------//

vector<int>parent, childs;
struct {
    void init(int n) {
        parent.resize(n+1);
        childs.resize(n+1);
        for(int i = 0; i <= n; i++) {
            parent[i] = i;
            childs[i] = 1;
        }
    }
    int findRoot(int v) {
        if(v == parent[v]) return v;
        //return findRoot(parent[v]);
        //optimized by, every node pointing to SuperParent
        return parent[v] = findRoot(parent[v]);
    }
    bool uniteSets(int u, int v) {
        u = findRoot(u);
        v = findRoot(v);
        if(u != v) {
            if(childs[u] < childs[v]) swap(childs[u], childs[v]);
            parent[v] = u;
            childs[u] += childs[v];
            return true;
        }
        return false;
    }
} DSU;
//-------------------------------------------------------------------------------------------------------------//


////////////////////////////////////////////////////////////////////////////////////////////////////////////////










//////////////////////////////////////// String Algorithms //////////////////////////////////////////////////



//----------------------------------------------String Hashing-------------------------------------------------//
/* notes: 
hash(s) = s[0] + (s[1] * p) + (s[i] * p^2) + ... + (s[n - 1] * p^n-1)
hash(s[i...j]) = (hash(s[0...j]) - hash(0...i - 1)) / (p^i);
1. used double hashing.
2. base = 31 and mod = 1e9+7 works if string's value = s[i] - 'a' + 1; (1 - 26)
3. call precal_pow() from main()
4. if getting WA even after 2 Hashes:
    i. use unpopular primes (!1e9 + 7)
    ii. use 3 hashes (might give TLE as % is very heavy operation)
    iii. use other algo's
5. some base = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499;
6. use int instead of long long and cast to ll with 1LL when doing any multiplication
7. using Long Long Int could give TLE, because ll takes 2x more memory then int and some calculations (multiplication, modulo) are very expensive.
*/
int N = 1e6 + 9; 
const int MOD1 = 127657753, MOD2 = 987654319;
const int base1 = 137, base2 = 277;
vector<pair<int, int>>pw(N), inv_pw(N);

int BE(int x, int y, int mod){
    int res = 1;
    x %= mod;
    if(x < 0) x += mod;
    while(y > 0){
        if(y & 1) res = (1LL * res * x) % mod; // MOD
        y >>= 1; // -> y /= 2;
        x = (1LL * x * x) % mod; // MOD
    } // MOD for larger numbers
    return res;
}

void precal_pow() {
    pw[0] = {1, 1};
    inv_pw[0] = {1, 1};
    int inv_base1 = BE(base1, MOD1 - 2, MOD1);
    int inv_base2 = BE(base2, MOD2 - 2, MOD2);
    for(int i = 1; i < N; i++) {
        pw[i].first = (1LL * pw[i - 1].first * base1) % MOD1;
        pw[i].second = (1LL * pw[i - 1].second * base2) % MOD2;
        inv_pw[i].first = (1LL * inv_pw[i - 1].first * inv_base1) % MOD1;
        inv_pw[i].second = (1LL * inv_pw[i - 1].second * inv_base2) % MOD2;
    }
}

struct Hashing {
    int n;
    string s; // 0 - indexed
    vector<pair<int, int>> hs; // 1 - indexed
    Hashing() {}
    Hashing(string _s) {
        n = size(_s);
        s = _s;
        hs.emplace_back(0, 0);
        for(int i = 0; i < n; i++) {
            pair<int, int> p;
            p.first = (hs[i].first + (1LL * pw[i].first * s[i]) % MOD1) % MOD1;
            p.second = (hs[i].second + (1LL * pw[i].second * s[i]) % MOD2) % MOD2;
            hs.push_back(p);
        }
    }
    pair<int, int> get_hash(int l, int r) { // 1 - indexed
        assert(1 <= l && l <= r && r <= n); // will cause RE
        pair<int, int> ans;
        ans.first = ((hs[r].first - hs[l - 1].first + MOD1) * 1LL * inv_pw[l - 1].first) % MOD1;
        ans.second = ((hs[r].second - hs[l - 1].second + MOD2) * 1LL * inv_pw[l - 1].second) % MOD2;
        return ans;
    }
    pair<int, int> get_hash() { // gets full string hash
        return get_hash(1, n);
    }
};
//-------------------------------------------------------------------------------------------------------------//



//--------------------------------------------------Reverse Hashing & palindrom -----------------------------------------------------//
int N = 2e6 + 9; 
const int MOD1 = 127657753, MOD2 = 987654319;
const int base1 = 137, base2 = 277;
vector<pair<int, int>>pw(N), inv_pw(N);
 
int BE(int x, int y, int mod){
    int res = 1;
    x %= mod;
    if(x < 0) x += mod;
    while(y > 0){
        if(y & 1) res = (1LL * res * x) % mod; // MOD
        y >>= 1; // -> y /= 2;
        x = (1LL * x * x) % mod; // MOD
    } // MOD for larger numbers
    return res;
}
 
void precal_pow() {
    pw[0] = {1, 1};
    inv_pw[0] = {1, 1};
    int inv_base1 = BE(base1, MOD1 - 2, MOD1);
    int inv_base2 = BE(base2, MOD2 - 2, MOD2);
    for(int i = 1; i < N; i++) {
        pw[i].first = (1LL * pw[i - 1].first * base1) % MOD1;
        pw[i].second = (1LL * pw[i - 1].second * base2) % MOD2;
        inv_pw[i].first = (1LL * inv_pw[i - 1].first * inv_base1) % MOD1;
        inv_pw[i].second = (1LL * inv_pw[i - 1].second * inv_base2) % MOD2;
    }
}
 
struct Hashing {
    int n;
    string s; // 0 - indexed
    vector<pair<int, int>> hs, rhs; // 1 - indexed
    Hashing() {}
    Hashing(string &_s) {
        n = size(_s);
        s = _s;
        hs.emplace_back(0, 0);
        rhs.emplace_back(0, 0);
        for(int i = 0; i < n; i++) {
            pair<int, int> p, rp;
            p.first = (hs[i].first + (1LL * pw[i].first * s[i]) % MOD1) % MOD1;
            p.second = (hs[i].second + (1LL * pw[i].second * s[i]) % MOD2) % MOD2;
            hs.push_back(p);
            
            rp.first = (rhs[i].first + (1LL * pw[i].first * s[n - i - 1]) % MOD1) % MOD1;
            rp.second = (rhs[i].second + (1LL * pw[i].second * s[n - i - 1]) % MOD2) % MOD2;
            rhs.push_back(rp);
        }
    }
    pair<int, int> get_hash(int l, int r) { // 1 - indexed
        assert(1 <= l && l <= r && r <= n); // will cause RE
        pair<int, int> ans;
        ans.first = ((hs[r].first - hs[l - 1].first + MOD1) * 1LL * inv_pw[l - 1].first) % MOD1;
        ans.second = ((hs[r].second - hs[l - 1].second + MOD2) * 1LL * inv_pw[l - 1].second) % MOD2;
        return ans;
    }
    pair<int, int> get_rev_hash(int _l, int _r) { // 1 - indexed
        int l = n - _r + 1; 
        int r = n - _l + 1;
        assert(1 <= l && l <= r && r <= n); // will cause RE
        pair<int, int> ans;
        ans.first = ((rhs[r].first - rhs[l - 1].first + MOD1) * 1LL * inv_pw[l - 1].first) % MOD1;
        ans.second = ((rhs[r].second - rhs[l - 1].second + MOD2) * 1LL * inv_pw[l - 1].second) % MOD2;
        return ans;
    }
    pair<int, int> get_hash() { // gets full string hash
        return get_hash(1, n);
    }
    bool is_palindrome(int l, int r) {
        return (get_hash(l, r) == get_rev_hash(l, r));
    }
};
//------------------------------------------------------------------------------------------------------------------//



//----------------------------------------------Manachers Algorithm-------------------------------------------------//
/* notes: 
1. modifies string by adding '#' bef&after every char;
2. p[i] = size of palindrome that exist in original string; here i is the center of that palindrome of modified string;
3. palindrome at ith(modified string's) position is: s.substr((index/2) - (pal_size/2), pal_size);
*/
struct Manachers { // O(n)
    vector<int> p;
    string s;
    Manachers(string &_s) {
        s = "#";  // build new string
        for(int i = 0; i < size(_s); i++) {
            s.push_back(_s[i]);  s.push_back('#');
        }
        int n = size(s);
        p.assign(n, 0);
        int l = 0, r = 0; //curr largest pal range
        for(int i = 0; i < n; i++) {
            int mirror = l + (r - i);
            if (i < r)  p[i] = min(r - i, p[mirror]);
            while(i + p[i] + 1 < n && i - p[i] - 1 >= 0 && s[i + p[i] + 1] == s[i - p[i] - 1]) {
                p[i]++;
            }
            if(i + p[i] > r) {  // if new largest pal found
                l = i - p[i]; r = i + p[i];
            }
        }
    }
    bool is_palindrome(int l, int r) { // O(1)
        l = 2 * l + 1; r = 2 * r + 1;
        int center = (l + r) / 2; // of modified string
        return p[center] >= (r - l + 1) / 2;
    }
};
//------------------------------------------------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////










/////////////////////////////////////////////////// SORTING ////////////////////////////////////////////////////

void merge(vector<int>&v, int l, int mid, int r) {
    int p = l, q = mid + 1, k = 0;
    vector<int>temp(r - l + 1);
    for(int i = l; i <= r; i++) {
        if(p > mid) temp[k++] = v[q++];
        else if(q > r) temp[k++] = v[p++];
        else if(v[p] < v[q]) temp[k++] = v[p++];
        else {
            temp[k++] = v[q++];
            // inv += mid - p + 1; // inv += q - p - (q - mid - 1) -> for inverse calculation
        }
    }
    for(int i = 0; i < k; i++) v[l++] = temp[i];
}

void mergeSort(vector<int>&v, int l, int r) {
    if(l < r) {
        int mid = (l + r) / 2;
        mergeSort(v, l, mid);
        mergeSort(v, mid+1, r);
        merge(v, l, mid, r);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void solve()
{
    int n = 10;
    vector<int>v(n);
    vector<int>newV(v.begin(), v.begin() + n); // copies first 10 elements of v to newV

    string name = "HelloWorld";
    set<char> my_name(name.begin(), name.end()); // stores all elements of string

    deque<int> dq; // it's just a duble ended vector


    /////////upper_bound, lower_bound functions///////////

    auto lower = lower_bound(all(v), n); //returns iterator
    // index of first occurrence of the value in the range
    int upper = upper_bound(all(v), n) - v.begin();
    // index of the first occurrence of a value that is greater than the searched value
    //If not present, than returns the next index

    ///////////////////////////////////////////////////////
    
    ////////Function inside main()///////// -> Lambda Function
    vector<int> v = {1, 2, 3};
    auto func = [&](vector<int>v) -> int { //[&] means all memory of main() is accessable
        int sum = 0;                       // () parameter can be empty
        for(auto i: v) sum += i;
        return sum;
    };
    cout<<func(v)<<endl; // returns 6;
    //////////////////////////////////////

    /////////////BITSET///////////////////
    // faster then vector<bool>
    // can be used to find if similar value in two different container (using &)
    bitset<8> b1;              // 00000000               // O(1)
    bitset<8> b2("10101010");  // From string            // O(N)
    b1[1] = 1;                 // Set 2nd bit            // O(1)
    b1.set(3);                 // Set 4th bit            // O(1)
    b1.reset(1);               // Clear 2nd bit          // O(1)
    b1.flip(0);                // Flip 1st bit           // O(1)

    int x = b1.count();        // Number of 1s           // O(N)
    bool any = b1.any();       // True if any bit is 1   // O(N)
    bool none = b1.none();     // True if all 0          // O(N)
    bool test = b1.test(2);    // True if 3rd bit is 1   // O(1)

    bitset<8> b3 = b1 & b2;    // AND                    // O(N)
    bitset<8> b4 = b1 | b2;    // OR                     // O(N)
    bitset<8> b5 = b1 ^ b2;    // XOR                    // O(N)

    string s = b1.to_string(); // Convert to string      // O(N)
    unsigned long val = b1.to_ulong(); // To ulong       // O(N)
    /////////////////////////////////////

    /////////sort vector of pair ////////
    vector<pair<int,int>> a;
    a.push_back({1,2});
    a.push_back({2,3});
    a.push_back({1,4});
    a.push_back({2,1});

    ////////////////////////Custome Sort////////////////////////
    sort(a.begin(),a.end(),[&](pair<int,int>x, pair<int,int>y){
        if(x.first==y.first){
            return (x.second<y.second);
        }
        return (x.first<y.first);
    });
    //////////////////////////////////////////////////////////

    ///////////isPerfectSquare///////////
    auto isPrefectSquare = [&](int value) -> bool { // lambda function
        return (int)sqrtl(value) * (int)sqrtl(value) == value;
    };
    ////////////////////////////////////

    ///////////Binary Search////////////
    auto check = [&](int mid) -> bool { 
        
        return false;
    };
    int lo = 0, hi = n - 1, ans, target;
    while(lo <= hi){
        int mid = lo + (hi - lo) / 2;
        if(check(mid)){
            ans = mid;
            hi = mid - 1;
        } else lo = mid + 1;
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

/* TIPS
    -> do as less division as possible
    -> float (least precise) < double < long double (most precise)
    -> check for integer overflow
    -> instead of a == b, use fabs(a - b) < epsilon(1e-9) in comparing double/float
    -> using a fixed number of iterations in BS (e.g., 100) is often more reliable than epsilon-based (e.g., absolute or    relative error of 10^6) termination
    -> asigning lo and hi as constant(e.g, 1e9) numbers may result in wrong answer and we should calculate bounds based on problems constraints for floating based precesion calculating BS
    -> check for out of bound in arrays (Mostly Results in RTE)
    -> any posibility of negative index?
    -> outer loop and inner loop shouldn't have same iterator variables
    -> To pass tricky test cases we can sometime shuffle the array
        mt19937 rng(chrono::steady_clock::now().time_since_epoch().count()); //generates seed
        shuffle(nums.begin(), nums.end(), rng); // O(n)
    -> sin, cos, asin(sin ^ -1), acos ... ... returns radian; multiply (pi/180) as a parameter to convert to degree;
    -> Heavy - Light trick, if string length > sqrt(max_string_length) : use hashing else : use trie
*/
/*  Send debug log to file
    ofstream fcerr("log.txt"); // paste before debug template
    // replace all cerr with fcerr
    fcerr.close(); // close the file before return 0;
*/