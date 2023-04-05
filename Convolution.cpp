//!!! NTT
const int mod = 998244353,g = 3;
int root,root_inv;
void fft(vector<int> & a, bool invert) {
int n = a.size();
for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit; if (i < j) swap(a[i], a[j]);
}
for (int len = 2; len <= n; len <<= 1) {
    int wlen = invert ? root_inv : root;
    for (int i = len; i < sz(a); i <<= 1)
        wlen = (int) (1LL * wlen * wlen % mod);
    for (int i = 0; i < n; i += len) {
        int w = 1;
        for (int j = 0; j < len / 2; j++) {
            int u = a[i + j], v = (int) (1LL * a[i + j + len / 2] * w % mod);
            a[i + j] = u + v < mod ? u + v : u + v - mod;
            a[i + j + len / 2] = u - v >= 0 ? u - v : u - v + mod;
            w = (int) (1LL * w * wlen % mod);}}}
if (invert) {
    int n_rev = modi(n, mod); for (int &x: a)
    x = (int) (1LL * x * n_rev % mod);}}
vector<int> multiply(vector<int> a, vector<int> b) {
    int n = 1; 
    while (n < sz(a) + sz(b)) n <<= 1;
    a.resize(n); b.resize(n);
    root = power(g, (mod - 1) / n, mod);
    root_inv = modi(root, mod);
    fft(a, false); fft(b, false);
    for (int i = 0; i < n; i++) {
        a[i] = (1LL * a[i] * b[i]) %= mod;}
    fft(a, true); return a;}

//!!DC DP
ll n;
// kth level , find for [l,r], optimal break lies btw [optl,optr]
void f(ll k,ll l, ll r, ll optl=1, ll optr=n){
    if(r<l)return;
    ll mid= (l+r)/2;
    ll op=optl;
    for(ll br=optl;br<=min(optr,mid -1);br++){
        ll nval=dp[br][k-1]+ Cost(br+1,mid);  //only change this line
        if(nval<dp[mid][k]){
            dp[mid][k]=nval;
            op=br;}
    }
    f(k,l,mid-1,optl,op);f(k,mid+1,r,op,optr);}

//!!CHT
using T = ll; // a/b rounded down
ll fdiv(ll a, ll b) { return a/b-((a^b)<0&&a%b); }

bool _Q = 0;
struct Line {
    T a, b; mutable T lst;
    /// friend str ts(const Line& L) { return ts(vl{L.a,L.b,L.lst}); }
    T eval(T x) const { return a*x+b; }
    bool operator<(const Line&o)const{return _Q?lst<o.lst:a<o.a;}
    T last_gre(const Line& o) const { assert(a <= o.a);
        // greatest x s.t. a*x+b >= o.a*x+o.b
        return lst=(a==o.a?(b>=o.b?INF:-INF):fdiv(b-o.b,o.a-a));}
};

struct LineContainer: multiset<Line> {
bool isect(iterator it) { auto n_it = next(it);
    if (n_it == end()) return it->lst = INF, 0;
    return it->last_gre(*n_it) >= n_it->lst; }
void add(T a, T b) { // remove lines after
    auto it = insert({a,b,0}); while (isect(it)) erase(next(it));
    if (it == begin()) return;
    if (isect(--it)) erase(next(it)), isect(it);
    while (it != begin()) { // remove lines before
        --it; if (it->lst < next(it)->lst) break;
        erase(next(it)); isect(it); }
}
T qmax(T x) { assert(!empty());
    _Q = 1; T res = lower_bound({0,0,x})->eval(x); _Q = 0;
    return res; }
};

struct LCdeque : deque<Line> {
void addBack(Line L) { // assume nonempty
    while (1) {
        auto a = back(); pop_back(); a.lst = a.last_gre(L);
        if (size() && back().lst >= a.lst) continue;
        push_back(a); break;
    }
    L.lst = INF; push_back(L);
}
void addFront(Line L) {
    while (1) {
        if (!size()) { L.lst = INF; break; }
        if (((L.lst = L.last_gre(front()))) >= front().lst) pop_front();
        else break;
    }
    push_front(L);
}
void add(T a, T b) { // line goes to one end of deque
    if (!size() || a <= front().a) addFront({a,b,0});
    else assert(a >= back().a), addBack({a,b,0});
}
int ord = 1; // 1 = x's come in increasing order, -1 = decreasing order
T query(T x) {
    if(size() == 0)return INF;
    assert(ord);
    if (ord == 1) {
        while (front().lst < x) pop_front();
        return front().eval(x);
    } else {
        while(size()>1&&prev(prev(end()))->lst>=x)pop_back();
        return back().eval(x);
    }
}
};
//LCDeque() => set order of query , use add to add a line Amortized O(N)
//LineContainer() =>use add and qmax, uAmortized O(Nlog(N))

//!!FST
void FST(vll& a, bool inv) {
for (ll n = sz(a), step = 1; step < n; step *= 2) {
    for (ll i = 0; i < n; i += 2 * step) for(ll j=i;j<i+step;j++) {
        ll &u = a[j], &v = a[j + step]; tie(u, v) =
                                                    inv ? pll(v - u, u) : pll(v, u + v); // AND
        // inv ? pll(v, u - v) : pll(u + v, u); // OR /// include-line
        // pll(u + v, u - v);                   // XOR /// include-line
    }
}
// if (inv) for (ll& x : a) x /= sz(a); // XOR only /// include-line
}
vll conv(vll a, vll b) {
    FST(a, 0); FST(b, 0);
    for(ll i=0;i<sz(a);i++) a[i] *= b[i];
    FST(a, 1); return a;
}

//!!!SOS convolution
#define FOR(i,a,b) for(ll i=a;i<b;i++)
const int B = 20; // Every input vector must need to be of size 1<<B
// O(B * 2 ^ B)
vll SOSDP(vll f) {
    FOR(i,0,B) FOR(mask,0,1LL<<B)
            if ((mask & (1 << i)) != 0)f[mask] += f[mask ^ (1 << i)];
    //for supermasks make ((mask & (1 << i)) ==0
    return f;
}// you can change the operator from + to min/gcd to find min/gcd of all f[submasks]

//odd negation means if bit has odd set bits f[bit]=-f[bit]
//it can also be done by take oddnegation(f),SOS(f),oddnegation(f)
vll INVSOSDP(vll f) {
    FOR(i,0,B) FOR(mask,0,1LL<<B)
            if ((mask & (1 << i)) != 0)f[mask] -= f[mask ^ (1 << i)];
    return f;
}
// f*g(s)=sum_{s' $ s} {f(s')*g(s\s')}
// O(B * B * 2 ^ B)
vll subset_sum_convolution(vll f, vll g) {
    vvll fhat(B + 1, vll (1 << B, 0));
    vvll ghat(B + 1, vll (1 << B, 0));
    // Make fhat[][] = {0} and ghat[][] = {0}
    FOR(mask,0,1LL<<B) {
        fhat[__builtin_popcount(mask)][mask] = f[mask];
        ghat[__builtin_popcount(mask)][mask] = g[mask];
    }
    // Apply zeta transform on fhat[][] and ghat[][]
    for (int i = 0; i < B; i++) {
        fhat[i]=SOSDP(fhat[i]);ghat[i]= SOSDP(ghat[i]);
    }
    vvll h(B + 1, vll (1 << B, 0));
    // Do the convolution and store into h[][] = {0}
    FOR(mask,0,1LL<<B) FOR(i,0,B)FOR(j,0,i+1)
                h[i][mask] += 1LL * fhat[j][mask] * ghat[i - j][mask] ;
    // Apply inverse SOS dp on h[][]
    for (int i = 0; i < B; i++) {
        h[i]=INVSOSDP(h[i]);
    }
    vll fog(1 << B, 0);
    FOR(mask,0,1LL<<B)  fog[mask] = h[__builtin_popcount(mask)][mask];
    return fog;}

//!! ax+by =c
int egcd(int a,int b,int &x,int &y){
    if(b==0){
        x=1;y=0;return a;
    }
    int x1,y1;
    int d=egcd(b,a%b,x1,y1);
    x=y1;y=x1-y1*(a/b);
    return d;
}
// if(g!=1) does not exist,ainv mod b is x in egcd(a,b,x,y) % b (make +ve) 
// ax=b (mod n) is (ainv modn *b)%n (make sure to divide a,b,n by gcd first)
//ax+by=c
bool any_sol(int a,int b,int c,int &x,int &y,int &g) {
    if(a==0 && b==0) {
        g = x = y = 0;
        return true;}
    g = egcd(abs(a), abs(b), x, y);
    if (c % g != 0) {
        return false;
    }
    x *= (c / g);y *= (c / g);
    if (a < 0) x = -x;if (b < 0) y = -y;
    return true;
}
// x+ cnt*(b/g)  , y- cnt*(a/g) (so you can generate, count solutions in some range)

//!!Hacks
//* Find the smallest i in [a,b] that maximizes f(i), assuming that f(a) < f(i) <= f(b)$.
//* To reverse which of the sides allows non-strict inequalities, change the < marked with (A) to <=, and reverse the loop at (B).
//* To minimize f, change it to >, also at (B).
ll ternSearch(ll a, ll b) {
    assert(a <= b);
    while (b - a >= 5) {
        ll mid = (a + b) / 2;
        if (f(mid) < f(mid+1)) a = mid; // (A)
        else b = mid+1;}
    for(ll i=a+1;i<(b+1),i++) if (f(a) < f(i)) a = i; // (B)
    return a;
}

// checks if a*b will overflow
bool is(ll a,ll b){
    if(a==0 || b==0) return false;
    if((a>LLONG_MAX/b) or ( b>LLONG_MAX/a)) return true;
    return false;
}
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);}
    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);}
    size_t operator()(pair<uint64_t,uint64_t> x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x.first + FIXED_RANDOM)^(splitmix64(x.second + FIXED_RANDOM) >> 1);}

    size_t operator()(vector<ll> x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        uint64_t v=sz(x);
        for(int i=0;i<sz(x);i++){
            v^=(splitmix64(x[i]+FIXED_RANDOM)>>(i%4));
        }
        return v;
    }
};

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
//ll val=rng();vll v;shuffle(all(v),rng);
ll rand(ll x,ll y){return uniform_int_distribution<ll>(x,y)(rng);}