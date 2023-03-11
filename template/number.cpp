const long long mod = 998244353;
template<long long mod>
class MINT {
public:
    long long x;
    MINT(const int & x) {this->x = x; }
    MINT(const long long & x) {this->x = x; }
    MINT() {this->x = 0;}
    constexpr long long val() {return x;}
    void norm() { x %= mod; x+= mod; x %= mod;}
    const long long operator()() const { return x; }
    MINT& operator += (const MINT& other) {x += other.x; if (x >= mod) x -= mod; return *this;}
    MINT& operator -= (const MINT& other) {x -= other.x; if (x < 0) x += mod; return *this;}
    MINT& operator *= (const MINT& other) {x *= other.x; x %= mod; return *this;}
    MINT& operator /= (const MINT& other) {x *= other.inv()(); x %= mod; return *this;}
    MINT inv() const {
        return power(mod - 2);
    }
    template<class U> MINT power(U n) const {
        MINT a = x, res = 1;
        while (n > 0) {
            if (n & 1) res *= a;
            a *= a;
            n >>= 1;
        }
        return res;
    }
    template<class T> MINT& operator += (const T& other) {return *this += MINT(other);}
    template<class T> MINT& operator -= (const T& other) {return *this -= MINT(other);}
    MINT& operator ++ () {return *this += 1;}
    MINT& operator -- () {return *this -= 1;}    
};
typedef MINT<mod> mint;
// typedef double mint;

ostream& operator << (ostream& out, mint x) {
    out << x();
    return out;
}

istream& operator >> (istream& in, mint &x) {
    in >> x.x;
    return in;
}

#ifdef DEBUG
void _debug(string s, mint x) {
    cerr << s << ": " << x() << endl;
}
#endif

template<class T, class U> T operator+ (const T &x, const U &y) {T ret = T(x); ret += T(y); return ret;}
template<class T, class U> T operator- (const T &x, const U &y) {T ret = T(x); ret -= T(y); return ret;}
template<class T, class U> T operator* (const T &x, const U &y) {T ret = T(x); ret *= T(y); return ret;}
template<class T, class U> T operator/ (const T &x, const U &y) {T ret = T(x); ret /= T(y); return ret;}
 
mint jc[1100030];
mint jc_1[1100030];
 
mint c[255][255];
void init() {
    jc[0] = jc_1[0] = 1;
    for (int i = 1; i < 1100030; ++i) {
        jc[i] = jc[i - 1] * i;
        jc_1[i] = mint(jc[i]).inv()();
    }
    memset(c, 0, sizeof(c));
    c[0][0] = 1;
    for (int i = 1; i < 255; i++) {
        c[i][i] = c[i][0] = 1;
        for (int j = 1; j < i; j++) {
            c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
        }
    } 
}

mint get_c(long long m, long long n) {
    if (n < 0 || n > m) return 0;
    if (m < 255) return c[m][n];
    return jc[m] * jc_1[n] * jc_1[m - n];
}
