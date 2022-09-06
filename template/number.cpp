long long dp[255][255];
long long jc[300030];
long long mod = 998244353;
long long c[255][255];
long long power(long long x, long long n) {
    if (n % 2 == 1) {
        long long r = power(x, n / 2);
        return r * r % mod * x % mod;
    }
    if (n == 0) {return 1;}
    long long r = power(x, n / 2);
    return r * r % mod;
}
long long li(long long x) {
    return power(x, mod - 2);
}
void init() {
    jc[0] = 1;
    for (int i = 1; i < 300030; ++i) {
        jc[i] = jc[i - 1] * i % mod;
    }
    memset(c, 0, sizeof(c));
    c[0][0] = 1;
    for (int i = 1; i < 255; i++) {
        c[i][i] = c[i][0] = 1;
        for (int j = 1; j < i; j++) {
            c[i][j] = c[i - 1][j - 1] + c[i - 1][j];
            c[i][j] %= mod;
        }
    } 
}
long long get_c(long long m, long long n) {
    if (n < 0 || n > m) return 0;
    if (m < 255) return c[m][n];
    return jc[m] * li(jc[n]) % mod * li(jc[m - n]) % mod;
}