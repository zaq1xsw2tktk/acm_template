long long mod = 998244353;
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